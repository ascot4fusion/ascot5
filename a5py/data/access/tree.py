"""Module for building a tree structure for accessing data.

The main purpose of the tree is to help the user to navigate and manipulate the
data using the public attributes of each component. Here's what references
different components have: ::

    Tree (root)
    ├── InputCategory (input category)
    │   └── Leaf (input variant)
    |       └── Public attributes and methods to access the actual data
    └ OutputLeaf (run variant)
      └── Public attributes and methods to access the actual data

What is not shown, and why this structure is not a proper tree, is that the run
variants contain public references to the input variants that were used.

The terms leafs and nodes are used within the context of the tree to emphasize
that the tree manages only the metadata and the tree structure. Internally,
everything is managed by `TreeManager` instance which contains all `Node` and
`Leaf` instances. Likewise, `Node` and `Leaf` instances have a reference to
the tree manager. Since the user cannot navigate from leaf to to it's parent
node, the leaf contains no reference to the parent. Instead, the leaf must use
the tree manager to trigger updates on the tree.
"""
from __future__ import annotations

from typing import Optional, TYPE_CHECKING

from a5py.exceptions import AscotDataException, AscotMeltdownError
from .hdf5io import TreeFileManager, DataAccess
from .leaves import Leaf, Status, InputVariant, OutputVariant

if TYPE_CHECKING:
    from .nodes import ImmutableNode


ROOT = "root"
"""Name of the root/output node."""


# Manager is allowed to access private parts of nodes and leafs
# pylint: disable=protected-access
class TreeManager(TreeFileManager):
    """Manages contents of the tree.

    This manager operates in the background and is responsible for keeping
    the contents of the tree consistent. The manager has references to all
    nodes and their leafs. To avoid cyclic calls, this manager is(only) called
    when:

    - The tree builds its nodes during initialization (initialize this
      instance).
    - New leaf is added to the tree (call enter_leaf).
    - Leaf destroys itself (call remove_leaf).

      - When node destroys a leaf, it should call the leaf's destroy method and
        not do else.

    - Leaf activates itself (call activate_leaf).
    - The note on a leaf is updated (call note_changed) since this may update
      the tag that is used to reference the leaf.
    """

    def __init__(
            self, hdf5file: Optional[tuple[str, bool]]=None,
            **nodes: ImmutableNode,
            ) -> None:
        """
        Parameters
        ----------
        hdf5file : tuple[str, bool], *optional*
            Name of the HDF5 file and if the file exists.

            If this value is given, the data will be stored and read from the
            file (not memory).
        **nodes : :class:`.ImmutableNode`
            Nodes belonging to this tree (including the root node).
        """
        self.inputs: list[InputVariant] = []
        """Leafs that belong to the input nodes."""

        self.outputs: list[OutputVariant] = []
        """Leafs that belong to the root node."""

        self.usedin: dict[InputVariant, list[OutputVariant]] = {}
        """Maps input leafs to the outputs leafs that have used them."""

        self.uses: dict[OutputVariant, list[InputVariant]] = {}
        """Maps output leaf to the inputs it has used."""

        self.nodes: dict[str, ImmutableNode] = nodes
        """All nodes and their names that make up the tree."""

        for node in nodes.values():
            with node._modify_attributes():
                node._treemanager = self

        self.filename = ""
        if hdf5file is not None:
            filename, file_exists = hdf5file
            super().__init__(
                root=ROOT, filename=filename, file_exists=file_exists,
                input_categories=[n for n in nodes if n != ROOT]
                )
            self._init_from_hdf5([n for n in nodes if n != ROOT])

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(inputs={self.inputs}, "
            f"outputs={self.outputs})>"
        )

    def _init_from_hdf5(self, input_categories: list[str]) -> None:
        """Initializes the tree structure from the HDF5 file.

        Parameters
        ----------
        input_categories : list[str]
            All input categories present in the tree.
        """
        for category in input_categories:
            active, names = self.read_node(category)
            for name in names:
                date, note, _ = self.read_variant(name, category)
                leaf = Leaf.create_leaf(
                    date=date, note=note, variant=name.split("_")[0]
                    )
                if not isinstance(leaf, InputVariant):
                    raise AscotMeltdownError("Input variant expected: " + name)
                leaf._name = name
                self.enter_leaf(leaf, category, save=False, new_entry=False)
                leaf._file = self.access_data(name, category)

            if active:
                self.get_leaf(active).activate()

        active, names = self.read_node(ROOT)
        for name in names:
            date, note, inputs = self.read_variant(name, ROOT)
            inputs = {} if inputs is None else inputs

            used_inputs = {}
            for category, variant in inputs.items():
                leaf = self.get_leaf(variant)
                if category not in self.nodes:
                    raise AscotMeltdownError(
                        "Output uses undefined input category: " + category
                        )
                if leaf not in self.nodes[category]:
                    raise AscotMeltdownError(
                        "Output uses input which is not present in the file."
                        )
                used_inputs[category] = leaf

            leaf = Leaf.create_leaf(
                date=date, note=note, variant=name.split("_")[0],
                inputs=used_inputs
                )
            if not isinstance(leaf, OutputVariant):
                raise AscotMeltdownError("Output variant expected: " + name)
            leaf._name = name
            self.enter_leaf(leaf, save=False, new_entry=False)
            leaf._load(self.access_data(name, ROOT))
        if active:
            self.get_leaf(active).activate()

    def save_leaf(self, leaf: Leaf) -> DataAccess:
        """Save the leaf to disk along with its metadata.

        This method just creates the group in file and stores attributes. To
        store actual data, a subsequent call to Leaf._save_data is required.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            Leaf whose metadata will be stored to disk.

        Raises
        ------
        :class:`.AscotDataException`
            If no file has been specified or if trying to save output before
            all inputs are saved.
        """
        if not self.filename:
            raise AscotDataException(
                "The data was not saved to disk because no file has been "
                "specified."
                )

        def get_inputs(output_leaf: OutputVariant) -> dict[str, InputVariant]:
            inputs = {}
            for input_leaf in self.uses[output_leaf]:
                for category, node in self.nodes.items():
                    if input_leaf in node:
                        inputs[category] = input_leaf
                        break
                if not input_leaf.status & Status.SAVED:
                    raise AscotDataException(
                        f"Cannot save output '{output_leaf.name}' before all "
                        f"inputs have been saved: {input_leaf.name}"
                        )

            return inputs

        for category, node in self.nodes.items():
            if leaf in node:
                break
        category = self.get_parentnode(leaf)
        if isinstance(leaf, OutputVariant) and leaf in self.nodes[ROOT]:
            inputs = get_inputs(leaf)
            name_by_category = {
                category: leaf.name for category, leaf in inputs.items()
                }
            self.write_variant(
                leaf.name, (leaf.date, leaf.note, name_by_category), ROOT,
                )
        else:
            self.write_variant(leaf.name, (leaf.date, leaf.note, None), category)
        if leaf is self.nodes[category].active:
            self.set_active(category, leaf.name)
        return self.access_data(leaf.name, category)

    #pylint: disable=too-many-arguments
    def enter_leaf(
            self, leaf: Leaf,
            category: str=ROOT,
            activate: bool=False,
            save: Optional[bool]=None,
            new_entry: bool=True,
            ) -> None:
        """Make the given leaf part of this tree.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf to add to this tree.
        category : str, *optional*
            Category to which the leaf belongs.

            Default value is the root node which applies for outputs, so this
            parameter is only required for inputs.
        activate : bool, *optional*
            Set the created leaf as active.
        save : bool, *optional*
            Store the new input variant in the HDF5 file (requires that the the
            tree has access to a file).

            If not explicitly set, the default behavior is to:

            - Store the input to a file if filename was specified when this tree
              was created.
            - Not to store the input if filename was not specified.
        new_entry : bool, *optional*
            This is a new leaf and not initialized from file.

        Raises
        ------
        :class:`.AscotMeltdownError`
            If the leaf already belongs to this tree or one or more of the
            used inputs (in case of an output variant) are not part of this
            tree.
        """
        if save is None:
            save = self.filename != ""

        if leaf._treemanager is not None:
            raise AscotMeltdownError(
                f"'{leaf.name}' already belongs to a tree."
                )

        def rename(leaf: Leaf) -> None:
            if not new_entry:
                return
            counter = 1
            for existing_leaf in (self.inputs + self.outputs):
                if existing_leaf.variant == leaf.variant:
                    counter = max(
                        counter, int(existing_leaf.name.split("_")[-1]) + 1
                        )
            leaf._name = f"{leaf.variant}_{counter}"

        if isinstance(leaf, OutputVariant):

            inputs, categories = [], []
            for cat in leaf._inputs:
                input_leaf = leaf[cat]
                if input_leaf not in self.inputs:
                    raise AscotMeltdownError(
                        f"Output {leaf.name} uses input {input_leaf.name} not "
                        f"present in the tree."
                        )
                inputs.append(input_leaf)
                categories.append(cat)

            rename(leaf)
            self.nodes[category]._add_leaf(leaf)
            self.uses[leaf] = []
            for input_leaf in inputs:
                self.usedin[input_leaf].append(leaf)
                self.uses[leaf].append(input_leaf)
            self.outputs.append(leaf)
        if isinstance(leaf, InputVariant):
            rename(leaf)
            self.nodes[category]._add_leaf(leaf)
            self.inputs.append(leaf)
            self.usedin[leaf] = []
        leaf._treemanager = self
        if activate:
            leaf.activate()
        if save:
            leaf._file = self.save_leaf(leaf)
            leaf._save_and_free_data()

    def destroy_leaf(self, leaf: Leaf, repack: Optional[bool]=None) -> None:
        """Destroy a single leaf.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf to destroy
        repack : bool, *optional*
            Repack the file afterwards.

        Raises
        ------
        :class:`.AscotDataException`
            If trying to remove input which is currently used.
        :class:`.AscotMeltdownError`
            If the leaf does not belong to this tree.
        """
        if isinstance(leaf, InputVariant) and leaf in self.inputs:
            if self.usedin[leaf]:
                names = [sim.name for sim in self.usedin[leaf]]
                raise AscotDataException(
                    f"Input '{leaf.name}' cannot be removed before the "
                    f"simulation(s) '{names}' has been removed since "
                    f"this input is being used there.")
            self.inputs.remove(leaf)
            node = self.get_parentnode(leaf)
            self.nodes[node]._remove_leaf(leaf)
        elif isinstance(leaf, OutputVariant) and leaf in self.outputs:
            self.outputs.remove(leaf)
            for input_leaf in self.uses[leaf]:
                self.usedin[input_leaf].remove(leaf)
            del self.uses[leaf]
            self.nodes[ROOT]._remove_leaf(leaf)
            node = ROOT
        else:
            raise AscotMeltdownError(f"Leaf '{leaf.name}' not found.")

        if self.filename and (leaf.status & Status.SAVED):
            self.remove_variant(leaf.name, node)
            try:
                active_leaf = self.nodes[node].active
                name = active_leaf.name
            except AscotDataException:
                name = ""
            self.set_active(node, name)
            if repack:
                self.repack()

    def activate_leaf(self, leaf: Leaf) -> None:
        """Set given leaf as active within the node it belongs to.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf to mark as active.
        """
        for category, node in self.nodes.items():
            if leaf in node:
                break
        else:
            category = ROOT
        self.nodes[category]._activate_leaf(leaf)
        if self.filename and leaf.status & Status.SAVED:
            self.set_active(category, leaf.name)

    def get_leaf(self, name: str) -> Leaf:
        """Retrieve leaf with the given name.

        Parameters
        ----------
        name : str
            Name of the desired leaf.

        Returns
        -------
        leaf : :class:`.Leaf`
            The leaf with the matching name.

        Raises
        ------
        :class:`.AscotDataException`
            If the leaf was not found.
        """
        for leaf in (self.inputs + self.outputs):
            if leaf.name == name:
                return leaf
        raise AscotDataException(f"Data '{name}' not found.")

    def get_parentnode(self, leaf: Leaf | str) -> str:
        """Retrieve the node given leaf belongs to.

        Parameters
        ----------
        leaf : `Leaf` or str
            Leaf or it's name.

        Returns
        -------
        name : str
            Name of the node.

        Raises
        ------
        ValueError
            If the leaf does not belong to any node.
        """
        for name, node in self.nodes.items():
            if leaf in node:
                return name
        raise ValueError("Leaf does not belong to any node.")

    def note_changed(self, leaf: Leaf) -> None:
        """Respond to notification that given leaf's note has changed.

        The note is used to create a tag which enables creating a custom
        reference to the leaf. Therefore, updating the note requires updating
        the references in the node it belongs to.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf whose note has changed.
        """
        for category, node in self.nodes.items():
            if leaf in node:
                break
        else:
            category = ROOT
        self.nodes[category]._organize()
        if self.filename and leaf.status & Status.SAVED:
            self.set_note(leaf.name, category, leaf.note)
