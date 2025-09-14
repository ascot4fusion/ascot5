"""Module for building a tree structure for accessing data.

The tree-like structure looks like this:

    Tree
    ├── bfield
    │   ├── B_TC_0123456789
    |   |   ├── bxyz
    │   │   └── Other *datasets*
    │   └── Other input *variants* in <variant>_<qid> format
    ├── Other input *categories*
    ├── run_2345678901
    │   ├── inistate
    │   │   ├── weight
    │   │   └── Other *datasets*
    │   └── Other *diagnostics*
    └── Other run *variants* in <variant>_<qid> format

The tree maintains the tree structure and manages the metadata within the Tree
object itself, input categories, and input and output variants. Only the latter
two access the actual data directly, which is why they are also referred as
"data".

The main purpose of the tree is to help the user to navigate and manipulate the
data using the public attributes of each component. Here's what references
different components have:

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

import warnings
from typing import Tuple, List, Dict, Optional

from a5py import utils
from a5py.exceptions import AscotIOException
from .hdf5 import HDF5Manager
from .leaf import MetaData, Leaf, Status, get_qid
from .nodes import ImmutableNode, InputCategory, OutputLeaf
from .dataholder import Format


ROOT = "root"
"""Name of the root/output node."""


class TreeManager(HDF5Manager):
    """Manages contents of the tree.

    This manager operates in the background and is responsible for keeping
    the contents of the tree consistent. The manager has references to all
    nodes and their leafs. To avoid cyclic calls, this manager is(only) called
    when:

    - The tree builds its nodes during initialization (call __init__).
    - New leaf is added to the tree (call enter_input or enter_run).
    - Leaf destroys itself (call remove_leaf).
      - When node destroys a leaf, it should call the leaf's destroy method and
        not do else.
    - Leaf activates itself (call activate_leaf).
    - The note on a leaf is updated (call note_changed) since this may update
      the tag that is used to reference the leaf.

    Attributes
    ----------
    inputs : List[Leaf]
        Leafs that belong to the input nodes.
    outputs : List[Leaf]
        Leafs that belong to the root node.
    usedin : Dict[MetaDataHolder, List[MetaDataHolder]]
        Maps input leafs to the outputs leafs that have used them.
    uses : Dict[MetaDataHolder, List[MetaDataHolder]]
        Maps output leaf to the inputs it has used.
    nodes : Dict[str, ImmutableNode]
        All nodes and their names that make up the tree.
    """

    def __init__(
            self, hdf5file: Optional[Tuple[str, bool]]=None,
            **nodes: ImmutableNode,
            ) -> None:
        """Initialize the tree.

        Parameters
        ----------
        hdf5file : Tuple[str, bool], optional
            Name of the HDF5 file and if the file exists.

            If this value is given, the data will be stored and read from the
            file (not memory).
        **nodes : ImmutableNode
            Nodes belonging to this tree (including the root node).
        """
        self.inputs: List[Leaf] = []
        self.outputs: List[Leaf] = []
        self.usedin: Dict[Leaf, List[Leaf]] = {}
        self.uses: Dict[Leaf, List[Leaf]] = {}
        self.nodes: Dict[str, ImmutableNode] = nodes
        for node in nodes.values():
            with node._modify_attributes():
                node._treemanager = self

        self.filename = None
        if hdf5file:
            filename, file_exists = hdf5file
            super().__init__(
                root=ROOT, filename=filename, file_exists=file_exists,
                input_categories=[n for n in nodes.keys() if n != ROOT]
                )
            self._init_from_hdf5([n for n in nodes.keys() if n != ROOT])

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(inputs={self.inputs}, "
            f"outputs={self.outputs})>"
        )

    def save(self, leaf):
        """Save the leaf to disk along with its metadata.

        This method just creates the group in file and stores attributes. To
        store actual data, a subsequent call to Leaf.save is required.

        Parameters
        ----------
        leaf : Leaf
            Leaf whose metadata will be stored to disk.
        """
        if self.filename is None:
            raise AscotIOException(
                "The data was not saved to disk because no file has been "
                "specified."
                )

        def get_inputs():
            inputs = {}
            for input_leaf in self.uses[leaf]:
                for category, node in self.nodes.items():
                    if input_leaf in node:
                        inputs[category] = input_leaf
                        break
                if not(input_leaf.status & Status.SAVED):
                    raise AscotIOException(
                        f"The data was created but not stored in the disk "
                        f"since it uses input with QID '{input_leaf.qid}' "
                        f"({category}) which is not stored in the file. Store "
                        f"the input first and only then this simulation data "
                        f"can be stored in the file.",
                        )

            return inputs

        for category, node in self.nodes.items():
            if leaf in node:
                break
        if leaf in self.nodes[ROOT]:
            inputs = get_inputs()
            qid_by_category = {
                category: leaf.qid for category, leaf in inputs.items()
                }
            self.write_output(leaf._metadata, qid_by_category)
            leaf._file = self.get_minimanager(
            ROOT, leaf.variant, leaf.qid
            )
        else:
            self.write_input(leaf._metadata, category)
        if leaf is self.nodes[category].active:
            self.set_active(category, leaf.qid)

    def _init_from_hdf5(self, input_categories):
        """Initializes the tree structure from the HDF5 file."""
        for category in input_categories:
            active, qids, variants = self.read_node(category)
            for qid, variant in zip(qids, variants):
                date, note = self.read_input(qid, variant, category)
                leaf = Leaf.create_leaf(MetaData(qid=qid, date=date, note=note, variant=variant))
                self.enter_input(leaf, category, save=False)
                leaf._file = self.get_minimanager(
                    category, variant, qid,
                )

            if active:
                self.get_leaf(active).activate()

        active, qids, variants = self.read_node(ROOT)
        for qid, variant in zip(qids, variants):
            date, note, inputqids = self.read_output(qid, variant)
            inputs= {}
            for leaf in [self.get_leaf(qid) for qid in inputqids]:
                for category, node in self.nodes.items():
                    if leaf in node:
                        inputs[category] = leaf
            leaf = Leaf.create_leaf(
                MetaData(qid=qid, date=date, note=note, variant=variant),
                inputs=inputs
                )
            self.enter_output(leaf, save=False)
            leaf._load(self.get_minimanager(
                ROOT, variant, qid,
            ))
        if active:
            self.get_leaf(active).activate()

    def enter_input(
            self,
            leaf: Leaf,
            category: str,
            activate: bool=False,
            dryrun: bool=False,
            save: Optional[bool]=None,
            ):
        """Add Leaf to a InputCategory within this tree.

        Parameters
        ----------
        leaf : `Leaf`
            Metadata of the input variant.
        activate : bool, optional
            Set the created leaf as active.
        dryrun : bool, optional
            Creates a leaf and checks if one with same QID exists on the tree,
            but doesn't include it into the tree or store it anywhere.
        save : bool, optional
            Store the new input variant in the HDF5 file (requires that the the
            tree has access to a file).

            If not explicitly set, the default behavior is to:
            - store the input to a file if filename was specified when this tree
              was created.
            - not to store the input if filename was not specified.

        Returns
        -------
        leaf : `Leaf`
            The input variant which was created.

        Raises
        ------
        AscotIOException
            If unable to add or store the leaf.
        """
        for existing_leaf in (self.inputs + self.outputs):
            if existing_leaf.qid == leaf.qid:
                raise AscotIOException(
                    f"Input variant with QID = {leaf.qid} already belongs to "
                    "this tree."
                )
        if dryrun:
            return leaf

        self.inputs.append(leaf)
        self.usedin[leaf] = []
        self.nodes[category]._add_leaf(leaf)
        leaf._treemanager = self
        if activate:
            leaf.activate()
        if save is None:
            save = self.filename is not None
        if save:
            if not self.filename:
                raise AscotIOException("No HDF5 file was provided.")
            leaf.save()

    def enter_output(
            self,
            leaf: OutputLeaf,
            activate: bool=False,
            dryrun: bool=False,
            save: Optional[bool]=None,
            ):
        """Add OutputLeaf to this tree.

        Parameters
        ----------
        meta : `MetaData`
            Metadata of the run variant.
        inputqids : List[str]
            QIDs of the input variants used in the new run.
        activate : bool, optional
            Set the created leaf as active.
        dryrun : bool, optional
            Creates a leaf and checks if one with same QID exists on the tree,
            but doesn't include it into the tree or store it anywhere.
        save : bool, optional
            Store the new run variant in the HDF5 file (requires that the the
            tree has access to a file).

            If not explicitly set, the default behavior is to:
            - store the input to a file if filename was specified when this tree
              was created.
            - not to store the input if filename was not specified.

        Returns
        -------
        leaf : `Leaf`
            The output variant which was created.

        Raises
        ------
        AscotIOException
            If unable to add or store the leaf.
        """
        inputqids = [leaf[cat].qid for cat in leaf._inputs]
        if save is None:
            save = self.filename is not None

        for existing_leaf in (self.inputs + self.outputs):
            if existing_leaf.qid == leaf.qid:
                raise AscotIOException(
                    f"Data with QID '{leaf.qid}' already belongs to this "
                    f"tree."
                )
        not_in_tree = [qid for qid in inputqids if qid in self.inputs]
        if not_in_tree:
            raise AscotIOException(
                f"Simulation with QID '{leaf.qid}' uses input(s) with QID "
                f"'{not_in_tree}' which is not present in the tree."
            )
        inputqid_is_output = [qid for qid in inputqids if qid in self.outputs]
        if inputqid_is_output:
            raise AscotIOException(
                f"Input with QID '{leaf.qid}' was used, but it corresponds to "
                f"a simulation - not input. Something's gone horribly wrong!"
            )
        inputs: List[Leaf] = []
        categories: List[str] = []
        for qid in inputqids:
            input_leaf = self.get_leaf(qid)
            for category, node in self.nodes.items():
                if input_leaf in node:
                    categories.append(category)
                    break
            inputs.append(input_leaf)
        if dryrun:
            return

        self.uses[leaf] = []
        for input_leaf in inputs:
            self.usedin[input_leaf].append(leaf)
            self.uses[leaf].append(input_leaf)
        self.outputs.append(leaf)
        self.nodes[ROOT]._add_leaf(leaf)
        leaf._treemanager = self
        if activate:
            leaf.activate()
        if save:
            self.save(leaf)

    def destroy_leaf(self, leaf: Leaf, repack: Optional[bool]=None) -> None:
        """Destroy a single leaf.

        Parameters
        ----------
        leaf : Leaf
            The leaf to destroy
        repack : bool, optional
            Repack the file afterwards.

        Raises
        ------
        AscotIOException
            If the leaf does not belong to this tree.
        """
        node = None
        if leaf in self.inputs:
            if self.usedin[leaf]:
                qids = [sim.qid for sim in self.usedin[leaf]]
                raise AscotIOException(
                    f"Input with QID '{leaf.qid}' cannot be removed before the "
                    f"simulation(s) with QID '{qids}' has been removed since "
                    f"this input is being used there.")
            self.inputs.remove(leaf)
            for category, node in self.nodes.items():
                if leaf in node:
                    self.nodes[category]._remove_leaf(leaf)
                    node = category
                    break
        elif leaf in self.outputs:
            self.outputs.remove(leaf)
            for input_leaf in self.uses[leaf]:
                self.usedin[input_leaf].remove(leaf)
            del self.uses[leaf]
            self.nodes[ROOT]._remove_leaf(leaf)
            node = ROOT
        else:
            raise AscotIOException(f"Leaf with QID '{leaf.qid}' not found.")

        if not self.filename is None and (leaf.status & Status.SAVED):
            self.remove_variant(leaf.qid, leaf.variant, node)
            qid = self.nodes[node]._active
            qid=""
            self.set_active(node, qid)
            if repack:
                self.repack()

    def activate_leaf(self, leaf: Leaf) -> None:
        """Set given leaf as active within the node it belongs to.

        Parameters
        ----------
        leaf : Leaf
            The leaf to mark as active.
        """
        for category, node in self.nodes.items():
            if leaf in node:
                break
        else:
            category = ROOT
        self.nodes[category]._activate_leaf(leaf)
        if not self.filename is None and leaf.status & Status.SAVED:
            self.set_active(category, leaf.qid)

    def get_leaf(self, qid: str) -> Leaf:
        """Retrieve leaf with the given QID.

        Parameters
        ----------
        qid : str
            QID of the desired leaf.

        Returns
        -------
        leaf : Leaf
            The leaf with the matching QID.
        """
        for leaf in (self.inputs + self.outputs):
            if leaf.qid == qid:
                return leaf
        raise AscotIOException(f"Data with qid = {qid} not found.")

    def get_parentnode(self, leaf: Leaf | str) -> str:
        """Retrieve the node given leaf belongs to.

        Parameters
        ----------
        leaf : Leaf or str
            Leaf or it's QID.
        """
        for name, node in self.nodes.items():
            if leaf in node:
                return name
        raise ValueError("Leaf does not belong to any node.")

    def note_changed(self, leaf: Leaf) -> None:
        """Notify that given leaf's note has changed.

        The note is used to create a tag which enables creating a custom
        reference to the leaf. Therefore, updating the note requires updating
        the references in the node it belongs to.

        Parameters
        ----------
        leaf : Leaf
            The leaf whose note has changed.
        """
        for category, node in self.nodes.items():
            if leaf in node:
                break
        else:
            category = ROOT
        self.nodes[category]._organize()
        if not self.filename is None and leaf.status & Status.SAVED:
            self.set_note(leaf.qid, leaf.variant, category, leaf.note)


class Tree(ImmutableNode):
    """The entry node for accessing simulation data."""

    def __init__(
            self, input_categories,
            hdf5file: Optional[Tuple[str, bool]]=None
            ) -> None:
        """Initialize an empty tree structure.

        Parameters
        ----------
        hdf5file : Tuple[str, bool], optional
            Filename of the HDF5 file, if used, and flag indicating if the file
            exists.
        """
        super().__init__()
        nodes = {ROOT:self}
        self._input_categories = input_categories
        for category in input_categories:
            self[category] = InputCategory()
            nodes[category] = self[category]
            self[category]._freeze()

        TreeManager(hdf5file=hdf5file, **nodes)

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(qids={self._qids}, tags={self._tags}, "
            f"active={self._active!r})>"
            )

    def _get_decorated_contents(self) -> str:
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        def print_category(category):
            return utils.decorate(f"{category.ljust(10)}", color="green")

        def print_howmanyinputs(number_of_inputs):
            if number_of_inputs > 1:
                return f" + {number_of_inputs-1} other(s)"
            if number_of_inputs == 1:
                return " (no other inputs)"
            return "*no inputs*\n\n"

        def print_name(leaf):
            if not leaf:
                return ""
            return (
                utils.decorate(
                    f"{leaf.variant.ljust(10)}{leaf.qid}", bold=True,
                )
                + f" {leaf.date}"
            )

        def print_note(leaf):
            if not leaf:
                return ""
            return f"\n{''.ljust(10)}\"{leaf.note}\"\n"

        def print_title(title):
            return utils.decorate(
                title, color="purple", underline=True, bold=True,
                )

        contents = ""
        contents += print_title("Inputs:")
        contents += utils.decorate(" [only active shown]\n", color="green")
        for category in self._input_categories:
            contents += print_category(category)
            leaf = None
            try:
                leaf = self[category].active
                contents += print_name(leaf)
            except AscotIOException:
                pass

            n_inp = len(self[category]._qids) # pylint: disable=protected-access
            contents += print_howmanyinputs(n_inp)
            contents += print_note(leaf)

        contents += print_title("\nSimulations:\n")
        for leaf in self:
            contents += print_name(leaf)
            if leaf == self.active:
                contents += utils.decorate(" [active]", color="green")
            contents += print_note(leaf)

        if not self._qids:
            contents += "No simulation results.\n"
        return contents

    @property
    def contents(self) -> str:
        """A string representation of the contents."""
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self) -> None:
        """Show on screen the metadata of the currently active inputs and list
        of all simulation runs.
        """
        print(self._get_decorated_contents())

    def ls(self, show: bool=False) -> str:
        """Get a string representation of the contents.

        Deprecated. Use `show_contents()` or `contents` instead.

        Parameters
        ----------
        show : str, optional
            If True, the contents are also printed on screen.

        Returns
        -------
        contents : str
            Multiline string decorated with ANSI escape sequences that list
            this node's meta data, all output data within this node, and inputs
            that were used.
        """
        warnings.warn(
            "'ls' will be removed in a future release. Use 'show_contents' "
            "instead.",
            DeprecationWarning, stacklevel=2)
        contents = self._get_decorated_contents()
        if show:
            print(contents)
        return contents

    def activate(self, data: Leaf | str) -> None:
        """Mark data as active.

        Parameters
        ----------
        data : `Leaf` or str
            Data to be activated (or it's QID or name).
        """
        if self._treemanager:
            if isinstance(data, str):
                data = self._treemanager.get_leaf(data)
            self._treemanager.activate_leaf(data)

    def destroy(
            self, *, repack: bool=True, data: Leaf | str=ROOT, **kwargs,
            ) -> None:
        """Remove data permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Repacking has some overhead but without it only the references to
            the data, not the data itself, are removed in the file.
        data : `Leaf` or str or None
            Data to be removed (or it's QID or name or entire category).

            The default value removes all results.
        """
        if data == ROOT:
            super().destroy(repack=repack, **kwargs)
        elif isinstance(data, str) and data in self:
            self[data].destroy(repack=repack)
        else:
            qid = get_qid(data)
            if self._treemanager:
                leaf = self._treemanager.get_leaf(qid)
                self._treemanager.destroy_leaf(leaf, repack=repack)
