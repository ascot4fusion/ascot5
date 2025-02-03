"""Contains classes that make up the tree data structure storing ASCOT5 data.

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

    Node (root)
    ├── Node (input category)
    │   └── Leaf (input variant)
    |       └── Public attributes and methods to access the actual data
    └ Leaf (run variant)
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

from typing import Any, Tuple, List, Dict, Union, Generator, Optional, Protocol

from ...exceptions import AscotIOException
from . import metadata
from .hdf5 import HDF5Manager
from .metadata import MetaData, MetaDataHolder
from .dataholder import Format
from .immutablestorage import ImmutableStorage


#pylint: disable=too-few-public-methods
class LeafFactory(Protocol):
    """Signature for leaf factory methods (used in type hints)."""

    def __call__(
            self,
            meta : MetaData,
            inputs : Optional[Dict[str, Leaf]] = None,
            diagnostics : Optional[List] = None,
            **kwargs: Any,
            ) -> Leaf:
        pass

class TreeManager():
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
    - The note on a leaf is updated (call note_changed).

    Attributes
    ----------
    leaf_factory : callable
        Function that uses the metadata to create a desired `Leaf` instance.
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
    storedin : Dict[MetaDataHolder, DataManager.Format]
        Where the data is being stored.
    hdf5manager : HDF5Manager
        Manages the contents of the HDF5 file if applicable.
    """

    def __init__(
            self,
            leaf_factory: LeafFactory,
            hdf5file: Optional[Tuple[str, bool]] = None,
            **nodes: ImmutableNode,
            ) -> None:
        """Initialize the tree.

        Parameters
        ----------
        leaf_factory : callable
            The factory method which creates new leafs.
        hdf5file : Tuple[str, bool], optional
            Name of the HDF5 file and if the file exists.

            If this value is given, the data will be stored and read from the
            file (not memory).
        **nodes : ImmutableNode
            Nodes belonging to this tree (including the root node).
        """
        self.leaf_factory = leaf_factory
        self.inputs: List[Leaf] = []
        self.outputs: List[Leaf] = []
        self.usedin: Dict[Leaf, List[Leaf]] = {}
        self.uses: Dict[Leaf, List[Leaf]] = {}
        self.nodes: Dict[str, ImmutableNode] = nodes
        for node in nodes.values():
            with node._modify_attributes():
                node._treemanager = self

        self.storedin: Dict[Leaf, Format] = {}
        self.hdf5manager: Optional[HDF5Manager] = None
        if hdf5file:
            filename, file_exists = hdf5file
            self.hdf5manager = HDF5Manager(filename, file_exists)

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(inputs={self.inputs}, "
            f"outputs={self.outputs})>"
        )

    def init_from_hdf5(self):
        """Initializes the tree structure from the HDF5 file."""
        for category in metadata.input_categories:
            active, qids, variants = self.hdf5manager.read_node(category)
            for qid, variant in zip(qids, variants):
                date, note = self.hdf5manager.read_input(qid, variant)
                leaf = self.enter_input(
                    MetaData(qid=qid, date=date, note=note, variant=variant),
                    store_hdf5=False,
                    )
                self.storedin[leaf] = Format.HDF5

            if active:
                self.get_leaf(active).activate()

        active, qids, variants = self.hdf5manager.read_node("root")
        for qid, variant in zip(qids, variants):
            date, note, inputqids = self.hdf5manager.read_run(qid, variant)
            leaf = self.enter_run(
                MetaData(qid=qid, date=date, note=note, variant=variant),
                [],
                inputqids,
                store_hdf5=False,
                )
            self.storedin[leaf] = Format.HDF5
        if active:
            self.get_leaf(active).activate()

    def enter_input(
            self,
            meta: MetaData,
            activate: bool=False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> Leaf:
        """Create and add a new input variant to this tree.

        This method uses the factory method to create the `Leaf` instance that
        is stored.

        Parameters
        ----------
        meta : `MetaData`
            Metadata of the input variant.
        activate : bool, optional
            Set the created leaf as active.
        dryrun : bool, optional
            Creates a leaf and checks if one with same QID exists on the tree,
            but doesn't include it into the tree or store it anywhere.
        store_hdf5 : bool, optional
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
        leaf = self.leaf_factory(meta)
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
        category = metadata.get_input_category(leaf.variant)
        self.nodes[category]._add_leaf(leaf)
        leaf._treemanager = self
        self.storedin[leaf] = Format.CSTRUCT
        if activate:
            leaf.activate()
        if store_hdf5 is None:
            store_hdf5 = self.hdf5manager is not None
        if store_hdf5:
            if not self.hdf5manager:
                raise AscotIOException("No HDF5 file was provided.")
            self.hdf5manager.write_input(leaf._metadata)
            self.storedin[leaf] = Format.HDF5
        return leaf

    def enter_run(
            self,
            meta: MetaData,
            diagnostics: List [str],
            inputqids: List [str],
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> Leaf:
        """Create and add a new run variant to this tree.

        Parameters
        ----------
        meta : `MetaData`
            Metadata of the run variant.
        diagnostics : List[str]
        inputqids : List[str]
            QIDs of the input variants used in the new run.
        activate : bool, optional
            Set the created leaf as active.
        dryrun : bool, optional
            Creates a leaf and checks if one with same QID exists on the tree,
            but doesn't include it into the tree or store it anywhere.
        store_hdf5 : bool, optional
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
        if store_hdf5 is None:
            store_hdf5 = self.hdf5manager is not None
        for existing_leaf in (self.inputs + self.outputs):
            if existing_leaf.qid == meta.qid:
                raise AscotIOException(
                    f"Data with QID '{meta.qid}' already belongs to this "
                    f"tree."
                )
        not_in_tree = [qid for qid in inputqids if qid in self.inputs]
        if not_in_tree:
            raise AscotIOException(
                f"Simulation with QID '{meta.qid}' uses input(s) with QID "
                f"'{not_in_tree}' which is not present in the tree."
            )
        inputqid_is_output = [qid for qid in inputqids if qid in self.outputs]
        if inputqid_is_output:
            raise AscotIOException(
                f"Input with QID '{meta.qid}' was used, but it corresponds to "
                f"a simulation - not input. Something's gone horribly wrong!"
            )
        inputs: List[Leaf] = []
        categories: List[str] = []
        for qid in inputqids:
            input_leaf = self.get_leaf(qid)
            categories.append(metadata.get_input_category(input_leaf.variant))
            inputs.append(input_leaf)
        input_objs = dict(zip(categories, inputs))
        leaf = self.leaf_factory(
            meta, inputs=input_objs, diagnostics=diagnostics
            )
        if dryrun:
            return leaf

        self.uses[leaf] = []
        for input_leaf in inputs:
            self.usedin[input_leaf].append(leaf)
            self.uses[leaf].append(input_leaf)
        self.outputs.append(leaf)
        self.nodes["root"]._add_leaf(leaf)
        leaf._treemanager = self
        self.storedin[leaf] = Format.CSTRUCT
        if activate:
            leaf.activate()
        if store_hdf5:
            if not self.hdf5manager:
                raise AscotIOException(
                    "The data was created but not stored in the disk since "
                    "no file has been specified."
                    )
            for category, input_leaf in zip(categories, inputs):
                if self.storedin[input_leaf] != Format.HDF5:
                    raise AscotIOException(
                        f"The data was created but not stored in the disk "
                        f"since it uses input with QID '{input_leaf.qid}' "
                        f"({category}) which is not stored in the file. Store "
                        f"the input first and only then this simulation data "
                        f"can be stored in the file.",
                        )
            qid_by_category = dict(zip(categories, inputqids))
            self.hdf5manager.write_run(leaf._metadata, qid_by_category)
            self.storedin[leaf] = Format.HDF5
        return leaf

    def destroy_leaf(self, leaf: Leaf, repack: Optional[bool] = None) -> None:
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
        if leaf in self.inputs:
            if self.usedin[leaf]:
                qids = [sim.qid for sim in self.usedin[leaf]]
                raise AscotIOException(
                    f"Input with QID '{leaf.qid}' cannot be removed before the "
                    f"simulation(s) with QID '{qids}' has been removed since "
                    f"this input is being used there.")
            self.inputs.remove(leaf)
            category = metadata.get_input_category(leaf.variant)
            self.nodes[category]._remove_leaf(leaf)
        elif leaf in self.outputs:
            self.outputs.remove(leaf)
            for input_leaf in self.uses[leaf]:
                self.usedin[input_leaf].remove(leaf)
            del self.uses[leaf]
            self.nodes["root"]._remove_leaf(leaf)
        else:
            raise AscotIOException(f"Leaf with QID '{leaf.qid}' not found.")

        if self.hdf5manager and self.storedin[leaf] == Format.HDF5:
            self.hdf5manager.remove_variant(leaf.qid, leaf.variant)
            if repack:
                self.hdf5manager.repack()
        del self.storedin[leaf]

    def activate_leaf(self, leaf: Leaf) -> None:
        """Set given leaf as active within the node it belongs to.

        Parameters
        ----------
        leaf : Leaf
            The leaf to mark as active.
        """
        try:
            category = metadata.get_input_category(leaf.variant)
        except ValueError:
            category = "root"
        self.nodes[category]._activate_leaf(leaf)
        if self.hdf5manager and self.storedin[leaf] == Format.HDF5:
            self.hdf5manager.set_active(category, leaf.qid)

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
        try:
            category = metadata.get_input_category(leaf.variant)
        except ValueError:
            category = "root"
        self.nodes[category]._organize()
        if self.hdf5manager and self.storedin[leaf] == Format.HDF5:
            self.hdf5manager.set_note(leaf.qid, leaf.variant, leaf.note)


class Leaf(MetaDataHolder):
    """Leaf of a tree which has no further children managed by the tree.

    Attributes
    ----------
    _treemanager: TreeManager
        Manager of the tree this leaf belongs to.
    """

    def __init__(self, **kwargs) -> None:
        """Initialize leaf that does not yet belong to a tree."""
        super().__init__(**kwargs)
        self._treemanager: Optional[TreeManager] = None

    @property
    def note(self) -> str:
        """Short note for the user to document this data."""
        # Redefined here to please mypy
        return super().note

    @note.setter
    def note(self, note: str) -> None:
        """Set the note."""
        # Changing note changes tag as well.
        self._note = note
        if self._treemanager:
            self._treemanager.note_changed(self)

    def activate(self) -> None:
        """Set this data as active.

        Active inputs are used when the simulation is run. Active data variants
        are also used during post-processing by default unless otherwise
        specified.
        """
        if self._treemanager:
            self._treemanager.activate_leaf(self)

    def destroy(self, *, repack: bool = True) -> None:
        """Remove this data permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file to reduce it's disk size.

            Deleting data in an HDF5 file with h5py only removes its references,
            not the actual data. Repacking copies the data to a new file and
            replaces the original, freeing up space. May take a moment for large
            files.
        """
        if self._treemanager:
            self._treemanager.destroy_leaf(self, repack=repack)


class ImmutableNode(ImmutableStorage):
    """Tree node which can store other nodes or leaves, and whose attributes
    cannot be altered once frozen.

    This class provides the following main functionalities:

    1. Attributes can be accessed in a dictionary-like manner, e.g.,
       `node.child` and `node["child"]` are equivalent.
    2. Freezing this instance prevents setting and removing any attributes.
    3. Leaves can be accessed using their name, QID, or a tag constructed from
       the user-given note.
    4. One leaf is always set as 'active' and is accessed via the `active`
       attribute.
    5. Leaves can be iterated over and are organized by their date of creation.
    6. Whether a leaf belongs to this node can be checked with `leaf in node`,
       where `leaf` is a `Leaf`, `qid`, or `name`.

    The attributes can be modified using a context:

    .. code-block:: python

       with node._modify_attributes():
           node.new_attribute = "new_value"

    Attributes
    ----------
    _qids : list of str
        QIDs of this node's leaves sorted by date starting from newest.
    _tags : list of str
        List of all tags that can be used to access the leaves in the same order
        as QIDs.
    _active : `Leaf`
        The currently active leaf.
    _treemanager : `TreeManager`
        The manager of the tree this node belongs to.
    """

    def __init__(self, **kwargs) -> None:
        """Initialize an empty unfrozen node.

        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)
        self._treemanager: Optional[TreeManager] = None
        self._qids: List = []
        self._tags: List = []
        self._active: Optional[Leaf] = None

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(qids={self._qids}, tags={self._tags}, "
            f"active={self._active!r}, frozen={self._frozen})>"
            )

    def __contains__(self, key: Union[str, Leaf]) -> bool:
        """Check whether this node contains the requested leaf.

        Parameters
        ----------
        key : str or `Leaf`
            The requested leaf or it's QID.

        Returns
        -------
        contains : bool
            True if this node contains the leaf.
        """
        if isinstance(key, Leaf):
            key = key.qqid
        return hasattr(self, key) or hasattr(self, f"q{key}")

    def __iter__(self) -> Generator[Leaf, None, None]:
        """Iterate over this node's leafs."""
        for qid in self._qids:
            yield self[f"q{qid}"]

    def _activate_leaf(self, leaf) -> None:
        """Set given leaf as active.

        Parameters
        ----------
        leaf : `Leaf`
            The leaf to be set as active.

        Raises
        ------
        AscotIOException
            If the leaf does not belong to this node.
        """
        if leaf not in self:
            raise AscotIOException(
                f"Data with QID = {leaf.qid} does not belong to this node."
                )
        with self._modify_attributes():
            self._active = leaf

    def _add_leaf(self, leaf: Leaf) -> None:
        """Add leaf to this node.

        Parameters
        ----------
        leaf : `Leaf`
            The leaf to be added.

        Raises
        ------
        AscotIOException
            If the leaf already belongs to this node.
        """
        if leaf in self:
            raise AscotIOException(
                f"Data with QID = {leaf.qid} already belongs to this node."
                )
        self._qids.append(leaf.qid)
        with self._modify_attributes():
            reference_by_name, reference_by_qid = leaf.name, leaf.qqid
            self[reference_by_qid] = leaf
            self[reference_by_name] = leaf

        self._organize()

    def _remove_leaf(self, leaf: Leaf) -> None:
        """Remove a leaf from this node.

        Parameters
        ----------
        leaf : `Leaf`
            The leaf to be removed.
        """
        if leaf not in self:
            raise AscotIOException(
                f"Data with QID = {leaf.qid} does not belong to this node."
                )
        self._qids.remove(leaf.qid)
        with self._modify_attributes():
            if self.active == leaf:
                if self._qids:
                    self._active = self[f"q{self._qids[0]}"]
                else:
                    self._active = None
            reference_by_name, reference_by_qid = leaf.name, leaf.qqid
            delattr(self, reference_by_qid)
            delattr(self, reference_by_name)

        self._organize()

    def _organize(self) -> None:
        """Organize this node and its references when its leafs has been
        modified.

        - If this node has no leafs, do nothing except set active leaf to
          None.
        - If this node has just one leaf, set it as active.
        - Sort `self._qids` by date.
        - Update references by tag. If multiple leafs have the same tag, update
          their tags with running index, i.e. `new_tag = tag_<index>`, counting
          from zero for the leaf with the most recent date.
        """
        def activate_if_first_leaf():
            """Activate the leaf if it is the first one in this node."""
            try:
                self.active
            except AscotIOException:
                self._active = self[f"q{self._qids[0]}"]

        def sort_qids_by_date():
            """Sort the collected list of qids by date."""
            dates = [leaf.date for leaf in self]
            self._qids = [
                qid for _, qid in sorted(zip(dates, self._qids), reverse=True)
                ]

        def update_references_by_tag():
            """Add references by tag with unique tags and remove the old ones.
            """
            for tag_to_be_removed in self._tags:
                delattr(self, tag_to_be_removed)
            self._tags = []

            dates = [leaf.date for leaf in self]
            unsorted_qids = self._qids
            unsorted_tags = [leaf._extract_tag() for leaf in self]

            counts = {}
            dates.reverse()
            for tag, _, qid in sorted(zip(unsorted_tags, dates, unsorted_qids)):
                if tag in counts:
                    counts[tag] += 1
                    new_tag = f"{tag}_{counts[tag]}"
                else:
                    counts[tag] = 0
                    if unsorted_tags.count(tag) > 1:
                        new_tag = f"{tag}_{counts[tag]}"
                    else:
                        new_tag = tag

                self._tags.append(new_tag)
                self[new_tag] = self[f"q{qid}"]

        with self._modify_attributes():
            if len(self._qids) == 0:
                self._active = None
                for tag_to_be_removed in self._tags:
                    delattr(self, tag_to_be_removed)
                self._tags = []
                return

            activate_if_first_leaf()
            sort_qids_by_date()
            update_references_by_tag()

    @property
    def active(self) -> Leaf:
        """The active data.

        Raises
        ------
        AscotIOException
            If there's no active data.
        """
        if self._active is None:
            raise AscotIOException(
                "No active data. Perhaps this node is empty?"
                )
        return self._active

    def destroy(self, *, repack: bool = False, **kwargs) -> None:
        """Remove all data belonging to this node permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Repacking has some overhead but without it only the references to
            the data, not the data itself, are removed in the file.
        """
        _ = kwargs # kwargs are required only for subclasses to have kw args
        if self._treemanager:
            for qid in list(self._qids):
                leaf = self[f"q{qid}"]
                # Repacking takes time so repack only at the last item
                if len(self._qids) > 1:
                    self._treemanager.destroy_leaf(leaf, repack=False)
                else:
                    self._treemanager.destroy_leaf(leaf, repack=repack)
