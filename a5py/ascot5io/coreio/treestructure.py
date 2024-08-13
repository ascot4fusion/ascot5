"""Module for building a tree data structure representing ASCOT5 input and output.

This module provides a tree consisting of `ImmutableNode` and `Leaf` instances,
rooted at a `Root` instance. The tree can be traversed using attributes or
dictionary-like access (e.g., `root.child` or `root['child']`). Nodes are
immutable to prevent accidental modifications.

Components
----------
- `Root`: Manages the tree structure, holding input data nodes and simulation
  runs. Inherits from `ImmutableNode`.
- `ImmutableStorage`: Base class for all nodes, enforcing immutability.
- `ImmutableNode`: Node that stores `Leaf` instances.
- `InputCategory`: Inherits from `ImmutableNode`, representing different input
  categories (e.g., bfield, efield). Each `InputCategory` node can contain
  multiple `Leaf` instances.
- `Leaf`: Alias for the `MetaDataHolder` class, holding metadata such as
  identification, creation date, and variant. `Leaf` instances do not have
  children managed by `Root`.
- `SimulationOutput`: Inherits from `Leaf` and `ImmutableStorage`, representing
  the output of simulations. Contains references to input `Leaf` instances,
  which technically makes this structure a graph rather than a strict tree.
"""
from __future__ import annotations

import warnings
from contextlib import contextmanager
from typing import Tuple, List, Dict, Any, Generator, Optional

from ...exceptions import AscotIOException
from ... import utils
from . import metadata
from .metadata import MetaDataHolder, MetaData
from .datamanager import DataManager
from .hdf5interface import HDF5Manager


class Leaf(MetaDataHolder):
    """Leaf of a tree which has no further children managed by the tree.
    """

    def __init__(self, **kwargs):
        """Leaf which has no parent initially."""
        super().__init__(**kwargs)
        self._parent: ImmutableNode = None
        self._data_manager: DataManager = DataManager()
        self._hdf5manager: Optional[HDF5Manager] = None

    def _adopt(self, parent: ImmutableNode) -> None:
        """Adopt this leaf to a tree.

        Parameters
        ----------
        parent : `ImmutableNode`
            The node where this leaf belongs to.

        Raises
        ------
        AscotIOException
            If the leaf is already adopted.
        """
        if self._parent:
            raise AscotIOException("Leaf is already adopted.")
        self._parent = parent

    @MetaDataHolder.note.setter
    def note(self, note):
        # Changing note changes tag as well.
        self._note = note
        self._parent and self._parent._organize()
        if self._hdf5manager:
            self._hdf5manager.update_leaf_note(
                self._parent._name, self._qid, self._variant, note
                )

    def activate(self) -> None:
        """Set this dataset as active.

        Active inputs are used when the simulation is run. Active datasets are
        also used during post-processing by default unless otherwise specified.
        """
        self._parent and self._parent._activate_leaf(self)

    def destroy(self, repack: bool = True) -> None:
        """Remove this dataset.

        This also removes any data on the disk.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Without repacking only references to the data, not the actual data,
            are removed. Repacking has some overhead for large files and
            requires third-party `h5repack` tool.
        """
        self._parent and self._parent._remove_leaf(self)
        self._parent = None
        if self._hdf5manager:
            self._hdf5manager.remove_leaf(
                self._parent._name, self._qid, self._variant
                )
            repack and self._hdf5manager.repack()

    def export_hdf5(self):
        """Stores the data to the HDF5 file.
        """
        if self._hdf5manager:
            raise AscotIOException("This data is already stored in a file.")
        self._hdf5manager = self._parent._hdf5manager
        self._hdf5manager.write_leaf(
            MetaData(self._qid, self._date, self._variant, self._note)
        )
        filename, path_in_file = self._hdf5manager.locate_leaf(self._metadata)
        self._data_manager.migrate(
            hdf5_filename=filename, path_within_hdf5=path_in_file,
            )


class ImmutableStorage():
    """Object which supports dictionary-like assignment and which can be made
    immutable.

    Attributes
    ----------
    frozen : bool
        Indicates whether the node is frozen, preventing attribute modification.
    """

    def __init__(self, **kwargs: Any) -> None:
        """Initialize an empty node which is unfrozen.

        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        self._frozen: bool = False
        super().__init__(**kwargs)

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return f"<{self.__class__.__name__}(frozen={self._frozen})>"

    def __setitem__(self, key: str, value: Any) -> None:
        """Add a new attribute this node in dictionary style.

        Parameters
        ----------
        key : str
            Name of the attribute.
        value
            Value of the attribute.

        Raises
        ------
        AscotIOException
            Raised if the node is frozen.
        """
        if self._frozen:
            raise AscotIOException(
                "The attributes of this class are immutable."
                )
        setattr(self, key, value)

    def __setattr__(self, key: str, value: Any) -> None:
        """Add a new attribute this node.

        Parameters
        ----------
        key: str
            Name of the attribute.
        value:
            Value of the attribute.

        Raises
        ------
        AscotIOException
            Raised if the node is frozen.
        """
        if key != "_frozen" and self._frozen:
            raise AscotIOException(
                "The attributes of this class are immutable."
                )
        super().__setattr__(key, value)

    def __getitem__(self, key: str) -> Any:
        """Retrieve attribute in dictionary-like manner.

        Parameters
        ----------
        key : str
            Name of the attribute.

        Returns
        -------
        value
            Value of the attribute.
        """
        return getattr(self, key)

    def _freeze(self) -> None:
        """Make this node immutable.
        """
        self._frozen = True

    def _unfreeze(self) -> None:
        """Make this node mutable.
        """
        self._frozen = False

    @contextmanager
    def _modify_attributes(self) -> Generator[ImmutableStorage, None, None]:
        """Open a context where attributes can be modified.
        """
        self._unfreeze()
        try:
            yield self
        finally:
            self._freeze()


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
    _root : `Root`
        The root to which this node belongs.
    _name : str
        Name of this node (optional).
    """

    def __init__(self, root: Root, name: str = "", **kwargs) -> None:
        """Initialize an empty node which is unfrozen.

        Parameters
        ----------
        root : `Root`
            The root to which this node belongs.
        name : str, optional
            Name of this node.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)
        self._name: str = name
        self._root: Root = root
        self._qids: List = []
        self._tags: List = []
        self._active: Leaf = None

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(root={self._name}, "
            f"root={self._root!r}, qids={self._qids}, tags={self._tags}, "
            f"active={self._active!r}, frozen={self._frozen})>"
            )

    def __contains__(self, key: str) -> bool:
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
        try:
            key = key.qqid
        except AttributeError:
            pass
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
                f"Leaf with QID {leaf.qid} does not belong to this node."
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
                f"Leaf with QID {leaf.qqid} already belongs to this node."
                )

        self._qids.append(leaf.qid)

        with self._modify_attributes():
            reference_by_name, reference_by_qid = leaf.name, leaf.qqid
            self[reference_by_qid] = leaf
            self[reference_by_name] = leaf

        self._organize()
        leaf._adopt(self)

    def _remove_leaf(self, leaf: Leaf) -> None:
        """Remove a leaf from this node.

        Parameters
        ----------
        leaf : `Leaf`
            The leaf to be removed.
        """
        if leaf not in self:
            raise AscotIOException(
                f"Leaf with QID {leaf.qqid} does not belong to this node."
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
        """Organize this node and its references when its children has been
        modified.

        - If this node has no children, do nothing except set active group to
          None.
        - If this node has just one child, set it as active.
        - Sort `self._qids` by date.
        - Update references by tag. If multiple groups have the same tag, update
          their tags with running index, i.e. `new_tag = tag_<index>`, counting
          from zero for the group with the most recent date.
        """
        def activate_if_first_child():
            """Activate the group if it is the first one in this node."""
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

            activate_if_first_child()
            sort_qids_by_date()
            update_references_by_tag()

    @property
    def active(self) -> Leaf:
        """The active group.

        Raises
        ------
        AscotIOException
            If there's no active group.
        """
        if self._active is None:
            raise AscotIOException(
                "No active dataset. Perhaps there's no data in this node?"
                )
        return self._active

    def destroy(self) -> None:
        """Remove all datasets belonging to this node.
        """
        for dataset in self:
            dataset.destroy()


class ImmutableNodeHDF5(ImmutableNode):
    """Extends the `ImmutableNode` to support data storage in HDF5 format.

    Attributes
    ----------
    _hdf5manager : HDF5 manager
        Manages the contents of the HDF5 file.
    """

    def __init__(self, root: Root, **kwargs):
        """Initialize node with no manager."""
        super().__init__(root, **kwargs)
        self._hdf5manager = None

    def _activate_leaf(self, leaf: Leaf):
        """
        """
        super()._activate_leaf(leaf)
        if self._hdf5manager:
            self._hdf5manager.set_active(self._name, self._active)

    def _add_leaf(self, leaf: Leaf, store_hdf5: bool = False):
        """
        """
        super()._add_leaf(leaf)
        if store_hdf5:
            self._hdf5manager.write_leaf(self._name, leaf.name)

    def _remove_leaf(self, leaf : Leaf):
        """
        """
        super()._remove_leaf(leaf)
        if self._hdf5manager:
            self._hdf5manager.write_leaf(self._name, leaf.name)

    def _set_hdf5manager(self, manager: HDF5Manager):
        """Associate this node with a file.
        """
        self._hdf5manager = manager

    def destroy(self, repack: bool = True) -> None:
        """Remove all datasets and associated data permanently from within this
        node.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Without repacking only references to the data, not the actual data,
            are removed. Repacking has some overhead for large files and
            requires third-party `h5repack` tool.
        """
        super().destroy()
        if self._hdf5manager and repack:
            self._hdf5manager.repack()


class InputCategory(ImmutableNodeHDF5):
    """Node that contains all inputs of the same category."""

    def _get_decorated_contents(self):
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        contents = ""
        for index, leaf in enumerate(self):
            contents += utils.decorate(f"{leaf.variant.ljust(10)} {leaf.qid}",
                                       bold=True)
            contents += f" {leaf.date}"
            if self.active == leaf:
                contents += utils.decorate(" [active]", color="green")

            contents += f"\n{self._tags[index]}\n{leaf.note}\n\n"

        if not contents:
            contents = "No data in this category.\n"
        return contents

    @property
    def contents(self):
        """A string representation of the contents.
        """
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self):
        """Show on screen the metadata of all inputs within this category and
        which input is active."""
        print(self._get_decorated_contents())

    def ls(self, show=False):
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

    def destroy(self, repack: bool = True) -> None:
        """Remove all inputs and data within this category permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Without repacking only references to the data, not the actual data,
            are removed. Repacking has some overhead for large files and
            requires third-party `h5repack` tool.
        """
        # This function overriden only to update the docstring.
        super().destroy(repack=repack)


class SimulationOutput(ImmutableStorage, Leaf):
    """Node that contains results of a single simulation.

    This node contains all the metadata associated with the simulation.

    References to different diagnostics, that contain the actual data, are
    stored as '<underscore><name of the diagnostic>'.

    References to the input datasets used in the simulation are stored by
    category. For instance, `node.bfield` is the magnetic field dataset used in
    the simulation. It should point to the same object which can be found within
    the 'bfield' category in the tree.

    If trying to access an input or diagnostic not used in the simulation, an
    exception is raised.
    """

    def __init__(
            self,
            inputs: List[str, Leaf],
            diagnostics: Dict[str, Any],
            **kwargs: Any,
            ) -> None:
        """Initialize simulation output node with given inputs and diagnostics.

        Parameters
        ----------
        root : `Root`
            The root this node belongs to.
        inputs : dict [str, `Leaf`]
            Dictionary with the inputs used in the simulation.

            Key is the name of the corresponding input category e.g. 'bfield'.
        diagnostics : dict [str, `Dataset`]
            Dictionary with the diagnostics used in the simulation.

            Key is the name of the corresponding diagnostics e.g. 'inistate'.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)

        for category in metadata.input_categories:
            if category in inputs:
                self[category] = inputs[category]
            else:
                self[category] = None

        for diagnostic in metadata.simulation_diagnostics:
            if diagnostic in diagnostics:
                self[f"_{diagnostic}"] = diagnostics[diagnostic]
            else:
                self[f"_{diagnostic}"] = None

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        inputs, diagnostics = [], []
        for category in metadata.input_categories:
            try:
                _ = self[category]
                inputs.append(category)
            except AscotIOException:
                pass
        for diagnostic in metadata.simulation_diagnostics:
            try:
                _ = self[f"_{diagnostic}"]
                diagnostics.append(diagnostic)
            except AscotIOException:
                pass
        return (
            f"<{self.__class__.__name__}("
            f"parent={self._parent!r}, "
            f"inputs={inputs}, "
            f"diagnostics={diagnostics}, "
            f"frozen={self._frozen})>"
            )

    def __getattribute__(self, key: str) -> Any:
        """Return attribute unless it refers to an input or diagnostic not
        present in the simulation.

        Parameters
        ----------
        key : str
            Name of the attribute or input category or diagnostic.

        Returns
        -------
        value : Any
            Value of the attribute.

        Raises
        ------
        AscotIOException
            If the queried input or diagnostic was not used in the simulation.
        """
        value = super().__getattribute__(key)
        if key in metadata.input_categories and not value:
            raise AscotIOException(
                f"Input '{key}' was not used in the simulation."
                )
        if key.startswith("_") and len(key) > 1:
            diag_key = key[1:]
            if diag_key in metadata.simulation_diagnostics and not value:
                raise AscotIOException(
                    f"Diagnostics '{diag_key}' was not present in "
                    f"the simulation."
                )

        return value

    def _get_decorated_contents(self) -> str:
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        contents = ""
        contents += utils.decorate(f"{self.variant.ljust(10)} {self.qid}",
                                   bold=True)
        contents += f" {self.date}"
        contents += f"\n{self.note}\n\n"

        contents += utils.decorate("Diagnostics:\n", color="purple",
                                   underline=True, bold=True)
        for diagnostic in metadata.simulation_diagnostics:
            try:
                self[f"_{diagnostic}"]
            except AscotIOException:
                continue
            contents += f"- {diagnostic}\n"

        contents += utils.decorate("\nInputs:\n", color="purple",
                                   underline=True, bold=True)
        for category in metadata.input_categories:
            try:
                leaf = self[category]
            except AscotIOException:
                continue

            contents += utils.decorate(category.ljust(8), color="green")
            contents += utils.decorate(f"{leaf.variant.ljust(10)} {leaf.qid}",
                                       bold=True)
            contents += f" {leaf.date}"
            contents += f"\n{"".ljust(8)}{leaf.note}\n"

        return contents

    @property
    def contents(self) -> str:
        """A string representation of the contents.
        """
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self) -> None:
        """Show on screen the metadata of this run among with all the inputs
        that were used and the active diagnostics.
        """
        print(self._get_decorated_contents())

    def ls(self, show:bool=False) -> str:
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


class Root(ImmutableNodeHDF5):
    """The entry node for accessing data in the HDF5 file.

    The root node spawns the treeview creating all other nodes and data groups.

    Datasets are managed via this class since there is an interdependency
    between inputs and results: an input cannot be removed if it was used by
    a simulation before the associated result is removed.

    Attributes
    ----------
    _hdf5manager : `HDF5Manager`
        Manages the contents of the HDF5 file if one is provided.
    """

    def __init__(self, hdf5_filename=None, **kwargs: Any) -> None:
        """Initialize an empty tree structure.

        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(root=self, name="root", **kwargs)
        for category in metadata.input_categories:
            self[category] = InputCategory(self, name=category)
            self[category]._freeze()

        hdf5manager: Optional[HDF5Manager]
        hdf5manager = HDF5Manager(hdf5_filename) if hdf5_filename else None

        if hdf5manager:
            self._set_hdf5manager(hdf5manager)
            self._init_from_hdf5()

        self._freeze()

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(qids={self._qids}, tags={self._tags}, "
            f"active={self._active!r}, frozen={self._frozen}, "
            f"hdf5manager={self._hdf5manager!r})>"
            )

    def _init_from_hdf5(self):
        """Initializes the tree structure from the HDF5 file."""
        for category in metadata.input_categories:
            node = self[category]
            active, qids, variants = self._hdf5manager.read_node(node._name)
            for qid, variant in zip(qids, variants):
                date, note = self._hdf5manager.read_input(qid, variant)
                self._add_input_dataset(
                    MetaData(qid=qid, date=date, note=note, variant=variant),
                    store_hdf5=False,
                    )

            if active:
                leaf, _ = self._locate_leaf(active)
                self.activate_dataset(leaf)

        active, qids, variants = self._hdf5manager.read_node(self._name)
        for qid, variant in zip(qids, variants):
            date, note, inputqids = self._hdf5manager.read_output(qid, variant)
            self._add_simulation_output(
                MetaData(qid=qid, date=date, note=note, variant=variant),
                [],
                inputqids,
                store_hdf5=False,
                )
        if active:
            leaf, _ = self._locate_leaf(active)
            self.activate_dataset(leaf)

    def _add_input_dataset(
            self,
            meta: MetaData,
            note: str = None,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ):
        """Add new input dataset to this tree.

        This method uses the factory method to create the `Leaf` instance that
        is stored.

        Parameters
        ----------
        dataset : `Leaf`
            Dataset to be added.
        note : str, optional
            Note of the dataset (overrides original if given).
        dryrun : bool, optional
            Only check if a dataset with same QID exists on the tree, but
            otherwise do nothing.
        store_hdf5 : bool, optional
            Store the dataset in the HDF5 file (requires that the this tree
            has access to a file).

            If not explicitly set, the default behavior is to:
            - store the input to a file if filename was specified when this tree
              was created.
            - not to store the input if filename was not specified.
        """
        try:
            self._locate_leaf(meta.qid)
        except AscotIOException:
            pass
        else:
            raise AscotIOException(
                "There is already a dataset with identical QID."
                )
        leaf = self._leaf_factory(meta)
        if note:
            leaf.note = note
        if not dryrun:
            category = metadata.get_input_category(leaf.variant)
            self[category]._add_leaf(leaf)

            if store_hdf5 is None:
                store_hdf5 = self._hdf5manager is not None

            if store_hdf5:
                if not self._hdf5manager:
                    raise AscotIOException("No HDF5 file was provided.")
                self._hdf5manager.write_leaf(leaf)

        return leaf

    def _add_simulation_output(
            self,
            meta: MetaData,
            diagnostics: List [str],
            inputqids: List [str],
            note: str = None,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> None:
        """Add simulation output to this tree.

        dataset : `Leaf`
            Dataset to be added.
        note : str, optional
            Note on the dataset (overrides original if given).
        dryrun : bool, optional
            Only check if a dataset with same QID exists on the tree, but
            otherwise do nothing.
        store_hdf5 : bool, optional
            Store the dataset in the HDF5 file (requires that the this tree
            has access to a file).
        """
        try:
            self._locate_leaf(meta.qid)
        except AscotIOException:
            pass
        else:
            raise AscotIOException(
                "There is already a dataset with identical QID."
                )
        inputs = {}
        for qid in inputqids:
            inputleaf, _ = self._locate_leaf(qid)
            category = metadata.get_input_category(inputleaf.variant)
            inputs[category] = inputleaf

        leaf = self._leaf_factory(
            meta, inputs=inputs, diagnostics=diagnostics
            )
        if note:
            leaf.note = note

        if not dryrun:
            if store_hdf5 is None:
                store_hdf5 = self._hdf5manager is not None
            if store_hdf5 and not self._hdf5manager:
                raise AscotIOException("No HDF5 file was provided.")
            self._add_leaf(leaf, store_hdf5=store_hdf5)
        return leaf

    @classmethod
    def _leaf_factory(cls, meta, **kwargs):
        """Create `Leaf` instances of different variety.

        Override this method to create specific `Leaf` instances that are stored
        in this tree. This function is automatically called when using the
        `_add_input_dataset` and `_add_simulation_output` methods. This base
        implementation simply creates a generic `Leaf` instance.
        """
        if meta.variant in metadata.run_variants:
            return SimulationOutput(**meta._asdict(), **kwargs)
        return Leaf(**meta._asdict())

    def _get_decorated_contents(self) -> str:
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        contents = ""
        contents += utils.decorate("Inputs:", color="purple",
                                   underline=True, bold=True)
        contents += utils.decorate(" [only active shown]\n", color="green")
        for category in metadata.input_categories:
            contents += utils.decorate(f"{category.ljust(10)}", color="green")

            number_of_inputs  = len(self[category]._qids)
            if number_of_inputs == 0:
                contents += "*no inputs*\n\n"
                continue

            leaf = self[category].active
            contents += utils.decorate(
                f"{leaf.variant.ljust(10)}{leaf.qid}", bold=True,
                )
            contents += f" {leaf.date}"
            if number_of_inputs > 1:
                contents += f" + {number_of_inputs-1} other(s)"
            else:
                contents += " (no other inputs)"
            contents += f"\n{"".ljust(10)}{leaf.note}\n"

        contents += utils.decorate("\nSimulations:\n", color="purple",
                                   underline=True, bold=True)
        for leaf in self:
            contents += utils.decorate(
                f"{leaf.variant.ljust(10)}{leaf.qid}", bold=True,
                )
            contents += f" {leaf.date}"
            if leaf == self.active:
                contents += utils.decorate(" [active]", color="green")
            contents += f"\n{"".ljust(10)}{leaf.note}\n"
        if not self._qids:
            contents += "No simulation results.\n"
        return contents

    def _locate_leaf(self, qid: str) -> Tuple[Leaf, ImmutableNode]:
        """Find leaf and its parent corresponding to the given QID.

        Parameters
        ----------
        qid : str
            QID of the leaf.

        Returns
        -------
        leaf : `Leaf`
            Leaf corresponding to the given QID.
        parent : `ImmutableNode`
            The parent of the located leaf.

        Raises
        ------
        AscotIOException
            If the leaf does not belong to this tree.
        """
        qqid = metadata.get_qid(qid, with_prefix=True)
        if qqid in self:
            leaf = self[qqid]
            return leaf, self

        for category in metadata.input_categories:
            if qqid in self[category]:
                leaf = self[category][qqid]
                return leaf, self[category]

        raise AscotIOException(f"Leaf {qid} not found in the tree.")

    @property
    def contents(self) -> str:
        """A string representation of the contents.
        """
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

    def activate_dataset(self, dataset: Leaf) -> None:
        """Set group as active and rebuild the tree.

        Parameters
        ----------
        dataset : `Leaf` or str
            Dataset to be activated or its name or QID.
        """
        qid = metadata.get_qid(dataset)
        leaf, parent = self._locate_leaf(qid)
        parent._activate_leaf(leaf)

    def destroy_dataset(self, dataset: Leaf, repack: bool=True) -> None:
        """Remove dataset and associated data permanently.

        Parameters
        ----------
        dataset : `Leaf` or str
            Dataset to be removed or its name or QID.
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Without repacking only references to the data, not the actual data,
            are removed. Repacking has some overhead for large files and
            requires third-party `h5repack` tool.
        """
        qid = metadata.get_qid(dataset)
        leaf, _ = self._locate_leaf(qid)
        leaf.destroy(repack=repack)

    def destroy(self, repack: bool = True) -> None:
        """Remove all results and associated data permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Without repacking only references to the data, not the actual data,
            are removed. Repacking has some overhead for large files and
            requires third-party `h5repack` tool.
        """
        # This function overriden only to update the docstring.
        super().destroy(repack=repack)
