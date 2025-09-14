"""Class defition for `ImmutableStorage`."""
from __future__ import annotations

import warnings
from contextlib import contextmanager
from typing import Any, Generator, List, Optional, Union, TYPE_CHECKING

from a5py import utils
from a5py.exceptions import AscotIOException

from .leaf import Leaf, Status
from .hdf5 import HDF5MiniManager

if TYPE_CHECKING:
    from .tree import TreeManager

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
        """Make this node immutable."""
        self._frozen = True

    def _unfreeze(self) -> None:
        """Make this node mutable."""
        self._frozen = False

    @contextmanager
    def _modify_attributes(self) -> Generator[ImmutableStorage, None, None]:
        """Open a context where attributes can be modified."""
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
        if not isinstance(key, str):
            key = key.qqid
        return hasattr(self, key) or hasattr(self, f"q{key}")

    def __iter__(self) -> Generator[Leaf, None, None]:
        """Iterate over this node's leafs."""
        for qid in self._qids:
            yield self[f"q{qid}"]

    def _activate_leaf(self, leaf) -> None:
        """Set given leaf as active.

        No reorganization required as only single reference is updated.

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

        Adds the QID to the QID list, and includes references by QID and
        name. Then reorganizes this node to update any other references.

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

        Removes the QID from the QID list, and removes references by QID and
        name. Then reorganizes this node to update any other references.

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
        def clear_references():
            """Remove all references."""
            self._active = None
            for tag_to_be_removed in self._tags:
                if not tag_to_be_removed is None:
                    delattr(self, tag_to_be_removed)
            self._tags = []

        def activate_if_first_leaf():
            """Activate the leaf if it is the first one in this node or if the
            previously active leaf was removed from this node."""
            try:
                self.active
            except AscotIOException:
                self._active = self[f"q{self._qids[0]}"]
            finally:
                if not self.active in self:
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
                if not tag_to_be_removed is None:
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
                if not new_tag is None:
                    self[new_tag] = self[f"q{qid}"]

        with self._modify_attributes():
            isempty = len(self._qids) == 0
            if isempty:
                clear_references()
            else:
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

    def destroy(self, *, repack: bool=False, **kwargs) -> None:
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


class InputCategory(ImmutableNode):
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

    def ls(self, show: bool=False):
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

    def destroy(self, repack: bool=True, **kwargs) -> None:
        """Remove all inputs and data within this category permanently.

        Parameters
        ----------
        repack : bool, optional
            Repack the HDF5 file reducing the size of the file on disk.

            Repacking has some overhead but without it only the references to
            the data, not the data itself, are removed in the file.
        """
        # This function overriden only to update the docstring.
        super().destroy(repack=repack, **kwargs)


@Leaf.register("output")
class OutputLeaf(ImmutableStorage, Leaf):
    """Leaf that contains data of a single simulation.

    Instances of this class contain all the metadata associated with the
    simulation. This class should be subclassed by various simulation variants
    to provide access to associated post-processing methods.

    References to the input datasets used in the simulation are stored by
    category. For instance, `node.bfield` is the magnetic field dataset used in
    the simulation. It should point to the same object which can be found within
    the 'bfield' category in the tree.

    If trying to access an input not used in the simulation, an exception is
    raised.
    """

    def __init__(
            self,
            qid: str,
            date: str,
            note: str,
            variant: str,
            inputs: dict[str, Leaf],
            **kwargs: Any,
            ) -> None:
        """Initialize simulation output node with given inputs.

        Parameters
        ----------
        qid : str
            Unique identifier for this data.
        date : str
            Date when this data was created.
        note : str
            Short note for the user to document this data.
        variant : str
            What data variant this instance represents.
        inputs : dict [str, `Leaf`]
            Inputs used in this simulation by category.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(
            qid=qid, date=date, note=note, variant=variant, **kwargs
            )
        self._inputs = []
        self._diagnostics = {}
        for category, leaf in inputs.items():
            self[category] = leaf
            self._inputs.append(category)

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        inputs = []
        for category in self._inputs:
            data = self[category]
            inputs.append(f"({category}:{data.qid})")
        return (
            f"<{self.__class__.__name__}("
            f"qid={self.qid}, "
            f"inputs={inputs}, "
            f"saved={not self._file is None})>"
            )

    def __getattribute__(self, key: str) -> Any:
        """Return attribute unless it refers to an input not present in the
        simulation.

        Parameters
        ----------
        key : str
            Name of the attribute or input category.

        Returns
        -------
        value : Any
            Value of the attribute.

        Raises
        ------
        AscotIOException
            If the queried input was not used in the simulation.
        """
        try:
            value = super().__getattribute__(key)
        except AttributeError:
            value = None
        if value is None and not key.startswith("_"):
            raise AscotIOException(
                f"Input '{key}' was not used in the simulation."
                )
        return value

    def _setup(self, params):
        """Create diagnostics and make them ready for simulation.

        Subclasses should override this method.
        """
        _ = params
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement '_setup'"
        )

    def _load(self, file):
        """Setup diagnostics from file.

        Subclasses should extend this method to create diagnostics based on the
        data on file. This implementation checks that the diagnostics are not
        already set and sets the filemanager.
        """
        if self.status is Status.SAVED:
            raise AscotIOException("Cannot load twice.")
        if len(self._diagnostics):
            raise AscotIOException("Diagnostics setup already.")

        self._file = file

    def _get_decorated_contents(self) -> str:
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        contents = ""
        contents += utils.decorate(f"{self.variant.ljust(10)} {self.qid}",
                                   bold=True)
        contents += f" {self.date}"
        contents += f"\n{self.note}\n"

        contents += utils.decorate("\nInputs:\n", color="purple",
                                   underline=True, bold=True)
        for category in self._inputs:
            leaf = self[category]

            contents += utils.decorate(category.ljust(8), color="green")
            contents += utils.decorate(f"{leaf.variant.ljust(10)} {leaf.qid}",
                                       bold=True)
            contents += f" {leaf.date}"
            contents += f"\n{''.ljust(8)}{leaf.note}\n"

        return contents

    @property
    def contents(self) -> str:
        """A string representation of the contents."""
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self) -> None:
        """Show on screen the metadata of this run among with all the inputs
        that were used.
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
