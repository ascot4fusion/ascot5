"""Contains classes that represent nodes of the tree data structure."""
from __future__ import annotations

from typing import Generator, Optional, Union, Any

from a5py import utils
from a5py.exceptions import AscotDataException, AscotMeltdownError

from .leaves import Leaf
from .tree import TreeManager, ROOT


class ImmutableNode(utils.ImmutableStorage):
    """Tree node which can store other nodes or leaves, and whose attributes
    cannot be altered once frozen.

    This class provides the following main functionalities:

    1. Attributes can be accessed in a dictionary-like manner, e.g.,
       ``node.child`` and ``node["child"]`` are equivalent.
    2. Freezing this instance prevents setting and removing any attributes.
    3. Leaves can be accessed using their name or a tag constructed from
       the user-given note.
    4. One leaf is always set as 'active' and is accessed via the ``active``
       property.
    5. Leaves within this node can be iterated and they are ordered by their
       date of creation.
    6. Whether a leaf belongs to this node can be checked with ``leaf in node``,
       where ``leaf`` is :class:`.Leaf`, its name, or
       its tag.

    The attributes can be modified using a context:

    .. code-block:: python

       with node._modify_attributes():
           node.new_attribute = "new_value"

    Initially node is empty and mutable.
    """

    def __init__(self, **kwargs: Any) -> None:
        """
        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)
        self._names: list[str] = []
        """Names of this node's leaves sorted by date starting from newest."""

        self._tags: list[str] = []
        """List of all tags that can be used to access the leaves.

        Note that not every leaf has a tag.
        """

        self._active: Optional[Leaf] = None
        """The currently active leaf."""

        self._treemanager: Optional[TreeManager] = None
        """The manager of the tree this node belongs to.

        Other than in testing environment, the nodes are always part of a tree.
        """

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(names={self._names}, "
            f"tags={self._tags}, active={self._active!r}, "
            f"frozen={self._frozen})>"
            )

    def __len__(self) -> int:
        """Return the number of leaves in this node."""
        return len(self._names)

    def __contains__(self, key: Union[str, Leaf]) -> bool:
        """Check whether this node contains the requested leaf.

        Parameters
        ----------
        key : str or :class:`.Leaf`
            The requested leaf or it's QID.

        Returns
        -------
        contains : bool
            True if this node contains the leaf.
        """
        if not isinstance(key, str):
            key = key.name
        return hasattr(self, key)

    def __iter__(self) -> Generator[Leaf, None, None]:
        """Iterate over this node's leafs."""
        for name in self._names:
            yield self[name]

    @property
    def active(self) -> Leaf:
        """The active data.

        Raises
        ------
        :class:`.AscotDataException`
            If there's no active data.
        """
        if self._active is None:
            raise AscotDataException(
                "No active data. Perhaps this node is empty?"
                )
        return self._active

    def _activate_leaf(self, leaf: Leaf) -> None:
        """Set given leaf as active.

        No reorganization required as only single reference is updated.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf to be set as active.

        Raises
        ------
        :class:`.AscotMeltdownError`
            If the leaf does not belong to this node.
        """
        if leaf not in self:
            raise AscotMeltdownError(
                f"Variant '{leaf.name}' does not belong to this node."
                )
        with self._modify_attributes():
            self._active = leaf

    def _add_leaf(self, leaf: Leaf) -> None:
        """Add leaf to this node.

        Adds the QID to the QID list, and includes references by QID and
        name. Then reorganizes this node to update any other references.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf to be added.

        Raises
        ------
        :class:`.AscotMeltdownError`
            If the leaf already belongs to this node or this node has another
            leaf with the same name.
        """
        if leaf in self:
            raise AscotMeltdownError(
                f"Variant '{leaf.name}' already belongs to this node."
                )
        if leaf.name in self._names:
            raise AscotMeltdownError(
                f"This node already has a variant named '{leaf.name}'."
                )
        self._names.append(leaf.name)
        with self._modify_attributes():
            self[leaf.name] = leaf

        self._organize()

    def _remove_leaf(self, leaf: Leaf) -> None:
        """Remove a leaf from this node.

        Removes the QID from the QID list, and removes references by QID and
        name. Then reorganizes this node to update any other references.

        Parameters
        ----------
        leaf : :class:`.Leaf`
            The leaf to be removed.

        Raises
        ------
        :class:`.AscotMeltdownError`
            If the leaf does not belong to this node.
        """
        if leaf not in self:
            raise AscotMeltdownError(
                f"Variant '{leaf.name}' does not belong to this node."
                )
        self._names.remove(leaf.name)
        with self._modify_attributes():
            delattr(self, leaf.name)

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
        def clear_references() -> None:
            """Remove all references."""
            self._active = None
            for tag_to_be_removed in self._tags:
                if tag_to_be_removed is not None:
                    delattr(self, tag_to_be_removed)
            self._tags = []

        def activate_if_first_leaf() -> None:
            """Activate the leaf if it is the first one in this node or if the
            previously active leaf was removed from this node.
            """
            try:
                self.active
            except AscotDataException:
                self._active = self[self._names[0]]
            finally:
                if self.active not in self:
                    self._active = self[self._names[0]]

        def sort_names_by_date() -> None:
            """Sort the collected list of qids by date."""
            dates = [leaf.date for leaf in self]
            self._names = [
                name for _, name in sorted(zip(dates, self._names), reverse=True)
                ]

        def update_references_by_tag() -> None:
            """Add references by tag with unique tags and remove the old ones.
            """
            for tag_to_be_removed in self._tags:
                if tag_to_be_removed is not None:
                    delattr(self, tag_to_be_removed)
            self._tags = []

            dates = [leaf.date for leaf in self]
            unsorted_names = self._names
            unsorted_tags = [Leaf.extract_tag(leaf.note)[0] for leaf in self]

            counts: dict[str, int] = {}
            dates.reverse()
            for tag, _, name in sorted(zip(unsorted_tags, dates, unsorted_names)):
                if tag is None:
                    continue
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
                if new_tag is not None:
                    self[new_tag] = self[name]

        with self._modify_attributes():
            isempty = len(self._names) == 0
            if isempty:
                clear_references()
            else:
                activate_if_first_leaf()
                sort_names_by_date()
                update_references_by_tag()

    def destroy(self, *, repack: bool=False, **kwargs: Any) -> None:
        """Remove all data belonging to this node permanently.

        Parameters
        ----------
        repack : bool, *optional*
            Repack the HDF5 file reducing the size of the file on disk.

            Repacking has some overhead but without it only the references to
            the data, not the data itself, are removed in the file.
        """
        _ = kwargs # kwargs are required only for subclasses to have kw args
        if self._treemanager:
            for name in list(self._names):
                leaf = self[name]
                # Repacking takes time so repack only at the last item
                if len(self._names) > 1:
                    self._treemanager.destroy_leaf(leaf, repack=False)
                else:
                    self._treemanager.destroy_leaf(leaf, repack=repack)


class InputCategory(ImmutableNode):
    """Node that contains all inputs of the same category."""

    def _get_decorated_contents(self) -> str:
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        contents = ""
        for index, leaf in enumerate(self):
            contents += utils.decorate(f"{leaf.name.ljust(15)}", bold=True)
            contents += f" {leaf.date}"
            if self.active == leaf:
                contents += utils.decorate(" [active]", color="green")

            contents += f"\n{self._tags[index]}\n{leaf.note}\n\n"

        if not contents:
            contents = "No data in this category.\n"
        return contents

    @property
    def contents(self) -> str:
        """A string representation of the contents."""
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self) -> None:
        """Show on screen the metadata of all inputs within this category and
        which input is active."""
        print(self._get_decorated_contents())

    def destroy(self, repack: bool=True, **kwargs: Any) -> None:
        """Remove all inputs and data within this category permanently.

        Parameters
        ----------
        repack : bool, *optional*
            Repack the HDF5 file reducing the size of the file on disk.

            Repacking has some overhead but without it only the references to
            the data, not the data itself, are removed in the file.
        """
        # This function overriden only to update the docstring.
        super().destroy(repack=repack, **kwargs)


class Tree(ImmutableNode):
    """The entry node for accessing simulation data."""

    def __init__(
            self, input_categories: list[str],
            hdf5file: Optional[tuple[str, bool]]=None
            ) -> None:
        """Initialize an empty tree structure.

        Parameters
        ----------
        hdf5file : Tuple[str, bool], *optional*
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

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(names={self._names}, tags={self._tags}, "
            f"active={self._active!r})>"
            )

    @property
    def contents(self) -> str:
        """A string representation of the contents."""
        return utils.undecorate(self._get_decorated_contents())

    def activate(self, data: Leaf | str) -> None:
        """Mark data as active.

        Parameters
        ----------
        data : :class:`.Leaf` or str
            Data to be activated (or it's QID or name).
        """
        if self._treemanager:
            if isinstance(data, str):
                data = self._treemanager.get_leaf(data)
            self._treemanager.activate_leaf(data)

    def destroy(
            self, *, repack: bool=True, data: Leaf | str=ROOT, **kwargs: Any,
            ) -> None:
        """Remove data permanently.

        Parameters
        ----------
        repack : bool, *optional*
            Repack the HDF5 file reducing the size of the file on disk.

            Repacking has some overhead but without it only the references to
            the data, not the data itself, are removed in the file.
        data : :class:`.Leaf` or str or None
            Data to be removed (or it's QID or name or entire category).

            The default value removes all results.
        """
        if data == ROOT:
            super().destroy(repack=repack, **kwargs)
        elif isinstance(data, str) and data in self:
            self[data].destroy(repack=repack)
        else:
            if self._treemanager:
                if isinstance(data, str):
                    data = self._treemanager.get_leaf(data)
                else:
                    data = self._treemanager.get_leaf(data.name)
                self._treemanager.destroy_leaf(data, repack=repack)

    def show_contents(self) -> None:
        """Show on screen the metadata of the currently active inputs and list
        of all simulation runs.
        """
        print(self._get_decorated_contents())

    def _get_decorated_contents(self) -> str:
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        def print_category(category: str) -> str:
            return utils.decorate(f"{category.ljust(10)}", color="green")

        def print_howmanyinputs(number_of_inputs: int) -> str:
            if number_of_inputs > 1:
                return f" + {number_of_inputs-1} other(s)"
            if number_of_inputs == 1:
                return " (no other inputs)"
            return "*no inputs*\n\n"

        def print_name(leaf: Leaf) -> str:
            return (
                utils.decorate(
                    f"{leaf.name.ljust(15)}", bold=True,
                )
                + f" {leaf.date}"
            )

        def print_note(leaf: Leaf) -> str:
            return f"\n{''.ljust(15)}\"{leaf.note}\"\n"

        def print_title(title: str) -> str:
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
            except AscotDataException:
                pass

            n_inp = len(self[category])
            if leaf is None:
                contents += ""
            else:
                contents += print_name(leaf)
            contents += print_howmanyinputs(n_inp)
            if leaf is not None:
                contents += print_note(leaf)

        contents += print_title("\nSimulations:\n")
        for leaf in self:
            contents += print_name(leaf)
            if leaf == self.active:
                contents += utils.decorate(" [active]", color="green")
            contents += print_note(leaf)

        if not self._names:
            contents += "No simulation results.\n"
        return contents
