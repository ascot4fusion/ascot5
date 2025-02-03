"""Module for building a tree structure for accessing ASCOT5 data.

The functionality of the tree is described and mostly implemented in
treeparts.py. Here we define the following classes:

- InputCategory: A node containing all input variants that belong to same
  category.
- RunVariant: Superclass for different run variants that contains the output
  data.
- Tree: The root node of the tree and the entry point.
"""
import warnings
from typing import Tuple, List, Dict, Any, Optional

from ...exceptions import AscotIOException
from ... import utils
from . import metadata
from .metadata import MetaData
from .treeparts import ImmutableNode, Leaf, TreeManager
from .immutablestorage import ImmutableStorage


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

    def destroy(self, repack: bool = True, **kwargs) -> None:
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


class RunVariant(ImmutableStorage, Leaf):
    """Leaf that contains data of a single simulation.

    Instances of this class contain all the metadata associated with the
    simulation. This class should be subclassed by various simulation variants
    to provide access to associated post-processing methods.

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
            inputs: Dict[str, Leaf],
            diagnostics: Dict[str, Any],
            **kwargs: Any,
            ) -> None:
        """Initialize simulation output node with given inputs and diagnostics.

        Parameters
        ----------
        inputs : dict [str, `Leaf`]
            Inputs used in this simulation by category.
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
                data = self[category]
                inputs.append(f"({category}:{data.qid})")
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
            f"qid={self.qid}, "
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
            contents += f"\n{''.ljust(8)}{leaf.note}\n"

        return contents

    @property
    def contents(self) -> str:
        """A string representation of the contents."""
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self) -> None:
        """Show on screen the metadata of this run among with all the inputs
        that were used and the active diagnostics.
        """
        print(self._get_decorated_contents())

    def ls(self, show: bool = False) -> str:
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


class Tree(ImmutableNode):
    """The entry node for accessing simulation data."""

    def __init__(self, hdf5file: Optional[Tuple[str, bool]] = None) -> None:
        """Initialize an empty tree structure.

        Parameters
        ----------
        hdf5file : Tuple[str, bool], optional
            Filename of the HDF5 file, if used, and flag indicating if the file
            exists.
        """
        super().__init__()
        nodes = {"root":self}
        for category in metadata.input_categories:
            self[category] = InputCategory()
            nodes[category] = self[category]
            self[category]._freeze()

        cls = type(self)
        TreeManager(cls._leaf_factory, hdf5file=hdf5file, **nodes)

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(qids={self._qids}, tags={self._tags}, "
            f"active={self._active!r}, frozen={self._frozen})>"
            )

    @classmethod
    def _leaf_factory(
        cls, meta: MetaData,
        inputs: Optional[Dict[str, Leaf]] = None,
        diagnostics: Optional[List] = None,
        **kwargs: Any,
        ) -> Leaf:
        """Create `Leaf` instances of different variety.

        Override this method to create specific `Leaf` instances that are stored
        in this tree. This function is automatically called when using the
        `_add_input_dataset` and `_add_simulation_output` methods. This base
        implementation simply creates a generic `Leaf` instance.
        """
        if(meta.variant in metadata.run_variants
            and (inputs is not None and diagnostics is not None) ):
            return RunVariant(
                **meta._asdict(), inputs=inputs, diagnostics=diagnostics,
                **kwargs,
            )
        return Leaf(**meta._asdict())

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
        for category in metadata.input_categories:
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

    def ls(self, show: bool = False) -> str:
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
            self, *, repack: bool = True, data: Leaf | str = "results", **kwargs,
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
        if data == "results":
            super().destroy(repack=repack, **kwargs)
        elif isinstance(data, str) and data in metadata.input_categories:
            self[data].destroy(repack=repack)
        else:
            qid = metadata.get_qid(data)
            if self._treemanager:
                leaf = self._treemanager.get_leaf(qid)
                self._treemanager.destroy_leaf(leaf, repack=repack)
