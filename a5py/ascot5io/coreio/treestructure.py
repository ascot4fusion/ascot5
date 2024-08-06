from __future__ import annotations
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
  which technically makes the structure a graph rather than a strict tree.

Note
----
`Leaf` is intended for internal use, while `MetaDataHolder` is for external
interfaces.
"""
import h5py
import random
import datetime
import warnings
from contextlib import contextmanager

import numpy as np

from a5py.exceptions import AscotIOException
from a5py import utils

default_tag = "TAG"
"""Default tag and description."""

input_categories = (
    "options", "bfield", "efield", "marker", "plasma", "neutral", "wall",
    "boozer", "mhd", "asigma", "nbi",
    )
"""All input categories."""

simulation_diagnostics = (
    "inistate", "endstate", "state", "orbits",
)
"""All simulation diagnostics."""

def generate_metadata():
    """Generate QID, date and default description/tag.

    Calls random number generator to create 32 bit string which is then
    converted as a QID string using left-padding with zeroes if necessary.

    Returns
    -------
    qid : str
        QID.
    data : str
        Date.
    desc : str
        Description.
    """
    qid = str( np.uint32( random.getrandbits(32) ) ).rjust(10, "0")
    date = utils.format2universaldate(datetime.datetime.now())
    description = default_tag

    return qid, date, description

def get_input_category(variant):
    """Return the input category corresponding to the given variant.

    Parameters
    ----------
    variant : str
        Variant of the data.

    Returns
    -------
    category : str
        Input category corresponding to the given variant.

    Raises
    ------
    ValueError
        If the variant is unknown.
    """
    if variant in ("B_TC", "B_GS", "B_2DS", "B_3DS", "B_3DST", "B_STS"):
        return "bfield"
    elif variant in ("E_TC", "E_1DS"):
        return "efield"
    elif variant in ("prt", "gc", "fl"):
        return "marker"
    elif variant in ("wall_2D", "wall_3D"):
        return "wall"
    elif variant in ("plasma_1D", "plasma_1DS", "plasma_1Dt"):
        return "plasma"
    elif variant in ("N0_1D", "N0_3D"):
        return "neutral"
    elif variant in ("Boozer"):
        return "boozer"
    elif variant in ("MHD_STAT", "MHD_NONSTAT"):
        return "mhd"
    elif variant in ("asigma_loc"):
        return "asigma"
    elif variant in ("NBI"):
        return "nbi"
    else:
        raise ValueError(f"Unknown variant: {variant}")


def get_qid(dataset, with_prefix=False):
    """Returns the QID from a given object or string.

    Parameters
    ----------
    dataset : `MetaDataHolder` or str
        Object that has QID, name of the object, or it's QID with or without
        the preceding 'q'.
    with_prefix : bool, optional
        Return the QID with the preceding 'q'.

    Returns
    -------
    qid : str
        A string with 10 digits.
    """
    if isinstance(dataset, str):
        if len(dataset) >= 10 and dataset[-10:].isdigit():
            qid = dataset[-10:]
        else:
            raise ValueError(
                f"This doesn't appear to be a valid QID: {dataset}"
                )
    else:
        qid = dataset.qid

    if with_prefix:
        qid = f"q{qid}"
    return qid


class MetaDataHolder():
    """Holds metadata and acts as a leaf for the tree data structure.

    Classes that have access to the actual data should inherit from this class.

    The metadata is immutable with the exception of `description`.

    Attributes
    ----------
    _root : `Root`
        Root node of the tree where this leaf belongs to.
    _ascot : `Ascot`
        The `Ascot` object for which the tree belongs to.
    _usedby : [MetaDataHolder]
        List of objects that reference this data.
    qid : str
        Unique identifier for this data.
    date : str
        Date when this data was created.
    description : str
        Short description for the user to document this data.
    variant : str
        What is the variant of the data this instance represents.
    """

    def __init__(self, qid, date, description, variant, **kwargs):
        """Initialize data container.
        """
        super().__init__(**kwargs)
        self._ascot = None
        self._parent = None

        self._qid = qid
        self._date = date
        self._usedby = []
        self._variant = variant
        self._description = description

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(root={self._parent!r}, "
            f"ascot={self._ascot!r}, qid={self._qid}, date={self._date}, "
            f"variant={self._variant})>"
            )

    def __setattr__(self, name, value):
        """Change the value of the attribute unless it is read-only.

        Parameters
        ----------
        name : str
            Name of the attribute.
        value : any
            Value to set.

        Raises
        ------
        AscotIOException
            If the attribute is read-only.
        """
        if name in ["qid", "qqid", "date", "variant", "name"]:
            raise AscotIOException(f"Attribute {name} is read only.")
        super().__setattr__(name, value)

    def _adopt(self, parent, ascot):
        """Adopt this leaf to a tree.

        Parameters
        ----------
        root : `Root`
            Root node of the tree where this leaf belongs to.
        ascot : `Ascot`
            The `Ascot` object for which the tree belongs to.

        Raises
        ------
        AscotIOException
            If the leaf is already adopted.
        """
        if self._parent:
            raise AscotIOException("Leaf is already adopted.")
        self._parent = parent
        self._ascot = ascot

    def _extract_tag(self):
        """Extracts a tag from description.

        Description is converted to a tag like this:

        1. The first word in the description is chosen, i.e., everything before
            the first whitespace.
        2. All special characters are removed from the first word.
        3. The word is converted to uppercase which becomes the tag.
        4. If the tag is invalid (empty string) or it starts with a number,
            the default tag is returned instead.

        Note that the returned value is not the actual tag used in the tree if
        there are other leafs with identical tags.

        Returns
        -------
        tag : str
            The tag.
        """
        tag_candidate = self.description.split(" ")[0]
        tag_candidate = "".join(ch for ch in tag_candidate if ch.isalnum())
        tag_candidate = tag_candidate.upper()
        if not len(tag_candidate) or tag_candidate[0] in "1234567890":
            return default_tag
        return tag_candidate

    @property
    def qid(self):
        """Unique identifier for this data."""
        return self._qid

    @property
    def qqid(self):
        """Unique identifier for this data with preceding 'q'."""
        return "q" + self._qid

    @property
    def date(self):
        """Date when this data was created."""
        return self._date

    @property
    def description(self):
        """Short description for the user to document this data."""
        return self._description

    @description.setter
    def description(self, description):
        """Set a short description to document this data."""
        self._description = description

    @property
    def variant(self):
        """What variant of data this object represents."""
        return self._variant

    @property
    def name(self):
        """Full name of the data."""
        return f"{self.variant}_{self.qid}"

    def activate(self):
        """Set this dataset as active.

        Active inputs are used when the simulation is run. Active datasets are
        also used during post-processing by default unless otherwise specified.
        """
        self._parent.activate_dataset(self.qid)

    def destroy(self, repack=True):
        """Remove this dataset.

        This also removes any data on the disk.

        Parameters
        ----------
        repack : bool, optional
            If True, repack the HDF5 file.

            Removing data from the HDF5 file only removes references to it and
            repacking is required to free the disk space. Repacking makes a copy
            of the HDF5 file and destroys the original, and it also requires
            3rd party tool `h5repack` which is why it's use is optional here.
        """
        self._parent._remove_leaf(self)
        self._parent = None
        self._ascot = None

Leaf = MetaDataHolder
"""Alias for `MetaDataHolder` for internal use.
"""

class ImmutableStorage():
    """Object which supports dictionary-like assignment and which can be made
    immutable.

    Attributes
    ----------
    frozen : bool
        Indicates whether the node is frozen, preventing attribute modification.
    """

    def __init__(self, **kwargs):
        """Initialize an empty node which is unfrozen.

        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        self._frozen = False
        super().__init__(**kwargs)

    def __repr__(self):
        """Return a string representation of this object."""
        return f"<{self.__class__.__name__}(frozen={self._frozen})>"

    def __setitem__(self, key, value):
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

    def __setattr__(self, key, value):
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

    def __getitem__(self, key):
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

    def _freeze(self):
        """Make this node immutable.
        """
        self._frozen = True

    def _unfreeze(self):
        """Make this node mutable.
        """
        self._frozen = False

    @contextmanager
    def _modify_attributes(self):
        """Open a context where attributes can be modified.
        """
        self._unfreeze()
        try:
            yield
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
       the user-given description.
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
    """

    def __init__(self, root, **kwargs):
        """Initialize an empty node which is unfrozen.

        Parameters
        ----------
        root : `Root`
            The root to which this node belongs.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)
        self._root = root
        self._qids = []
        self._tags = []
        self._active = None

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(root={self._root!r}, "
            f"qids={self._qids}, tags={self._tags}, active={self._active!r}, "
            f"frozen={self._frozen})>"
            )

    def __contains__(self, key):
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

    def __iter__(self):
        """Iterate over this node's leafs."""
        for qid in self._qids:
            yield self[f"q{qid}"]

    def _add_leaf(self, leaf):
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

    def _remove_leaf(self, leaf):
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
                if len(self._qids):
                    self._active = self[f"q{self._qids[0]}"]
                else:
                    self._active = None
            reference_by_name, reference_by_qid = leaf.name, leaf.qqid
            delattr(self, reference_by_qid)
            delattr(self, reference_by_name)

        self._organize()

    def _organize(self):
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
            if not len(self._qids):
                self._active = None
                for tag_to_be_removed in self._tags:
                    delattr(self, tag_to_be_removed)
                self._tags = []
                return

            activate_if_first_child()
            sort_qids_by_date()
            update_references_by_tag()

    @property
    def active(self):
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

    def destroy(self, repack=True):
        """Remove all datasets belonging to this node.

        Parameters
        ----------
        repack : bool, optional
            If True, repack the HDF5 file.

            Removing data from the HDF5 file only removes references to it and
            repacking is required to free the disk space. Repacking makes a copy
            of the HDF5 file and destroys the original, and it also requires
            3rd party tool `h5repack` which is why it's use is optional here.
        """
        for dataset in self:
            dataset.destroy(repack=repack)

class InputCategory(ImmutableNode):
    """Node that contains all inputs of the same category.
    """

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

            contents += f"\n{self._tags[index]}\n{leaf.description}\n\n"

        if not len(contents):
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

class SimulationOutput(ImmutableStorage, MetaDataHolder):
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

    def __init__(self, inputs, diagnostics, **kwargs):
        """Initialize simulation output node with given inputs and diagnostics.

        Parameters
        ----------
        root : `Root`
            The root this node belongs to.
        inputs : dict [str, `Dataset`]
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

        for category in input_categories:
            if category in inputs:
                self[category] = inputs[category]
            else:
                self[category] = None

        for diagnostic in simulation_diagnostics:
            if diagnostic in diagnostics:
                self[f"_{diagnostic}"] = diagnostics[diagnostic]
            else:
                self[f"_{diagnostic}"] = None

    def __repr__(self):
        """Return a string representation of this object."""
        inputs, diagnostics = [], []
        for category in input_categories:
            try:
                self[category]
                inputs.append(category)
            except AscotIOException:
                pass
        for diagnostic in simulation_diagnostics:
            try:
                self[f"_{diagnostic}"]
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

    def __getattribute__(self, key):
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
        if key in input_categories and not value:
            raise AscotIOException(
                f"Input '{key}' was not used in the simulation."
                )
        if key.startswith("_") and len(key) > 1:
            diag_key = key[1:]
            if diag_key in simulation_diagnostics and not value:
                raise AscotIOException(
                    f"Diagnostics '{diag_key}' was not present in "
                    f"the simulation."
                )

        return value

    def _get_decorated_contents(self):
        """Get a string representation of the contents decorated with ANSI
        escape sequences.
        """
        contents = ""
        contents += utils.decorate(f"{self.variant.ljust(10)} {self.qid}",
                                   bold=True)
        contents += f" {self.date}"
        contents += f"\n{self.description}\n\n"

        contents += utils.decorate("Diagnostics:\n", color="purple",
                                   underline=True, bold=True)
        for diagnostic in simulation_diagnostics:
            try:
                self[f"_{diagnostic}"]
            except AscotIOException:
                continue
            contents += f"- {diagnostic}\n"

        contents += utils.decorate("\nInputs:\n", color="purple",
                                   underline=True, bold=True)
        for category in input_categories:
            try:
                leaf = self[category]
            except AscotIOException:
                continue

            contents += utils.decorate(category.ljust(8), color="green")
            contents += utils.decorate(f"{leaf.variant.ljust(10)} {leaf.qid}",
                                       bold=True)
            contents += f" {leaf.date}"
            contents += f"\n{"".ljust(8)}{leaf.description}\n"

        return contents

    @property
    def contents(self):
        """A string representation of the contents.
        """
        return utils.undecorate(self._get_decorated_contents())

    def show_contents(self):
        """Show on screen the metadata of this run among with all the inputs
        that were used and the active diagnostics.
        """
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

class Root(ImmutableNode):
    """The entry node for accessing data in the HDF5 file.

    The root node spawns the treeview creating all other nodes and data groups.

    Datasets are managed via this class since there is an interdependency
    between inputs and results: an input cannot be removed if it was used by
    a simulation before the associated result is removed.

    Attributes
    ----------
    _hdf5_filename : `Ascot`
        The `Ascot` object whose HDF5 file is what this tree structure
        represents.
    """

    def __init__(self, hdf5_filename=None, **kwargs):
        """Initialize the tree structure based on the available data.

        Parameters
        ----------
        hdf5_filename : str
            The HDF5 file from which the tree is constructed and where new data
            will be written.

            If None, no data is read from the HDF5 file.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(root=self, **kwargs)
        self._hdf5_filename = hdf5_filename
        for category in input_categories:
            self[category] = InputCategory(self)

        self._freeze()

        #input_from_ids, output_from_ids = self._read_groups_from_ids()
        #input_from_hdf5, output_from_hdf5 = self._read_groups_from_hdf5()

        #input_groups = input_from_hdf5 + input_from_ids
        #output_groups = output_from_hdf5 + output_from_ids

    def __repr__(self):
        """Return a string representation of this object."""
        return (
            f"<{self.__class__.__name__}(qids={self._qids}, tags={self._tags}, "
            f"active={self._active!r}, frozen={self._frozen})>"
            )

    def _build(self):
        """(Re-)build node structure.
        """
        input_from_ids, output_from_ids = self._read_groups_from_ids()
        input_from_hdf5, output_from_hdf5 = self._read_groups_from_hdf5()
        input_from_memory, output_from_memory = self._read_groups_from_memory()

        hdf5_filename = self._hdf5_filename
        for v in list(vars(self)):
            delattr(self, v)

        super().__init__(root=self)
        for inp in fileapi.INPUTGROUPS:
            self[inp] = InputNode(self)
        self._hdf5_filename = hdf5_filename

        input_groups = input_from_hdf5 + input_from_ids + input_from_memory
        output_groups = output_from_hdf5 + output_from_ids + output_from_memory

            # if "results" in h5:
            #     results = h5["results"]
            #     for run in results.keys():

            #         # Fetch those input groups that correspond to this run.
            #         inputqids = fileapi.get_inputqids(h5, results[run])

            #         # Make a result group
            #         runnode = self._create_resultgroup(
            #             results[run].name, h5, inputqids)
            #         self._add_child(run, runnode)

            #     self._finalize(h5)

        self._freeze()

    def _read_groups_from_hdf5(self):
        """Read all input and output groups from a HDF5 file and return them as
        lists.

        Returns
        -------
        input_groups : list [TreeData]
            List of input groups.
        output_groups : list [TreeData]
            List of output groups.
        """
        input_groups, output_groups = [], []
        if self._hdf5_filename is None:
            return input_groups, output_groups

        with h5py.File(self._hdf5_filename, "r") as h5:
            for inp in h5.keys():
                if inp in fileapi.INPUTGROUPS:
                    qids = fileapi.get_qids(h5, inp)
                    for qid in qids:
                        g = fileapi.get_group(h5, qid)
                        date = fileapi.get_date(h5, g)
                        desc = fileapi.get_desc(h5, g)
                        address = _Address.from_hdf5(
                            hdf5_filename=self._hdf5_filename,
                            path_within_hdf5=g.name,
                            )
                        input_groups.append(TreeData())

        return input_groups, output_groups

    def _read_groups_from_ids(self):
        """Read all input and output groups from IDS and return them as lists.

        Returns
        -------
        input_groups : list [TreeData]
            List of input groups.
        output_groups : list [TreeData]
            List of output groups.
        """
        input_groups = []
        output_groups = []
        return input_groups, output_groups

    def _get_group(self, name):
        """Fetch group based on its QID or name.

        Parameters
        ----------
        group : str
            Name or QID of the group to be fetched.

        Returns
        -------
        group : DataGroup
            The group that was fetched.
        """
        qid = "q" + fileapi.get_qid(name)
        if qid in self:
            return self[qid]
        for parent in fileapi.INPUTGROUPS:
            if parent in self and qid in self[parent]:
                return self[parent][qid]

    def _add_input(
            self,
            dataset,
            ascot,
            description=None,
            dryrun=False,
            store_hdf5=None,
            ):
        """Add input dataset to this tree."""
        if description is not None:
            dataset.description = description
        if not dryrun:
            category = get_input_category(dataset.variant)
            dataset._adopt(self[category], ascot)
            self[category]._add_leaf(dataset)
            if store_hdf5:
                dataset._write_hdf5()

    def _add_run(
            self,
            dataset,
            ascot,
            description=None,
            store_hdf5=None,):
        """Add simulation output to this tree.
        """
        for category in input_categories:
            try:
                dataset[category]
            except AscotIOException:
                continue

            inputdata = dataset[category]
            if inputdata not in self[category]:
                raise AscotIOException(
                    f"Simulation {category} input cannot be found in the tree."
                    )

        if description is not None:
            dataset.description = description
        dataset._adopt(self, ascot)
        self._add_leaf(dataset)
        if store_hdf5:
            dataset._write_hdf5()

    def _locate_leaf(self, qid):
        """Find leaf and its parent corresponding to the given QID.

        Parameters
        ----------
        qid : str
            QID of the leaf.

        Returns
        -------
        leaf : `MetaDataHolder`
            Leaf corresponding to the given QID.
        parent : `ImmutableNode`
            The parent of the located leaf.

        Raises
        ------
        AscotIOException
            If the leaf does not belong to this tree.
        """
        qqid = get_qid(qid, with_prefix=True)
        if qqid in self:
            leaf = self[qqid]
            return leaf, self

        for category in input_categories:
            if qqid in self[category]:
                leaf = self[category][qqid]
                return leaf, self[category]

        raise AscotIOException(f"Leaf {qid} not found in the tree.")

    def destroy_dataset(self, dataset, repack=True):
        """Remove group and associated data permanently.

        Parameters
        ----------
        dataset : `MetaDataHolder` or str
            Dataset to be removed or its name or QID.
        repack : bool, optional
            If True, repack the HDF5 file.
        """
        qid = get_qid(dataset)
        leaf, parent = self._locate_leaf(qid)
        leaf.destroy(repack=repack)
        # with h5py.File(fn, "a") as f:
        #     fileapi.remove_group(f, group)

        # if repack:
        #     fntemp = fn + "_repack"
        #     subprocess.call(["h5repack", fn, fntemp], stdout=subprocess.DEVNULL)
        #     subprocess.call(["mv", fntemp, fn])

        # self._build(fn)

    def activate_dataset(self, group):
        """Set group as active and rebuild the tree.

        Parameters
        ----------
        group : str
            Name or QID of the group to be activated.
        """
        for inp in fileapi.INPUTGROUPS:
            if group in self[inp]:
                self[inp]._active = group
        #fn = self._ascot.file_getpath()
        #with h5py.File(fn, "a") as f:
        #    fileapi.set_active(f, group)


    def destroy(self, repack=True):
        """Remove all results from the HDF5 file.

        Parameters
        ----------
        repack : bool, optional
            If True, repack the HDF5 file.

            Removing data from the HDF5 file only removes references to it and
            repacking is required to free the disk space. Repacking makes a copy
            of the HDF5 file and destroys the original, and it also requires
            3rd party tool `h5repack` which is why it's use is optional here.
        """
        super().destroy(repack=repack)

    def ls(self, show=True):
        """Get a string representation of the contents.

        Parameters
        ----------
        show : str, optional
            If True, the contents are also printed on screen.

        Returns
        -------
        contents : str
            Multiline string decorated with ANSI escape sequences that list
            all results and their meta data and currently active inputs.
        """
        out = ""
        out += _FancyText.title("Inputs:")
        out += _FancyText.active(" [only active shown]\n")
        for inp in fileapi.INPUTGROUPS:
            if not inp in self:
                continue

            ngrp  = len(self[inp]._qids)
            qid   = self[inp].active.get_qid()
            date  = self[inp].active.get_date()
            desc  = self[inp].active.get_desc()
            gtype = self[inp].active.get_type()
            out += _FancyText.active(inp.ljust(8))
            out += _FancyText.header(gtype.ljust(10) + " " + qid)
            out += " " + date
            out += "\n" + desc + "\n"
            out += "+ " + str(ngrp-1) + " other(s)\n"

        out += _FancyText.title("Results:\n")
        for q in self._qids:
            date  = self["q"+q].get_date()
            desc  = self["q"+q].get_desc()
            gtype = self["q"+q].get_type()
            out += _FancyText.header(gtype.ljust(10) + " " + q)
            out += " " + date
            if q == self.active.get_qid():
                out += _FancyText.active(" [active]")

            out += "\n" + desc + "\n"

        if show:
            print(out)
        return out
