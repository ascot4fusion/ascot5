from __future__ import annotations
"""Module for building the treeview that shows ASCOT5 HDF5 file contents.

Treeview consists of nodes acting as containers for other nodes and data. The
tree is spanned by the `RootNode`.
"""
import subprocess
from contextlib import contextmanager
import h5py

from a5py.exceptions import AscotIOException

from . import fileapi
from .treedata import DataGroup

class _FancyText():
    """Helper class to decorate output of `ls` commands so it is easier to read.
    """
    _green     = "\033[32m"
    _purple    = "\033[35m"
    _reset     = "\033[0m"
    _bold      = "\033[01m"
    _underline = "\033[04m"

    @staticmethod
    def title(string):
        """Convert string to title (underlined bold purple text).

        Parameters
        ----------
        string : str
            String to be converted.

        Returns
        -------
        title : str
            The converted string.
        """
        return _FancyText._bold + _FancyText._underline + _FancyText._purple \
            + string + _FancyText._reset

    @staticmethod
    def header(string):
        """Convert string to header (bold text).

        Parameters
        ----------
        string : str
            String to be converted.

        Returns
        -------
        header : str
            The converted string.
        """
        return _FancyText._bold + string + _FancyText._reset

    @staticmethod
    def active(string):
        """Convert string to active input (green text).

        Parameters
        ----------
        string : str
            String to be converted.

        Returns
        -------
        active : str
            The converted string.
        """
        return _FancyText._green + string + _FancyText._reset

class _Address():
    """Data class that contains information on how to access the data.
    """

    from enum import Enum

    class Format(Enum):

        IMAS_IDS = 1
        HDF5_FILE = 2
        IN_MEMORY = 3

    def __init__(self,
                 hdf5_filename=None,
                 path_within_hdf5=None,
                 imas_ids=None,
                 pointer_to_c_struct=None,
                 **kwargs
                 ) -> None:
        super().__init__(**kwargs)
        self.imas_ids = imas_ids
        self.hdf5_filename = hdf5_filename
        self.path_within_hdf5 = path_within_hdf5
        self.pointer_to_c_struct = pointer_to_c_struct

        self.format = None
        inconsistent_address = False
        if imas_ids is not None:
            self.format = _Address.Format.IMAS_IDS
        if hdf5_filename is not None and path_within_hdf5 is not None:
            self.format = _Address.Format.HDF5_FILE
        if pointer_to_c_struct is not None:
            self.format = _Address.Format.IN_MEMORY
        no_address = self.format is None

        if inconsistent_address:
            raise ValueError("Inconsistent address specified.")
        if no_address:
            raise ValueError("No address specified.")


    @classmethod
    def from_hdf5(cls, hdf5_filename, path_within_hdf5) -> _Address:
        """Create an Address instance from an HDF5 file.

        Parameters
        ----------
        hdf5_filename : str
            Path to HDF5 file.

        Returns
        -------
        _Address
            An instance of _Address.
        """
        return cls(
            hdf5_filename=hdf5_filename,
            path_within_hdf5=path_within_hdf5,
            )

    @classmethod
    def from_imas(cls, imas_ids) -> _Address:
        """Create an Address instance from IMAS IDS.

        Parameters
        ----------
        imas_ids : list of str
            List of IMAS IDS.

        Returns
        -------
        _Address
            An instance of _Address.
        """
        return cls(imas_ids=imas_ids)

class TreeData():
    """Data container that also has meta data (QID, date, description).
    """

    def __init__(self, address, qid, date, description, type_, **kwargs):
        """Initialize data container.

        Parameters
        ----------
        address : :class:`_Address`
            Address of the data.
        qid : str
            Unique identifier for this data.
        date : str
            Date when this data was created.
        description : str
            Short description for the user to document this data.
        type_ : str
            What type of data this object represents.
        """
        self._qid = qid
        self._date = date
        self._type = type_
        self._address = address
        self._description = description

    def __setattr__(self, name, value):
        """Prevent changing read-only attributes.

        Parameters
        ----------
        name : str
            Name of the attribute.
        value : any
            Value to set.
        """
        if name in ["qid", "qqid", "date", "type", "name"]:
            raise AscotIOException(
                f"Attribute {name} is read only."
                )
        super().__setattr__(name, value)

    def _adopt(self, ascot, root):
        """Make this data as a part of the tree.

        Parameters
        ----------
        ascot : :class:`Ascot5`
            The ASCOT5 object this data container belongs to.
        root : :class:`RootNode`
            The root node this data container belongs to.
        """
        self._root = root
        self._ascot = ascot

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
        """Add short description to document this data."""
        self._description = description

    @property
    def type(self):
        """What type of data this object represents."""
        return self._type

    @property
    def name(self):
        """Name of the data"""
        return f"{self.type}_{self.qid}"

    def activate(self):
        """Set this group as active.

        Active inputs are used when the simulation is run. Active groups are
        also used during post-processing.
        """
        self._root.activate_group(self.qid)

    def destroy(self, repack=True):
        """Remove this group from the HDF5 file.

        Parameters
        ----------
        repack : bool, optional
            If True, repack the HDF5 file.

            Removing data from the HDF5 file only removes references to it and
            repacking is required to free the disk space. Repacking makes a copy
            of the HDF5 file and destroys the original, and it also requires
            3rd party tool `h5repack` which is why it's use is optional here.
        """
        self._root.destroy_group(self.qid, repack)

class _Node():
    """Base class which all tree nodes inherit.

    This base class provides all nodes with two functionalities:

    1. Its attributes can be accessed in dictionary-like manner, e.g.
       `node.child` and `node["child"]` are equivalent.
    2. Freeze (and unfreeze) this instance preventing (allowing) adding
       or removing attributes. Once the tree is constructed, all nodes
       should be frozen.

    Attributes
    ----------
    _frozen : bool
        When True, attributes cannot be added or removed for this node.
    """

    def __init__(self, **kwargs):
        """Initialize a mutable (unfrozen) node.

        Parameters
        ----------
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        self._frozen = False
        super().__init__(**kwargs)

    def _freeze(self):
        """Make this node immutable (frozen).
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

    def __contains__(self, key):
        """Check whether this node contains the attribute.

        Parameters
        ----------
        key : str
            Name of the attribute.

        Returns
        -------
        contains : bool
            True if this node contains the attribute.
        """
        return hasattr(self, key)

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

class _ParentNode(_Node):
    """Node that has data groups as children.

    This node stores QIDs of its children and provides a way of accessing them
    through active-flag, QID, name, and tag.

    Attributes
    ----------
    _qids : list [str]
        QIDs of this node's children.
    _tags : list [str]
        Unsorted list of all children's tags.
    _root : `RootNode`
        The `RootNode` this node belongs to.
    _active : `DataGroup`
        The currently active group.
    """

    def __init__(self, root, **kwargs):
        """Initialize a parent node that initially has no children.

        Parameters
        ----------
        root : `RootNode`
            The `RootNode` this node belongs to.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)
        self._qids  = []
        self._tags  = []
        self._root  = root
        self._active = None

    def __iter__(self):
        """Iterate over this node's children.

        Returns
        -------
        child : `DataGroup`
        """
        for qid in self._qids:
            yield self["q" + qid]

    def _add(self, group):
        """Add an input or a result group to this node.

        Parameters
        ----------
        group : `DataGroup`
            The child data.
        """
        self._qids.append(group.qid)

        with self._modify_attributes():
            reference_by_name, reference_by_qid = group.name, group.qqid
            self[reference_by_qid] = group
            self[reference_by_name] = group

        self._organize()

    def _remove(self, group):
        """Remove an input or a result group from this node.

        Parameters
        ----------
        group : `DataGroup`
            The child data.
        """
        self._qids.remove(group.qid)

        with self._modify_attributes():
            if self.active == group:
                if len(self._qids):
                    self._active = self["q" + self._qids[0]]
                else:
                    self._active = None
            reference_by_name, reference_by_qid = group.name, group.qqid
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
            """Activate the group if it is the first in this node."""
            try:
                self.active
            except AscotIOException as no_active_group:
                self._active = self["q" + self._qids[0]]

        def sort_qids_by_date():
            """Sort the collected list of qids by date."""
            dates = [self["q" + qid].date for qid in self._qids]
            self._qids = [
                qid for _, qid in sorted(zip(dates, self._qids), reverse=True)
                ]

        def update_references_by_tag():
            """Add references by tag with unique tags and remove the old ones.
            """
            for tag_to_be_removed in self._tags:
                delattr(self, tag_to_be_removed)
            self._tags = []

            dates = [self["q" + qid].date for qid in self._qids]
            descriptions = [self["q" + qid].description for qid in self._qids]
            unsorted_qids = self._qids
            unsorted_tags = [
                RootNode._description2tag(desc) for desc in descriptions
                ]

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
                self[new_tag] = self["q" + qid]

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

    @staticmethod
    def _description2tag(description):
        """Convert description to a valid tag.

        Description is converted to a tag like this:

        1. The first word in the description is chosen, i.e., everything before
           the first whitespace.
        2. All special characters are removed from the first word.
        3. The word is converted to uppercase which becomes the tag.
        4. If the tag is invalid (empty string) or it starts with a number,
           the default tag is returned instead.

        Parameters
        ----------
        description : str
            Description of the data.

        Returns
        -------
        tag : str
            The tag.
        """
        tag_candidate = description.split(" ")[0]
        tag_candidate = "".join(ch for ch in tag_candidate if ch.isalnum())
        tag_candidate = tag_candidate.upper()
        if not len(tag_candidate) or tag_candidate[0] in "1234567890":
            return "TAG"
        return tag_candidate

    @property
    def active(self):
        """Get the active group.

        Returns
        -------
        active : `DataGroup`
            The active group.

        Raises
        ------
        AscotIOException
            If there's no active group.
        """
        if self._active is None:
            raise AscotIOException(
                "No active group. Perhaps there's no data in this parent group?"
                )
        return self._active

    def destroy(self, repack=True):
        """Remove this group from the HDF5 file.

        Parameters
        ----------
        repack : bool, optional
            If True, repack the HDF5 file.

            Removing data from the HDF5 file only removes references to it and
            repacking is required to free the disk space. Repacking makes a copy
            of the HDF5 file and destroys the original, and it also requires
            3rd party tool `h5repack` which is why it's use is optional here.
        """
        self._root._destroy_group(self._path, repack=repack)

    def get_contents(self):
        """Return metadata for all childs sorted by date.

        Returns
        -------
        qids : array_like, str
            QIDs.
        dates : array_like, str
            Dates.
        descs : array_like, str
            Descriptions.
        types : array_like, str
            Data types.
        """
        nchild = len(self._qids)
        qids  = [None] * nchild
        dates = [None] * nchild
        descs = [None] * nchild
        types = [None] * nchild
        for i in range(nchild):
            qid = self._qids[i]
            qids[i]  = qid
            dates[i] = self["q"+qid].get_date()
            descs[i] = self["q"+qid].get_desc()
            types[i] = self["q"+qid].get_type()

        return qids, dates, descs, types

class InputNode(_ParentNode):
    """Node that represents an input parent group.

    Input parent groups are "bfield", "efield", etc. in the HDF5 file that can
    contain several input data groups.
    """

    def __init__(self, root, **kwargs):
        """Create an input node and its children.

        Parameters
        ----------
        root : `RootNode`
            The `RootNode` this node belongs to.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(root=root, **kwargs)
        #parent = h5[path]
        #for name, group in parent.items():
        #    inputtype = fileapi.get_type(name)
        #    inputobj  = root._create_datagroup(inputtype, group.name)
        #    self._add_child(name, inputobj)

        #self._finalize(h5)
        self._freeze()

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
            all data groups within this node.
        """
        out = ""
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

class ResultNode(_Node, DataGroup):
    """Node that represents a results group.

    Result groups contain the data groups containing simulation results. Result
    node has its own meta data (QID, etc.) and its children can only be
    referenced by name e.g. `resultnode.inistate`. Result node also stores QIDs
    of the input groups used to obtain that result.
    """

    def __init__(self, root, inputqids, **kwargs):
        """Create a result node and its children, and store references to
        the input data.

        Parameters
        ----------
        root : `RootNode`
            The `RootNode` this node belongs to.
        path : str
            Path to this node within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.
        inputqids : dict [str, str]
            Dictionary containing the name of the input parent group
            (e.g. "bfield") and the QID of the input used for this result.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(root=root, **kwargs)

        # Put references to the input data
        for name, qid in inputqids.items():
            self[name] = root[name]["q"+qid]

        # Store output data
        #for name, group in h5[path].items():
        #    self["_" + name] = root._create_datagroup(name, group.name)

        self._freeze()

    def get_contents(self):
        """Return names of all childs.

        Returns
        -------
        names : list [str]
            List of names.
        """
        out = []
        for i in fileapi.OUTPUTGROUPS:
            if i in self:
                out.append(i)

        return out

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
            this node's meta data, all output data within this node, and inputs
            that were used.
        """
        qid   = self.get_qid()
        date  = self.get_date()
        desc  = self.get_desc()
        gtype = self.get_type()
        out = ""
        out += _FancyText.header(gtype.ljust(10) + " " + qid)
        out += " " + date
        out += "\n" + desc + "\n"

        out += _FancyText.title("Contents:\n")
        for cnt in self.get_contents():
            out += cnt + "\n"

        out += _FancyText.title("Input:\n")
        for inp in fileapi.INPUTGROUPS:
            if not inp in self:
                continue
            qid   = self[inp].get_qid()
            date  = self[inp].get_date()
            desc  = self[inp].get_desc()
            gtype = self[inp].get_type()
            out += _FancyText.active(inp.ljust(8))
            out += _FancyText.header(gtype.ljust(10) + " " + qid)
            out += " " + date
            out += "\n" + desc + "\n"

        if show:
            print(out)
        return out

class RootNode(_ParentNode):
    """The entry node for accessing data in the HDF5 file.

    The root node spawns the treeview creating all other nodes and data groups.

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
        for inp in fileapi.INPUTGROUPS:
            self[inp] = InputNode(self)

        input_from_ids, output_from_ids = self._read_groups_from_ids()
        input_from_hdf5, output_from_hdf5 = self._read_groups_from_hdf5()

        input_groups = input_from_hdf5 + input_from_ids
        output_groups = output_from_hdf5 + output_from_ids

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

    def destroy_group(self, group, repack=True):
        """Remove group and associated data permanently.

        Parameters
        ----------
        group : str
            Name or QID of the group to be removed.
        repack : bool, optional
            If True, repack the HDF5 file.
        """
        group = f"q{group[-10:]}"
        for inp in fileapi.INPUTGROUPS:
            if group in self[inp]:
                self[inp]._remove(group)
        # with h5py.File(fn, "a") as f:
        #     fileapi.remove_group(f, group)

        # if repack:
        #     fntemp = fn + "_repack"
        #     subprocess.call(["h5repack", fn, fntemp], stdout=subprocess.DEVNULL)
        #     subprocess.call(["mv", fntemp, fn])

        # self._build(fn)

    def activate_group(self, group):
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

    def create_BTC(
            self,
            bxyz,
            jacobian,
            rhoval,
            description=None,
            activate=None,
            dryrun=False,
            store_hdf5=True,
            ):
        from a5py.ascot5io.bfield import B_TC
        obj = B_TC(bxyz, jacobian, rhoval)
        self.bfield._add(obj)
        obj._adopt(None, self)
        return obj

    def get_runs(self, inputqid):
        """Fetch QIDs of runs using given input.

        Parameters
        ----------
        inputqid : `qid`
            QID of the input.

        Returns
        -------
        qids : list [str]
            List of run QIDs.
        """
        # Find the parent group first
        parent  = None
        qid = "q" + inputqid
        for p in fileapi.INPUTGROUPS:
            if p in self and qid in self[p]:
                parent = p
        if parent is None:
            return []

        # Find the runs
        runqids = []
        for qid in self._qids:
            if parent in self["q"+qid] and \
               self["q"+qid][parent].get_qid() == inputqid:
                runqids.append(qid)

        return runqids

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
        self._destroy_group("results", repack=repack)

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
