"""Module for building the treeview that shows ASCOT5 HDF5 file contents.

Treeview consists of nodes acting as containers for other nodes and data. The
tree is spanned by the `RootNode`.
"""
import subprocess
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
                "Node is frozen and attributes cannot be set")

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
                "Node is frozen and attributes cannot be set")
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
    _root : `RootNode`
        The `RootNode` this node belongs to.
    _path : str
        Path to this node within the HDF5 file.
    active : `DataGroup`
        The currently active group.
    """

    def __init__(self, root, path, **kwargs):
        """Initialize a parent node that initially has no children.

        Parameters
        ----------
        root : `RootNode`
            The `RootNode` this node belongs to.
        path : str
            Path to this node within the HDF5 file.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(**kwargs)
        self._qids  = []
        self._root  = root
        self._path  = path
        self.active = None

    def _add_child(self, key, datagroup):
        """Add a child to this node.

        This should be called only when the HDF5 file is open and the tree is
        being built. After all childs have been added, `_finalize` must be
        called to organize this node.

        Parameters
        ----------
        key : str
            Name of the child in format <type>_<QID>.
        datagroup : `DataGroup`
            The child data.
        """
        qid = fileapi.get_qid(key)
        self._qids.append(qid)

        # Add reference by name
        self[key] = datagroup

        # Add reference by QID
        self["q" + qid] = datagroup

    def _finalize(self, h5):
        """Organize contents and add references.

        This method should be called once after all childs have been added and
        the HDF5 file is still open.

        Parameters
        ----------
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.
        """
        # Find active group and set it
        qid = fileapi.get_activeqid(h5, h5[self._path])
        self.active = self["q"+qid]

        # List all dates but remove the active
        dates = []
        for q in self._qids:
            grp = fileapi.get_group(h5, q)
            dates.append(fileapi.get_date(h5, grp))

        # Sort QIDs by date starting from most recent
        self._qids = [x for _, x in sorted(zip(dates, self._qids))][::-1]

        # Add references by tag
        def desc2tag(desc):
            """Convert desc to tag.
            """
            # Cut from first whitespace
            desc = desc.split(" ")[0]

            # Maximum length is 10 characters
            #if len(desc) > 10:
            #    desc = desc[:10]

            # Remove all special characters
            desc = "".join(c for c in desc if c.isalnum())

            # Make all caps
            desc = desc.upper()

            # Use default if empty first character is number
            if desc == "" or desc[0] in "1234567890":
                desc = "TAG"

            return desc

        replacedtags = [] # Tags that appear more than once are deleted
        for q in self._qids:
            tag = desc2tag(fileapi.get_desc(h5, q))
            if tag in self:
                # This tag is repeated

                if tag not in replacedtags:
                    # First time we notice a tag is repeated, rename the tag
                    # already present to "tag_0" and add the tag to list of tags
                    # we remove afterwards. (We cannot remove the tag here as
                    # then we wouldn't notice repeated entries).
                    replacedtags.append(tag)
                    self[tag + "_0"] = self[tag]

                # Find next available index and create "tag_i".
                i = 0
                tag0 = tag
                while tag in self:
                    tag = tag0 + "_" + str(i)
                    i += 1

            self[tag] = self["q"+q]

        # Remove repeated tags as they have their renamed variants in place.
        for tag in replacedtags:
            delattr(self, tag)

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

    def __init__(self, root, path, h5, **kwargs):
        """Create an input node and its children.

        Parameters
        ----------
        root : `RootNode`
            The `RootNode` this node belongs to.
        path : str
            Path to this node within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(root=root, path=path, **kwargs)
        parent = h5[path]
        for name, group in parent.items():
            inputtype = fileapi.get_type(name)
            inputobj  = root._create_datagroup(inputtype, group.name)
            self._add_child(name, inputobj)

        self._finalize(h5)
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

    def __init__(self, root, path, h5, inputqids, **kwargs):
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
        super().__init__(root=root, path=path, **kwargs)

        # Put references to the input data
        for name, qid in inputqids.items():
            self[name] = root[name]["q"+qid]

        # Store output data
        for name, group in h5[path].items():
            self["_" + name] = root._create_datagroup(name, group.name)

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
    _ascot : `Ascot`
        The `Ascot` object whose HDF5 file is what this tree structure
        represents.
    """

    def __init__(self, ascot, **kwargs):
        """Initialize the tree structure based on the HDF5 file.

        This opens the HDF5 file obtained from `Ascot.file_getpath` and builds
        the tree structure. Nothing is done if the filename is `None`

        Parameters
        ----------
        ascot : `Ascot`
            The `Ascot` object whose HDF5 file is what this tree structure
            represents.
        **kwargs
            Arguments passed to other constructors in case of multiple
            inheritance.
        """
        super().__init__(root=self, path="/", **kwargs)
        self._ascot = ascot

        fn = self._ascot.file_getpath()
        if fn is not None:
            h5py.File(fn, "r") # Try opening the file
            self._build(fn)

    def _build(self, fn):
        """(Re-)build node structure from file.

        Use this liberally to update the treeview every time the contents have
        changed.

        Parameters
        ----------
        fn : str
            The HDF5 file from which the tree is constructed.
        """
        ascot = self._ascot
        for v in list(vars(self)):
            delattr(self, v)

        super().__init__(root=self, path="/results")
        self._ascot = ascot
        if fn is None:
            return

        with h5py.File(fn, "r") as h5:

            # Initialize input groups.
            for inp in h5.keys():
                if inp in fileapi.INPUTGROUPS:
                    self[inp] = self._create_inputgroup(inp, h5)

            if "results" in h5:
                results = h5["results"]
                for run in results.keys():

                    # Fetch those input groups that correspond to this run.
                    inputqids = fileapi.get_inputqids(h5, results[run])

                    # Make a result group
                    runnode = self._create_resultgroup(
                        results[run].name, h5, inputqids)
                    self._add_child(run, runnode)

                self._finalize(h5)

        self._freeze()

    def _destroy_group(self, group, repack=True):
        """Remove a group from the HDF5 file and rebuild the tree.

        Parameters
        ----------
        group : str
            Name of the group to be removed.
        repack : bool, optional
            If True, repack the HDF5 file.
        """
        fn = self._ascot.file_getpath()
        with h5py.File(fn, "a") as f:
            fileapi.remove_group(f, group)

        if repack:
            fntemp = fn + "_repack"
            subprocess.call(["h5repack", fn, fntemp], stdout=subprocess.DEVNULL)
            subprocess.call(["mv", fntemp, fn])

        self._build(fn)

    def _activate_group(self, group):
        """Set group as active and rebuild the tree.

        Parameters
        ----------
        group : str
            Name or QID of the group to be activated.
        """
        fn = self._ascot.file_getpath()
        with h5py.File(fn, "a") as f:
            fileapi.set_active(f, group)

        self._build(fn)

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

    def _create_inputgroup(self, path, h5):
        """Create an input group to be added to the treeview.

        Parameters
        ----------
        path : str
            Path to the input group within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.

        Returns
        -------
        group : `InputNode`
            The input group that was created
        """
        pass

    def _create_resultgroup(self, path, h5, inputgroups):
        """Create a result group to be added to the treeview.

        Parameters
        ----------
        path : str
            Path to the result node within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.
        inputqids : dict [str, str]
            Dictionary containing the name of the input parent group
            (e.g. "bfield") and the QID of the input used for this result.

        Returns
        -------
        group : `ResultNode`
            The result group that was created.
        """
        pass

    def _create_datagroup(self, grouptype, group):
        """Create data group based on the given type.

        Parameters
        ----------
        grouptype : str
            Type of the group as it appears in `HDF5TOOBJ`.
        path : str
            Path to the data in the HDF5 file that corresponds to the group.

        Returns
        -------
        `DataContainer`
            The data group that was created.
        """
        pass

    @staticmethod
    def _create_file(fn):
        """Create an empty HDF5 file.
        """
        f = h5py.File(fn, "w")
        f.close()

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
