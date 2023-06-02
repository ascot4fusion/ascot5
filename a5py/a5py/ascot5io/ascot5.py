"""
Module for building treeview showing ASCOT5 HDF5 file contents.
"""
import h5py
import warnings
import subprocess

from collections import OrderedDict

from . ascot5file import get_qid, get_activeqid, get_desc, get_date, get_type,\
    set_desc, get_inputqids, remove_group, set_active

from a5py.exceptions import AscotFrozenException

from a5py.ascot5io.B_TC       import B_TC
from a5py.ascot5io.B_GS       import B_GS
from a5py.ascot5io.B_2DS      import B_2DS
from a5py.ascot5io.B_3DS      import B_3DS
from a5py.ascot5io.B_3DST     import B_3DST
from a5py.ascot5io.B_STS      import B_STS
from a5py.ascot5io.E_TC       import E_TC
from a5py.ascot5io.E_1DS      import E_1DS
from a5py.ascot5io.E_3D       import E_3D
from a5py.ascot5io.E_3DS      import E_3DS
from a5py.ascot5io.E_3DST     import E_3DST
from a5py.ascot5io.mrk_prt    import mrk_prt
from a5py.ascot5io.mrk_prt_shined    import mrk_prt_shined
from a5py.ascot5io.mrk_gc     import mrk_gc
from a5py.ascot5io.mrk_fl     import mrk_fl
from a5py.ascot5io.wall_2D    import wall_2D
from a5py.ascot5io.wall_3D    import wall_3D
from a5py.ascot5io.plasma_1D  import plasma_1D
from a5py.ascot5io.plasma_1DS import plasma_1DS
from a5py.ascot5io.N0_3D      import N0_3D
from a5py.ascot5io.boozer     import Boozer
from a5py.ascot5io.mhd        import MHD
from a5py.ascot5io.options    import Opt
from a5py.ascot5io.nbi        import nbi

from a5py.ascot5io.state      import State
from a5py.ascot5io.orbits     import Orbits
from a5py.ascot5io.transcoef  import Transcoef
from a5py.ascot5io.dist_5D    import Dist_5D
from a5py.ascot5io.dist_6D    import Dist_6D
from a5py.ascot5io.dist_rho5D import Dist_rho5D
from a5py.ascot5io.dist_rho6D import Dist_rho6D

from a5py.ascot5io.ascot5file import INPUT_PARENTS
from a5py.ascot5io.runmethods import RunMethods

class textcolor:
    """
    Colors for making the ls command easier to read.
    """
    reset='\033[0m'
    bold='\033[01m'
    underline='\033[04m'
    black='\033[30m'
    red='\033[31m'
    green='\033[32m'
    orange='\033[33m'
    blue='\033[34m'
    purple='\033[35m'
    cyan='\033[36m'
    lightgrey='\033[37m'
    darkgrey='\033[90m'
    lightred='\033[91m'
    lightgreen='\033[92m'
    yellow='\033[93m'
    lightblue='\033[94m'
    pink='\033[95m'
    lightcyan='\033[96m'
    title  = bold + underline + purple
    header = bold
    active = green

def create_inputobject(key, root, h5group):
    """
    Create an input object based on the HDF5 group name

    Whenever you add a new input type, add it here and it then can be accessed
    via ASCOT object.
    """
    name_and_object = {
        "B_TC" : B_TC, "B_GS" : B_GS, "B_2DS" : B_2DS, "B_3DS" : B_3DS,
        "B_3DST" : B_3DST, "B_STS" : B_STS,
        "E_TC" : E_TC, "E_1DS" : E_1DS, "E_3D" : E_3D, "E_3DS" : E_3DS,
        "E_3DST" : E_3DST,
        "prt" : mrk_prt, "prt_shined" : mrk_prt_shined, "gc" : mrk_gc,
        "fl" : mrk_fl,
        "wall_2D" : wall_2D, "wall_3D" : wall_3D,
        "plasma_1D" : plasma_1D, "plasma_1DS" : plasma_1DS,
        "N0_3D" : N0_3D,
        "Boozer" : Boozer, "MHD_STAT" : MHD, "MHD_NONSTAT" : MHD,
        "opt" : Opt,
        "nbi" : nbi
    }

    if key not in name_and_object:
        warnings.warn("Unknown input group " + key)
        return None

    return name_and_object[key](root, h5group)

class _ContainerNode(_Node):
    """
    Node that has data groups as children.

    This class is inherited by RootNode (which has run groups as child data
    groups) and _InputNode.
    """

    _MAX_DESC = 20

    def __init__(self):
        """
        Initialize a mutable node.
        """
        super().__init__()
        self._qids  = []
        self._descs = []
        self._dates = []
        self._types = []

    def _init_store_qidgroup(self, h5file, h5group, dataobject):
        groupname = h5group.name.split("/")[-1]
        qid = get_qid(groupname)
        self._qids.append("q" + qid)
        self._descs.append( get_desc(h5file, h5group) )
        self._dates.append( get_date(h5file, h5group) )
        self._types.append( get_type(h5group) )

        # Make sure reference by description is not too long
        max_ind = min( len(self._descs[-1]), self._MAX_DESC )

        referencename = self._remove_illegal_chars(groupname)
        descreference = self._descs[-1][:max_ind]

        self[referencename] = dataobject
        self["q" + qid]     = dataobject

        if descreference != "No description.":
            self[descreference] = dataobject

    def _init_organize(self):
        # Organize qids, descriptions, dates and field names by active status
        # and date (active one first, then sorted by date from newest to oldest)
        index = self._qids.index("q" + self.active.get_qid())

        sortedqids  = [ self._qids.pop(index)  ]
        sorteddates = [ self._dates.pop(index) ]
        sortedtypes = [ self._types.pop(index) ]
        sorteddescs = [ self._descs.pop(index) ]

        if len(self._qids) > 0:
            # This sorts elements in y by sorted x
            sortedqids  += \
                    [x for _, x in sorted(zip(self._dates, self._qids))]
            sortedtypes += \
                    [x for _, x in sorted(zip(self._dates, self._types))]
            sorteddescs += \
                    [x for _, x in sorted(zip(self._dates, self._descs))]
            sorteddates += sorted(self._dates)

        self._qids  = sortedqids
        self._dates = sorteddates
        self._types = sortedtypes
        self._descs = sorteddescs

    def destroy(self, repack=True):
        """
        Remove the group from the hdf5 file.
        """
        with h5py.File(self._file, "a") as f:
            remove_group(f, self._path)

        if repack:
            subprocess.call(["h5repack", self._file, "repack_" + self._file])
            subprocess.call(["mv", "repack_" + self._file, self._file])

class _InputNode(_ContainerNode):
    """
    Node that represents an input parent group.

    Create an instance that represents the given input data.

    This function is called for all input fields when ascot5.Ascot object is
    initialized. These objects should be light-weight when they are
    initialized and all actual input data reading should be left to be done
    once the methods of these objects are called.

    Args:
        type_: String from which the type of the input data is recognized
        h5pygroup: Input data's h5py group.

    Returns:
        DataGroup object representing the given input data.
    """

    def __init__(self, root, parent):
        """
        Initialize this node by initializing input objects and storing them.
        """
        super().__init__()
        self._root = root
        self._name = parent.name.split("/")[-1]

        for key in parent.keys():
            type_ = get_type(parent[key].name.split("/")[-1])
            inputobj = create_inputobject(type_, root, parent[key])
            if inputobj is not None:
                self._init_store_qidgroup(parent.file, parent[key], inputobj)

        self.active = self["q"+get_activeqid(parent.file, parent)]

        self._init_organize()

        self._freeze()

    def __str__(self):
        """
        Get a table showing fields, qids, dates, and descriptions as a string.
        """
        string = ""
        for i in range(0, len(self._qids)):
            string += textcolor.header + self._types[i].ljust(10) + " "  + \
                      self._qids[i][1:] + textcolor.reset + " "  +       \
                      self._dates[i]
            if i == 0:
                string += textcolor.active + "  [active]" + textcolor.reset

            string += "\n" + self._descs[i]
            if i < ( len(self._qids) - 1 ) :
                string += "\n"

        return string

    def get_inputinfo(self, sortbydate=False):
        """
        Return tuple (QIDs, types, descs, dates) for all childs.

        Results are either sorted by date or active first and then by date
        (latter is default).
        """
        qids  = self._qids.copy()
        types = self._types.copy()
        descs = self._descs.copy()
        dates = self._dates.copy()
        if sortbydate:
            if len(self._qids) > 0:
                sortedqids  = \
                    [x[1:] for _, x in sorted(zip(dates, qids))]
                qids = sortedqids
                sortedtypes = \
                    [x for _, x in sorted(zip(dates, types))]
                types = sortedtypes
                sorteddescs = \
                    [x for _, x in sorted(zip(dates, descs))]
                descs = sorteddescs
                sorteddates = sorted(dates)
                dates = sorteddates

        return (qids, types, descs, dates)

    def destroy(self, repack=True):
        """
        Remove group from the HDF5 file.
        """
        print(self._name)
        self._root._destroy(self._name, repack)

class _RunNode(_Node, RunMethods):
    """
    Create an instance that represents the given run group data.

    This function is called for all run groups when ascot5.Ascot object is
    initialized. These objects should be light-weight when they are initialized
    and all actual output data reading should be left to be done once the
    methods of these objects are called.

    This is different to ascot5._create_input_group() because that only
    initializes a single data object while this one initializes all objects
    (whose data is present in the run group) and returns a node containing
    them. The reason is that input fields have QIDs and other metadata
    while the groups within the run group doesn't.

    Args:
        rungroup : Run group's h5py group.
        inputgroups : tuple inputgroups' names and h5pygroups

    Returns:
        ascot5._StandardNode object representing the given run group data.
    """

    def __init__(self, root, rungroup, inputgroups):
        super().__init__()

        self._root = root
        self._qid  = get_qid(rungroup)
        self._date = get_date(rungroup.file, rungroup)
        self._desc = get_desc(rungroup.file, rungroup)

        # Put references to the input data
        for inp in range(0, len(inputgroups)):
            self[inputgroups[inp][0]] = inputgroups[inp][1]

        for key in rungroup:
            key = rungroup[key].name.split("/")[-1]
            if key == "inistate":
                self[key] = State(root, rungroup)
            if key == "endstate":
                self[key] = State(root, rungroup)
            if key == "orbit":
                self[key] = Orbits(root, rungroup)
            if key == "dist5d":
                self[key] = Dist_5D(root, rungroup)
            if key == "dist6d":
                self[key] = Dist_6D(root, rungroup)
            if key == "dist5drho":
                self[key] = Dist_rho5D(root, rungroup)
            if key == "dist6drho":
                self[key] = Dist_rho6D(root, rungroup)
            if key == "transcoef":
                self[key] = Transcoef(root, rungroup)

        # Store the filename and path in the file
        self._file = rungroup.file.filename
        self._path = rungroup.name
        self.init_runmethods()

        self._freeze()

    def __str__(self):
        string = textcolor.title + "Run:\n"  \
                 + textcolor.reset + textcolor.header \
                 + "        " + self._qid + textcolor.reset + " " \
                 + self._date + "\n" \
                 + "        " + self._desc + "\n"

        string += textcolor.title + "\nInput:\n" \
                  + textcolor.reset
        for inp in INPUT_PARENTS:
            if inp in self:
                string += textcolor.active + inp.ljust(14) + textcolor.reset  \
                          + textcolor.header + self[inp].get_type().ljust(10) \
                          + " " + self[inp].get_qid() +  textcolor.reset    \
                          + "    " + self[inp].get_date()  \
                          + "\n              "                   \
                          + self[inp].get_desc() + "\n"
            else:
                string += textcolor.red + inp.ljust(14) + textcolor.reset \
                          + textcolor.header + "*not present*\n" \
                          + textcolor.reset + "\n"

        string += textcolor.title + "\nOutput:\n" \
                  + textcolor.reset
        for key in vars(self):
            if key not in INPUT_PARENTS and key[0] != "_":
                string += key + "\n"

        return string

class RootNode(_ContainerNode):
    """
    Entry node for accessing data in the HDF5 file.

    Initializing this node builds rest of the tree. This object and its child
    nodes act as container objects for the data in the HDF5 file. At top level
    (this node) run nodes containing simulation results can be accessed as well
    as the parent nodes (bfield, efield, etc.) that in turn contain the actual
    input groups.

    The data can be accessed as
    <pre>root.bfield.B_2DS_1234567890</pre>
    or, equivalently,
    <pre>root["bfield"]["B_2DS_1234567890"]</pre>
    In each input group, one input is always set as "active" (meaning it would
    be used for the next simulation) and it can be accessed as
    <pre>root.bfield.active</pre>
    QID can be used as a reference as well
    <pre>root.bfield.q1234567890</pre>

    Run groups are accessed in a similar fashion, e.g.
    <pre>root.run_1234567890</pre>
    and the data within is accessed with
    <pre>root.run_1234567890["orbits"]</pre>
    However, the easiest way to access the simulation output is via the methods
    in the RunNode. The active run (by default the most recent simulation) can
    be accessed with
    <pre>root.active</pre>
    and its inputs as
    <pre>root.active.bfield</pre>

    Most of these examples also work with dictionary-like reference but here we
    use only the attribute-like referencing for brevity.

    You can even use field description to refer to it
    <pre>a5.That_PRL_run</pre>
    However, there are few rules to this:
    - If the description is over 20 characters long, only the first 20
      characters are used in referencing.
    - Spaces and hyphens are turned to underscores and dots are removed. Avoid
      using
      any special characters.
    - If two or more fields have identical descriptions (like the default _),
      there is no quarantee to which one the description refers to.

    Finally, you can print the contents of a node with
    <pre>a5.ls()</pre>
    The list is ordered so that the first item is active qid and the rest are
    sorted by date they were created from newest to oldest. You can use
    the index at which the field appears in this list to reference it, so you
    can e.g. refer to the active field as
    <pre>a5.bfield[0]</pre>

    Note: deleting or adding attributes to the nodes won't modify the contents
    of the HDF5 file. For that, use the methods found in this root node or its
    children.
    """

    def __init__(self, ascot):
        """
        Initialize the whole node structure recursively and create data objects.

        Args:
          ascot : Ascot
            Ascot object this node belongs to.
        """
        super().__init__()
        self._ascot = ascot

        fn = self._ascot.file_getpath()
        if fn is not None:
            h5py.File(fn, "r") # Try opening the file
            self._build(fn)

    def __str__(self):
        """
        Overview of inputs and results in the HDF5 file in a string format.
        """
        string = textcolor.title + "Inputs:\n" + textcolor.reset
        for inp in INPUT_PARENTS:
            if(hasattr(self, inp)):
                g = self[inp]
                string += textcolor.active + inp.ljust(8) + textcolor.reset \
                          + textcolor.header + g._types[0].ljust(10) + " "   \
                          + g._qids[0][1:] + textcolor.reset + " " \
                          + g._dates[0] \
                          + "\n        " + g._descs[0] + "\n"

        string += textcolor.title + "\nResults:\n" + textcolor.reset
        for i in range(0, len(self._qids)):
            string += textcolor.header + "run" + " " + self._qids[i][1:] + " " \
                      + textcolor.reset + self._dates[i] + "\n" + self._descs[i]
            if i < ( len(self._qids) - 1 ) :
                string += "\n"

        return string

    def _build(self, fn):
        """
        (Re-)build node structure from file.
        """
        ascot = self._ascot
        for v in list(vars(self)):
            delattr(self, v)

        super().__init__()
        self._ascot = ascot
        if fn is None:
            return

        with h5py.File(fn, "r") as h5:

            # Initialize input groups.
            for inp in h5.keys():
                if( inp in INPUT_PARENTS ):
                    self[inp] = _InputNode(self, h5[inp])

            if "results" in h5:
                for run in h5["results"].keys():

                    # Fetch those input groups that correspond to this run.
                    inputqids = get_inputqids(
                        h5["results"][run].file, h5["results"][run],
                        ignore = ["nbi","marker_shined"])

                    inputgroups = []
                    for inp in range(len(inputqids)):

                        if hasattr(self, INPUT_PARENTS[inp]):
                            groups = getattr(self, INPUT_PARENTS[inp])
                            inputgroups.append((INPUT_PARENTS[inp],
                                                groups["q" + inputqids[inp]]))

                    # Make a run node and store it.
                    runnode = _RunNode(self, h5["results"][run], inputgroups)
                    self._init_store_qidgroup(h5, h5["results"][run], runnode)

                self.active = self["q"+get_activeqid(h5, h5["results"])]
                self._init_organize()

        self._freeze()

    def _remove_group(self, group, repack=True):
        """
        Remove group from file.
        """
        fn = self._ascot.file_getpath()
        with h5py.File(fn, "a") as f:
            remove_group(f, group)

        if repack:
            fntemp = fn + "_repack"
            subprocess.call(["h5repack", fn, fntemp])
            subprocess.call(["mv", fntemp, fn])

        self._build()

    def _activate_group(self, group):
        """
        Set group as active.
        """
        fn = self._ascot.file_getpath()
        with h5py.File(fn, "a") as f:
            set_active(f, group)

        self._build()

    def create_input(inputtype, inputdata):
        """
        Create input and write the data to the HDF5 file.

        Args:
          inputtype : str
            Type of the input e.g. "B_2DS" or "options".
          inputdata : dict
            Dictionary containing all the data that is needed to create the
            requested input type.
        """
        self._build()

    def create_premade(write=True, **kwargs):
        """
        Create inputs based on predefined simulations and write the data.
        """
        
        self._build()

    def destroy_runs(self):
        """
        Remove every run from the HDF5 file.
        """
        self._remove_from_file("results", repack=True)

    def get_runs(self, inputgroup=None):
        """
        Fetch QIDs of all runs.

        Args:
          inputgroup :
        """
        # Find the parent group
        for parent in INPUT_PARENTS:
            if parent in self and "q" + inputqid in self[parent]:
                break

        runqids = []
        for qid in self._qids:
            if parent in self[qid] and self[qid][parent].get_qid() == inputqid:
                runqids.append(qid[1:])

        return runqids

    def get_parents(self):
        """
        Return ordered Dict of the parent groups.

        Dict contains name of the group and the corresponding node or None if
        the node is not present in the file.
        """
        parents = OrderedDict()
        with h5py.File(self._hdf5fn, "r") as h5:
            for p in INPUT_PARENTS:
                if p in h5:
                    parents[p] = self[p]
                else:
                    parents[p] = None

        return parents

    def get_resultsinfo(self, sortbydate=False):
        """
        Return tuple (QIDs, types, descs, dates) for all childs.

        Results are either sorted by date or active first and then by date
        (latter is default).
        """
        qids  = self._qids.copy()
        types = self._types.copy()
        descs = self._descs.copy()
        dates = self._dates.copy()
        if sortbydate:
            if len(self._qids) > 0:
                sortedqids  = \
                    [x[1:] for _, x in sorted(zip(dates, qids))]
                qids = sortedqids
                sortedtypes = \
                    [x for _, x in sorted(zip(dates, types))]
                types = sortedtypes
                sorteddescs = \
                    [x for _, x in sorted(zip(dates, descs))]
                descs = sorteddescs
                sorteddates = \
                    [x[:-1] for x in sorted(dates)]
                dates = sorteddates

        return (qids, types, descs, dates)

class _MetaDataMixin():
    """
    Mixin class for all objects that contain meta data (qid, date, type
    """

class _Node():
    """
    Base class which all tree nodes inherit.

    This base class provides all nodes with two functionalities:

    1. Its attributes can be accessed in dictionary-like manner, e.g.
       node.child and node["child"] are equivalent.
    2. Freeze (and unfreeze) this instance preventing (allowing) adding
       or removing attributes. Once the tree is constructed, all nodes
       should be frozen.

    Attributes:
      _frozen : bool
        When True, attributes cannot be added or removed for this node.
    """

    def __init__(self):
        """
        Initialize a mutable (unfrozen) node.
        """
        self._frozen = False

    def _freeze(self):
        """
        Make this node immutable (frozen).
        """
        self._frozen = True

    def _unfreeze(self):
        """
        Make this node mutable.
        """
        self._frozen = False

    def __setitem__(self, key, value):
        """
        Add a new attribute this node in dictionary style.
        """
        if self._frozen:
            raise AscotFrozenException()

        cleankey = self._remove_illegal_chars(key)
        setattr(self, cleankey, value)

    def __setattr__(self, key, value):
        """
        Add a new attribute this node.

        Args:
          key: str
            Name of the attribute.
          value: any
            Value of the attribute.
        """
        if key != "_frozen" and self._frozen:
            raise AscotFrozenException()
        else:
            cleankey = self._remove_illegal_chars(key)
            super().__setattr__(cleankey, value)

    def __contains__(self, key):
        """
        Called when quering if key in node.

        Args:
          key : str
            Name of the attribute.
        Returns:
          bool
            True if this node contains the attribute.
        """
        try:
            getattr(self, key)
            return True
        except AttributeError:
            return False

    def __getitem__(self, key):
        """
        Called when accessing attributes in dictionary-like manner.

        Args:
          key : str
            Name of the attribute.
        Returns:
          any
            Value of the attribute.
        """
        return getattr(self, key)

    @staticmethod
    def _remove_illegal_chars(key):
        """
        Remove illegal characters from argument so it becomes a valid attribute

        Args:
            String to be cleaned
        Returns:
            String from which illegal characters are replaced or removed
        """
        key = key.replace(" ", "_")
        key = key.replace("-", "_")
        key = key.replace(".","")
        return key
