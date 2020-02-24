"""
Main module for reading ASCOT5 HDF5 files.

To use this module, initialize an Ascot object as
<pre>a5 = ascot5.Ascot("/path/to/ascot_hdf5_file.h5")</pre>
This object acts as an container object or Matlab-like struct, meaning one can
use it to inspect the Ascot HDF5 file e.g. as
<pre>a5.bfield.B_2D_1234567890</pre>
or, equivalently,
<pre>a5["bfield"]["B_2D_1234567890"]</pre>

One can also access the simulation results e.g. as
<pre>a5.run_1234567890["orbits"]</pre>

The lowest level objects in the hierarchy are objects that represents that
specific type of input or output, each with their own methods. These methods can
be used e.g. as
<pre>a5["bfield"]["B_2D_1234567890"].plot_fluxsurfaces()</pre>
or
<pre>a5.run_1234567890["orbits"].plot_2D("R", "z", endstate="wall")</pre>

The lowest level objects contain no data, but they can read the data that is
stored in the HDF5 file. To read the data into a Python dictionary, call
>a5.run_1234567890["inistate"].read()

The intermediate and top level objects (i.e. a5, a5.bfield, a5.run_1234567890 in
the examples above) are node-objects that are only used to navigate the HDF5
file. As such they only contain metadata and have no access to the HDF5 file
once the ascot5.Ascot object has been initialized.

The top level node contains input parent groups (bfield, efield, options,
neutral, marker, plasma, and wall) and the run groups that hold simulation
results.

There are several ways to use the ascot5.Ascot object to navigate. One can refer
to the active run as
<pre>a5.active</pre>
or active input field as
<pre>a5.bfield.active</pre>
To get the input field that was active in the given run:
<pre>a5.active.bfield</pre>
QID of the active field
<pre>a5.bfield.activeqid</pre>
Use QID as a reference to the field
<pre>a5.bfield.q1234567890</pre>

These examples also work with dictionary-like reference but here we use only the
attribute-like referencing for brevity.

You can even use field description to refer to it
<pre>a5.That_PRL_run</pre>
However, there are few rules to this:
- If the description is over 20 characters long, only the first 20 characters
  are used in referencing.
- Spaces and hyphens are turned to underscores and dots are removed. Avoid using
  any special characters.
- If two or more fields have identical descriptions (like the default _), there
  is no quarantee to which one the description refers to.

Finally, you can print the contents of a node with
<pre>a5.ls()</pre>
The list is ordered so that the first item is active qid and the rest are sorted
by date they were created from newest to oldest. You can use the index at which
the field appears in this list to reference it, so you can e.g. refer to the
active field as
<pre>a5.bfield[0]</pre>

Note: Methods and functions in this module cannot be used to modify the data in
the HDF5 file. You can (if you try hard enough) modify the object and its
attributes but then you are modifying only the object and not the HDF5 file. If
you do modify the attributes, then the functionality of this module is no longer
quaranteed. However, you can always re-initialize the object from the HDF5 file.

File: ascot5.py
"""

import h5py
import warnings

from . ascot5file import get_qid, get_activeqid, get_desc, get_date, get_type
from . ascot5file import get_inputqids

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

from a5py.ascot5io.state      import State
from a5py.ascot5io.orbits     import Orbits
from a5py.ascot5io.transcoef  import Transcoef
from a5py.ascot5io.dist_5D    import Dist_5D
from a5py.ascot5io.dist_6D    import Dist_6D
from a5py.ascot5io.dist_rho5D import Dist_rho5D
from a5py.ascot5io.dist_rho6D import Dist_rho6D

from a5py.ascot5io.ascot5file import INPUT_PARENTS


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

def create_inputobject(key, h5group):
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
        "prt" : mrk_prt, "gc" : mrk_gc, "fl" : mrk_fl,
        "wall_2D" : wall_2D, "wall_3D" : wall_3D,
        "plasma_1D" : plasma_1D, "plasma_1DS" : plasma_1DS,
        "N0_3D" : N0_3D,
        "Boozer" : Boozer, "MHD_STAT" : MHD, "MHD_NONSTAT" : MHD,
        "opt" : Opt
    }

    if key not in name_and_object:
        warnings.warn("Unknown input group " + key)
        return None

    return name_and_object[key](h5group)


def create_outputobject(key, h5group, runnode):
    """
    Create an output object based on the HDF5 group name

    Whenever you add a new output type, add it here and it then can be accessed
    via ASCOT object.
    """
    name_and_object = {
        "inistate" : State, "endstate" : State, "orbit" : Orbits,
        "dist5d" : Dist_5D, "dist6d" : Dist_6D, "distrho5d" : Dist_rho5D,
        "distrho6d" : Dist_rho6D, "transcoef" : Transcoef
    }

    if key not in name_and_object:
        warnings.warn("Unknown output group " + key)
        return None

    return name_and_object[key](h5group, runnode)


class _Node():
    """
    Class which lets its attributes be accessed in a dictionary-like manner.

    Instances of this class can be made (almost) immutable.
    """

    def __init__(self):
        """
        Initialize a mutable node.
        """
        self._frozen = False

    def _freeze(self):
        """
        Make this node immutable.
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

        Args:
            key: Name of the attribute
            value: Value of the attribute
        """
        if self._frozen:
            print("Frozen - new entries are not accepted.")
            return

        cleankey = self._remove_illegal_chars(key)
        setattr(self, cleankey, value)

    def __setattr__(self, key, value):
        """
        Add a new attribute this node.

        Args:
            key: Name of the attribute
            value: Value of the attribute
        """
        if key != "_frozen" and self._frozen:
            print("Frozen - new entries are not accepted.")

        else:
            cleankey = self._remove_illegal_chars(key)
            super().__setattr__(cleankey, value)

    def ls(self):
        """
        Print a string representation of this node.
        """
        print(str(self))

    def __contains__(self, key):
        try:
            getattr(self, key)
            return True
        except AttributeError:
            return False

    def __getitem__(self, key):
        """
        Allows accessing attributes dictionary-like.

        Args:
            key: Attribute name or index
        Returns:
            Attribute value or None if not found or invalid index
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

class _ContainerNode(_Node):
    """
    Node from which contents can be accessed via QID or description.
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

    def __getitem__(self, key):
        """
        Allows accessing attributes dictionary-like and by index.

        Args:
            key: Attribute name or index
        Returns:
            Attribute value or None if not found or invalid index
        """

        # Is item a direct reference or reference by index?
        if type(key) is str:
            # Direct reference
            return super().__getitem__(key)
        else:
            if key >= len(self._qids):
                print("Index out of bounds. Maximum index is "
                      + len(self._qids))
                return None

            return super().__getitem__(self._qids[key])

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

    def _init_store_activegroup(self, h5file, parent):
        self.activeqid = "q" + get_activeqid(h5file, parent)
        self.active    = self[self.activeqid]

    def _init_organize(self):
        # Organize qids, descriptions, dates and field names by active status
        # and date (active one first, then sorted by date from newest to oldest)
        index = self._qids.index(self.activeqid)

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
        AscotData object representing the given input data.
    """

    def __init__(self, parent):
        """
        Initialize this node by initializing input objects and storing them.
        """
        super().__init__()

        for key in parent.keys():
            type_ = get_type(parent[key].name.split("/")[-1])
            inputobj = create_inputobject(type_, parent[key])
            if inputobj is not None:
                self._init_store_qidgroup(parent.file, parent[key], inputobj)

        self._init_store_activegroup(parent.file, parent)

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

class _RunNode(_Node):
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

    def __init__(self, rungroup, inputgroups):
        super().__init__()

        self._qid  = get_qid(rungroup)
        self._date = get_date(rungroup.file, rungroup)
        self._desc = get_desc(rungroup.file, rungroup)

        # Put references to the input data
        for inp in range(0, len(inputgroups)):
            self[inputgroups[inp][0]] = inputgroups[inp][1]

        for key in rungroup:
            key = rungroup[key].name.split("/")[-1]
            outputobj = create_outputobject(key, rungroup[key], self)
            if outputobj is not None:
                self[key] = outputobj

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
            string += textcolor.active + inp.ljust(8) + textcolor.reset  \
                      + textcolor.header + self[inp].get_type().ljust(10) \
                      + " " + self[inp].get_qid() +  textcolor.reset    \
                      + "    " + self[inp].get_date()  \
                      + "\n        "                   \
                      + self[inp].get_desc() + "\n"

        string += textcolor.title + "\nOutput:\n" \
                  + textcolor.reset
        for key in vars(self):
            if key not in INPUT_PARENTS and key[0] != "_":
                string += key + "\n"

        return string

    def get_date(self):
        return self._date

    def get_qid(self):
        return self._qid

    def get_desc(self):
        return self._desc


class Ascot(_ContainerNode):
    """
    Top node used for exploring the HDF5 file.

    This node holds all input nodes and run nodes. Different run nodes (i.e.
    run groups) can be accessed by their QID or description whereas input nodes
    can only be accessed by their name e.g. bfield.
    """

    def __init__(self, fn):
        """
        Initialize the whole node structure recursively and create data objects.
        """
        super().__init__()
        self._hdf5fn = fn
        self.reload()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return

    def __str__(self):
        """
        Overview of inputs and results in the HDF5 file in a string format.
        """
        string = textcolor.title + "Inputs:\n" \
                 + textcolor.reset
        for inp in INPUT_PARENTS:
            if(hasattr(self, inp)):
                g = self[inp]
                string += textcolor.active + inp.ljust(8) + textcolor.reset \
                          + textcolor.header + g._types[0].ljust(10) + " "   \
                          + g._qids[0][1:] + textcolor.reset + " " \
                          + g._dates[0] \
                          + "\n        " + g._descs[0] + "\n"

        string += textcolor.title + "\nResults:\n" \
                  + textcolor.reset
        for i in range(0, len(self._qids)):
            string += textcolor.header + "run" + " " + self._qids[i][1:] + " " \
                      + textcolor.reset + self._dates[i] + "\n" + self._descs[i]
            if i < ( len(self._qids) - 1 ) :
                string += "\n"

        return string

    def reload(self):
        fn = self._hdf5fn
        for v in list(vars(self)):
            delattr(self, v)

        super().__init__()
        self._hdf5fn = fn

        with h5py.File(self._hdf5fn, "r") as h5:

            # Initialize input groups.
            for inp in h5.keys():
                if( inp in INPUT_PARENTS ):
                    self[inp] = _InputNode(h5[inp])

            if "results" in h5:
                for run in h5["results"].keys():

                    # Fetch those input groups that correspond to this run.
                    inputqids   = get_inputqids(h5["results"][run].file,
                                                h5["results"][run])

                    inputgroups = []
                    for inp in range(0, len(INPUT_PARENTS)):
                        if hasattr(self, INPUT_PARENTS[inp]):
                            groups = getattr(self, INPUT_PARENTS[inp])
                            inputgroups.append((INPUT_PARENTS[inp],
                                                groups["q" + inputqids[inp]]))

                    # Make a run node and store it.
                    runnode = _RunNode(h5["results"][run], inputgroups)
                    self._init_store_qidgroup(h5, h5["results"][run], runnode)

                self._init_store_activegroup(h5, h5["results"])
                self._init_organize()

        self._freeze()
