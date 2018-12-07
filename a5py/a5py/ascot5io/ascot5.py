"""
Main module for reading ASCOT5 HDF5 files.

To use this module, initialize an Ascot object as
>a5 = ascot5.Ascot("/path/to/ascot_hdf5_file.h5")

This object acts as an container object or Matlab-like struct, meaning one can
use it to inspect the Ascot HDF5 file e.g. as
>a5.bfield.B_2D_1234567890

or, equivalently,
>a5["bfield"]["B_2D_1234567890"]

One can also access the simulation results e.g. as
>a5.run_1234567890["orbits"]

The lowest level objects in the hierarchy are objects that represents that
specific type of input or output, each with their own methods. These methods can
be used e.g. as
>a5["bfield"]["B_2D_1234567890"].plot_fluxsurfaces()

or
>a5.run_1234567890["orbits"].plot_2D("R", "z", endstate="wall")

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
>a5.active

or active input field as
>a5.bfield.active

These examples also work with dictionary-like reference but here we use only the
attribute-like referencing for brevity.

To get the input field that was active in the given run:
>a5.active.bfield

QID of the active field
>a5.bfield.activeqid

Use QID as a reference to the field
>a5.bfield.q1234567890

List of QIDs of all fields
>a5.bfield.qids

The list is ordered so that the first item is active qid and the rest are sorted
by date they were created from newest to oldest. You can use the index at which
the field appears in this list to reference it, so you can e.g. refer to the
active field as
>a5.bfield[0]

You can even use field description to refer to it
>a5.That_PRL_run

However, there are few rules to this:
- If the description is over 20 characters long, only the first 20 characters
  are used in referencing.
- Spaces and hyphens are turned to underscores and dots are removed. Avoid using
  any special characters.
- If two or more fields have identical descriptions (like the default _), there
  is no quarantee to which one the description refers to.

Finally, you can print the contents of a node with
>a5.ls

Note: Methods and functions in this module cannot be used to modify the data in
the HDF5 file. You can (if you try hard enough) modify the object and its
attributes but then you are modifying only the object and not the HDF5 file. If
you do modify the attributes, then the functionality of this module is no longer
quaranteed. However, you can always re-initialize the object from the HDF5 file.

File: ascot5.py
"""

import h5py

from . ascot5file import get_qid, get_activeqid, get_desc, get_date, get_type
from . ascot5file import get_inputqids

from a5py.ascot5io.B_TC       import B_TC
from a5py.ascot5io.B_GS       import B_GS
from a5py.ascot5io.B_2DS      import B_2DS
from a5py.ascot5io.B_3DS      import B_3DS
from a5py.ascot5io.E_TC       import E_TC
from a5py.ascot5io.mrk_prt    import mrk_prt
from a5py.ascot5io.mrk_gc     import mrk_gc
from a5py.ascot5io.mrk_fl     import mrk_fl
from a5py.ascot5io.wall_2D    import wall_2D
from a5py.ascot5io.wall_3D    import wall_3D
from a5py.ascot5io.plasma_1D  import plasma_1D
from a5py.ascot5io.N0_3D      import N0_3D
from a5py.ascot5io.options    import options

from a5py.ascot5io.state      import State
from a5py.ascot5io.orbits     import Orbits
from a5py.ascot5io.dist_5D    import Dist_5D
from a5py.ascot5io.dist_6D    import Dist_6D
from a5py.ascot5io.dist_rho5D import Dist_rho5D
from a5py.ascot5io.dist_rho6D import Dist_rho6D

def _create_input_group(type_, h5pygroup):
    """
    Create an instance that represents the given input data.

    This function is called for all input fields when ascot5.Ascot object is
    initialized. These objects should be light-weight when they are initialized
    and all actual input data reading should be left to be done once the methods
    of these objects are called.

    Args:
        type_: String from which the type of the input data is recognized
        h5pygroup: Input data's h5py group.

    Returns:
        AscotData object representing the given input data.
    """

    # We simply determine the type of the input data, initialize corresponding
    # object and return.
    inputobj = None
    if type_ == "B_TC":
        inputobj = B_TC(h5pygroup)

    if type_ == "B_GS":
        inputobj = B_GS(h5pygroup)

    if type_ == "B_2DS":
        inputobj = B_2DS(h5pygroup)

    if type_ == "B_3DS":
        inputobj = B_3DS(h5pygroup)

    if type_ == "E_TC":
        inputobj = E_TC(h5pygroup)

    if type_ == "particle":
        inputobj = mrk_prt(h5pygroup)

    if type_ == "guiding_center":
        inputobj = mrk_gc(h5pygroup)

    if type_ == "field_line":
        inputobj = mrk_fl(h5pygroup)

    if type_ == "wall_2D":
        inputobj = wall_2D(h5pygroup)

    if type_ == "wall_3D":
        inputobj = wall_3D(h5pygroup)

    if type_ == "plasma_1D":
        inputobj = plasma_1D(h5pygroup)

    if type_ == "N0_3D":
        inputobj = N0_3D(h5pygroup)

    if type_ == "options":
        inputobj = options(h5pygroup)

    return inputobj

def _create_run_group(h5pygroup):
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
        h5pygroup: Run group's h5py group.

    Returns:
        ascot5._StandardNode object representing the given run group data.
    """

    # Make a new node instance and unfreeze it so that we can put some data in
    rungroup = _StandardNode()
    rungroup._unfreeze()

    # rungroup is a node where each node is different type of output data object
    for key in h5pygroup:
        if key == "inistate":
            rungroup[key] = State(h5pygroup["inistate"])

        if key == "endstate":
            rungroup[key] = State(h5pygroup["endstate"])

        if key == "orbits":
            rungroup[key] = Orbits(h5pygroup["orbits"])
        if key == "dists":
            for d in h5pygroup["dists"]:
                if d == "dist_5D":
                    rungroup[d] = Dist_5D(h5pygroup[key][d])
                if d == "dist_6D":
                    rungroup[d] = Dist_5D(h5pygroup[key][d])
                if d == "dist_rho5D":
                    rungroup[d] = Dist_rho5D(h5pygroup[key][d])
                if d == "dist_rho6D":
                    rungroup[d] = Dist_rho6D(h5pygroup[key][d])


    # Put references to the input data, these qids are replaced with actual
    # references where this function was called from.
    qids = get_inputqids(h5pygroup.file, h5pygroup)
    rungroup["options"] = "q" + qids[0]
    rungroup["bfield"]  = "q" + qids[1]
    rungroup["efield"]  = "q" + qids[2]
    rungroup["marker"]  = "q" + qids[3]
    rungroup["plasma"]  = "q" + qids[4]
    rungroup["neutral"] = "q" + qids[5]
    rungroup["wall"]    = "q" + qids[6]

    # Freeze and return
    rungroup._freeze()
    return rungroup

class _StandardNode():
    """
    Class which lets its attributes be accessed in a dictionary-like manner.

    Instances of this class can be made (almost) immutable.
    """

    ## How many characters are used in description based reference
    _MAX_DESC = 20

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
            print("Frozen. No new entries can be added manually.")
            return

        cleankey =_remove_illegal_chars(key)
        setattr(self, cleankey, value)

    def __setattr__(self, key, value):
        """
        Add a new attribute this node.

        Args:
            key: Name of the attribute
            value: Value of the attribute
        """
        if key != "_frozen" and self._frozen:
            print("Frozen. No new entries can be added manually.")

        else:
            cleankey =_remove_illegal_chars(key)
            super().__setattr__(cleankey, value)

    def __getitem__(self, key):
        """
        Access attributes in a dictionary style.

        Args:
            key: Name of the attribute to be fetched
        """
        return getattr(self, key)

    def ls(self):
        """
        Print a string representation of this node.
        """
        print(str(self))

class _InputNode(_StandardNode):
    """
    Node for accessing input data.
    """

    def __init__(self, parent):
        """
        Initialize this node by initializing input objects and storing them.
        """
        super().__init__()
        self.qids  = []
        self.descs = []
        self.dates = []
        self.types = []

        for key in parent.keys():
            qid = get_qid(key)
            self.qids.append("q" + qid)
            self.descs.append( get_desc(parent.file, parent[key]) )
            self.dates.append( get_date(parent.file, parent[key]) )
            self.types.append( get_type(key) )

            cleankey =_remove_illegal_chars(key)
            self[cleankey] = _create_input_group(self.types[-1],
                                                 parent[key])
            self["q" + qid] = self[cleankey]

            # Make sure reference by description is not too long
            max_ind = min( len(self.descs[-1]), self._MAX_DESC )
            self[ self.descs[-1][:max_ind] ] = self[cleankey]

        self.activeqid = "q" + get_activeqid(parent.file, parent)
        self.active = self[self.activeqid]

        # Organize qids, descriptions, dates and field names by active status
        # and date (active one first, then sorted by date from newest to oldest)
        index       = self.qids.index( self.activeqid )
        sortedqids  = [ self.qids.pop(index)  ]
        sorteddates = [ self.dates.pop(index) ]
        sortedtypes = [ self.types.pop(index) ]
        sorteddescs = [ self.descs.pop(index) ]

        if len(self.qids) > 0:
            # This sorts elements in y by sorted x
            sortedqids  += \
                    [x for _, x in sorted(zip(self.dates, self.qids))]
            sortedtypes += \
                    [x for _, x in sorted(zip(self.dates, self.types))]
            sorteddescs += \
                    [x for _, x in sorted(zip(self.dates, self.descs))]
            sorteddates += sorted(self.dates)

        self.qids  = sortedqids
        self.dates = sorteddates
        self.types = sortedtypes
        self.descs = sorteddescs

        self._freeze()

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
            if key >= len(self.qids):
                print("Index out of bounds. Maximum index is "
                      + len(self.qids))
                return None

            return super().__getitem__(self.qids[key])

    def __str__(self):
        """
        Get a table showing fields, qids, dates, and descriptions as a string.
        """
        string = ""
        for i in range(0, len(self.qids)):
            string += self.types[i]    + " "  + \
                      self.qids[i][1:] + " "  + \ # Without the 'q' prefix.
                      self.dates[i]    + "\n" + \
                      self.descs[i]
            if i < ( len(self.qids) - 1 ) :
                string += "\n"

        return string


class Ascot(_StandardNode):
    """
    Top node used for exploring the HDF5 file.
    """

    def __init__(self, fn):
        """
        Initialize the whole node structure recursively and create data objects.
        """
        super().__init__()
        self._hdf5fn = fn
        self.qids    = []
        self.descs   = []
        self.dates   = []

        with h5py.File(self._hdf5fn, "r") as h5:

            for key in h5.keys():
                if( key == "bfield" or key == "efield" or key == "options" or \
                    key == "marker" or key == "neutral" or key == "plasma" or \
                    key == "wall" ):
                    self[key] = _InputNode(h5[key])

            if "results" in h5:
                for key in h5["results"].keys():
                    qid = get_qid(key)
                    self.qids.append("q" + qid)
                    self.descs.append( get_desc(h5, h5["results"][key]) )
                    self.dates.append( get_date(h5, h5["results"][key]) )

                    cleankey =_remove_illegal_chars(key)
                    self[cleankey] = _create_run_group(h5["results"][key])
                    self["q" + qid]  = self[cleankey]

                    # Convert stored qids as references in rungroup
                    for i in ["options", "bfield", "efield", "marker", "plasma",
                              "neutral", "wall"]:
                        self[cleankey]._unfreeze()
                        self[cleankey][i] = self[i][self[cleankey][i]]
                        self[cleankey]._freeze()

                    # Make sure reference by description is not too long
                    max_ind = min( len(self.descs[-1]), self._MAX_DESC )
                    self[ self.descs[-1][:max_ind] ] = self[cleankey]

                self.activeqid = "q" + get_activeqid(h5, h5["results"])
                self.active = self[self.activeqid]

                # Organize qids, descriptions, dates and field names by active
                # status and date (active one first, then sorted by date from
                # newest to oldest)
                index       = self.qids.index( self.activeqid )
                sortedqids  = [ self.qids.pop(index)  ]
                sorteddates = [ self.dates.pop(index) ]
                sorteddescs = [ self.descs.pop(index) ]

                if len(self.qids) > 0:
                    # This sorts elements in y by sorted x
                    sortedqids  += \
                        [x for _, x in sorted(zip(self.dates, self.qids))]
                    sorteddescs += \
                        [x for _, x in sorted(zip(self.dates, self.descs))]
                    sorteddates += sorted(self.dates)

                self.qids  = sortedqids
                self.dates = sorteddates
                self.descs = sorteddescs

        self._freeze()

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
            if key >= len(self.qids):
                print("Index out of bounds. Maximum index is "
                      + len(self.qids))
                return None

            return super().__getitem__(self.qids[key])

    def __str__(self):
        """
        Overview of inputs and results in the HDF5 file in a string format.
        """
        string = "Inputs:\n"
        if(hasattr(self, "options")): string += "options\n"
        if(hasattr(self, "bfield" )): string += "bfield\n"
        if(hasattr(self, "efield" )): string += "efield\n"
        if(hasattr(self, "plasma" )): string += "plasma\n"
        if(hasattr(self, "neutral")): string += "neutral\n"
        if(hasattr(self, "wall"   )): string += "wall\n"

        string += "\nResults:\n"
        for i in range(0, len(self.qids)):
            string += "run" + " " + self.qids[i][1:] + " " + \
                      self.dates[i] + "\n" + self.descs[i]
            if i < ( len(self.qids) - 1 ) :
                string += "\n"

        return string

def _remove_illegal_chars(key):
    """
    Remove illegal characters from argument so it can be used as an attribute.

    Args:
        String to be cleaned
    Returns:
        String from which illegal characters are replaced or removed
    """
    key = key.replace(" ", "_")
    key = key.replace("-", "_")
    key = key.replace(".","")
    return key
