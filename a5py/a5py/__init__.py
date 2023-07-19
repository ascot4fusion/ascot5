"""Package for processing ASCOT5 data and generating inputs.
"""
import os
import numpy as np
import unyt

from .ascot5io import Ascot5IO
from .ascotpy  import Ascotpy

from .ascotpy.libascot import _LIBASCOT
from .exceptions import *

# Define the unit system ascot uses and add our own unit types
unyt.define_unit("markers", 1*unyt.Dimensionless)
unyt.define_unit("e", -unyt.electron_charge)
u=unyt.UnitSystem("ascot", "m", "atomic_mass_unit", "s", angle_unit="deg")
u["energy"] = "eV"
u["charge"] = "e"
u["magnetic_field"] = "T"

class Ascot(Ascotpy):
    """Primary tool for processing ASCOT5 data.

    Ascot object represents the data in a single HDF5 file. The class provides
    methods to:

    - add, modify, and delete the data in HDF5 and set descriptions.
    - access to input and output data, and premade plots for visualization.
    - calculate derived quantities from the data.
    - Python interface to ASCOT5 C routines to evaluate data as ASCOT5 does.
    - capability to simulate individual markers using the Python interface.

    In a nutshell, this object is all you need to work with ASCOT5 data unless
    one is doing something exotic.

    The data is accessed via :attr:`~a5py.Ascot.data` attribute and its
    attributes that form a treeview. For example, active magnetic field is
    accessed as :code:`Ascot.data.bfield.active`. See
    :class:`.Ascot5IO` for details.

    The active run is accessed as `Ascot5.data.active` and its methods contain
    most of the plotting and data access routines.

    Attributes
    ----------
    data : :class:`.Ascot5IO`
        Container for the HDF5 data.
    """

    def __init__(self, inputfile=None, create=False):
        """Initialize Ascot instance.

        Parameters
        ----------
        inputfile : str, optional
            Name of the HDF5 file or `None` to create an empty instance.
        """
        super().__init__()

        self._inputfile = None
        self.data       = None

        if create:
            if os.path.isfile(inputfile):
                raise AscotIOException("Cannot create file: file exists.")
            Ascot5IO._create_file(inputfile)

        self.file_load(inputfile)

    def file_getpath(self):
        """Return name of the HDF5 file from which this instance reads the data.

        Returns
        -------
        path : str
            Absolute path or `None` if this is an empty instance.
        """
        return self._inputfile

    def file_load(self, inputfile):
        """Open a new HDF5 file and free resources from the old one.

        Parameters
        ----------
        inputfile : str
            Path to the new HDF5 file.
        """
        # Free any resources that are currently being used (by the old file)
        if _LIBASCOT: self.input_free()

        if inputfile == None:
            # No input file so this is now an empty instance
            self._inputfile = None
            self.data = Ascot5IO(self)
            if _LIBASCOT:
                self._sim.hdf5_in = "".encode('UTF-8')
        else:
            inputfile = os.path.abspath(inputfile)
            try:
                self._inputfile = inputfile
                self.data = Ascot5IO(self)
                if _LIBASCOT:
                    self._sim.hdf5_in = inputfile.encode('UTF-8')
            except:
                self.inputfile = None
                self.data = Ascot5IO(self)
                if _LIBASCOT:
                    self._sim.hdf5_in = "".encode('UTF-8')

                # Pass exception from failed RootNode initialization
                raise

    def input_init(self, run=False, bfield=False, efield=False, plasma=False,
                   wall=False, neutral=False, boozer=False, mhd=False,
                   switch=False):
        """Initialize input so it can be accessed via Python interface.

        This method can be used in three ways:

        1. Without any arguments in which case all inputs are initialized. The
           initialized input is the one that has been set as active.
        2. Provide one of the input arguments (bfield, efield, plasma, wall,
           neutral, boozer, mhd) in which case only that input is initialized.
           The argument can either be QID of the input to be initialized, or
           `True` in which case the active input of that type is initialized.
        3. Set run = `True` or QID of a run group, in which case the inputs of
           that run group are initialized. If you only want a specific input
           to be initialized, set that argument as True, e.g. bfield = True
           would only initialize the magnetic field used by the given run.

        If the input is already initialized, nothing is done. In case there is
        already different input of same type initialized, an error is raised.

        If the input argument is a dict, it is used instead of reading the data from hdf5.

        Parameters
        ----------
        run     : str or bool, optional
            The run group or True to use the one that is active.
        bfield  : str or bool or dict, optional
            The magnetic field input or True to use the one that is active.
        efield  : str or bool or dict, optional
            The electric field input or True to use the one that is active.
        plasma  : str or bool or dict, optional
            The plasma input or True to use the one that is active.
        wall    : str or bool or dict, optional
            The wall input or True to use the one that is active.
        neutral : str or bool or dict, optional
            The neutral input or True to use the one that is active.
        boozer  : str or bool or dict, optional
            The boozer input or True to use the one that is active.
        mhd     : str or bool or dict, optional
            The MHD input or True to use the one that is active.
        switch : bool, optional
            If ``True``, no error is raised if there is already a different
            input initialized and instead the old input is free'd and the new
            one initialized.
        """
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")

        args = locals()
        if run:
            # Init inputs using the run; first check whether to init everything
            if not ( isinstance(bfield, bool)  and isinstance(efield, bool) and
                     isinstance(plasma, bool)  and isinstance(wall, bool)   and
                     isinstance(neutral, bool) and isinstance(boozer, bool) and
                     isinstance(mhd, bool) ):
                raise ValueError(
                    "Cannot specify input explicitly when run=True")
            initall = not any([bfield, efield, plasma, wall, neutral, boozer,
                               mhd])

            # Init either the active run or that which was requested
            run = self.data.active if isinstance(run, bool) else\
                  self.data["q"+run]

            # Find QIDs of the inputs to be initialized
            for inp in ["bfield", "efield", "plasma", "wall", "neutral",
                        "boozer", "mhd"]:
                if args[inp] and not inp in self.data:
                    raise AscotIOException("Input \"" + inp + "\" not present.")
                elif args[inp]:
                    args[inp] = getattr(run, inp).get_qid()
                elif initall and inp in self.data:
                    args[inp] = getattr(run, inp).get_qid()
                else:
                    args[inp] = None

        else:
            # Init given inputs; first check whether to init everything
            initall = not any([bfield, efield, plasma, wall, neutral, boozer,
                               mhd])

            for inp in ["bfield", "efield", "plasma", "wall", "neutral",
                        "boozer", "mhd"]:
                if args[inp] and not inp in self.data:
                    # Requested data not present
                    raise AscotIOException("Input \"" + inp + "\" not present.")
                elif initall and inp in self.data:
                    # Init all inputs that are present
                    args[inp] = getattr(self.data, inp).active.get_qid()
                elif isinstance(args[inp], str):
                    # Argument is QID already but check that it is in the data
                    parent = getattr(self.data, inp)
                    if not ("q"+args[inp]) in parent:
                        raise AscotIOException(
                            "Input \"" + inp + "/" + args[inp]
                            + "\" not present.")
                elif isinstance(args[inp], dict):
                    # Argument is a dictionary, presumably in the correct format
                    # It is simply passed forward to _init()
                    pass
                elif args[inp]:
                    # Argument is True, and data is present
                    args[inp] = getattr(self.data, inp).active.get_qid()
                else:
                    args[inp] = None

        self._init(
            self.data, bfield=args["bfield"], efield=args["efield"],
            plasma=args["plasma"], wall=args["wall"], neutral=args["neutral"],
            boozer=args["boozer"], mhd=args["mhd"], switch=switch)

    def input_free(self, bfield=False, efield=False, plasma=False, wall=False,
                   neutral=False, boozer=False, mhd=False):
        """Free input used by the Python interface.

        Arguments toggle which input fields are free'd. If called without
        providing any arguments, then all input fields are free'd.

        Parameters
        ----------
        bfield  : bool, optional
            Toggle to free the magnetic field.
        efield  : bool, optional
            Toggle to free the electric field input.
        plasma  : bool, optional
            Toggle to free the plasma input.
        wall    : bool, optional
            Toggle to free the wall input.
        neutral : bool, optional
            Toggle to free the neutral input.
        boozer  : bool, optional
            Toggle to free the boozer input.
        mhd     : bool, optional
            Toggle to free the MHD input.
        """
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        freeall = not any([bfield, efield, plasma, wall, neutral, boozer,
                           mhd])
        if freeall:
            self._free(bfield=True, efield=True, plasma=True, wall=True,
                       neutral=True, boozer=True, mhd=True)
        else:
            self._free(bfield=bfield, efield=efield, plasma=plasma, wall=wall,
                       neutral=neutral, boozer=boozer, mhd=mhd)
