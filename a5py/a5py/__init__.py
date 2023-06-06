"""
Contains Ascot class definition and other package wide stuff.
"""
import os
import numpy as np

from .ascot5io import Ascot5IO
from .ascotpy  import Ascotpy

from .ascotpy.libascot import LIBRARY_AVAILABLE
from .exceptions import *

# Define the unit system ascot uses and add our own unit types
import unyt

unyt.define_unit("markers", 1*unyt.Dimensionless)
unyt.define_unit("e", -unyt.electron_charge)

u=unyt.UnitSystem("ascot", "m", "atomic_mass_unit", "s", angle_unit="deg")
u["energy"] = "eV"
u["charge"] = "e"
u["magnetic_field"] = "T"



class Ascot(Ascotpy):
    """
    Primary tool for processing ASCOT5 data.

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
    :class:`~a5py.ascot5io.RootNode` for details.

    The active run is accessed as `Ascot5.data.active` and its methods contain
    most of the plotting and data access routines.

    Attributes
    ----------
        data : :class:`~a5py.ascot5io.RootNode`
            Container for the HDF5 data.
    """

    def __init__(self, inputfile=None, create=False):
        """
        Initialize Ascot instance.

        Args:
            inputfile : str, optional
               Name of the HDF5 file or None to create an empty instance.
        """
        super().__init__()

        self._inputfile = None
        self.data       = None

        if create:
            Ascot5IO._create_file(inputfile)

        self.file_load(inputfile)

    def file_getpath(self):
        """
        Return name of the HDF5 file from which this instance reads the data.

        Returns:
            Absolute path or None if this is an empty instance.
        """
        return self._inputfile

    def file_load(self, inputfile):
        """
        Open a new HDF5 file and free resources from the old one.

        Args:
            inputfile : str
                Name of the new HDF5 file.
        """
        # Free any resources that are currently being used (by the old file)
        self.input_free()

        if inputfile == None:
            # No input file so this is now an empty instance
            self._inputfile = None
            self.data = Ascot5IO(self)
            if LIBRARY_AVAILABLE:
                self._sim.hdf5_in = "".encode('UTF-8')
        else:
            inputfile = os.path.abspath(inputfile)
            try:
                self._inputfile = inputfile
                self.data = Ascot5IO(self)
                if LIBRARY_AVAILABLE:
                    self._sim.hdf5_in = inputfile.encode('UTF-8')
            except:
                self.inputfile = None
                self.data = Ascot5IO(self)
                if LIBRARY_AVAILABLE:
                    self._sim.hdf5_in = "".encode('UTF-8')

                # Pass exception from failed RootNode initialization
                raise

    def input_init(self, run=False, bfield=False, efield=False, plasma=False,
                   wall=False, neutral=False, boozer=False, mhd=False):
        """
        Initialize input so it can be accessed via Python interface.

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

        Args:
            run     : bool or str ,optional
                The run group or True to use the one that is active.
            bfield  : bool, optional
                The magnetic field input or True to use the one that is active.
            efield  : bool, optional
                The electric field input or True to use the one that is active.
            plasma  : bool, optional
                The plasma input or True to use the one that is active.
            wall    : bool, optional
                The wall input or True to use the one that is active.
            neutral : bool, optional
                The neutral input or True to use the one that is active.
            boozer  : bool, optional
                The boozer input or True to use the one that is active.
            mhd     : bool, optional
                The MHD input or True to use the one that is active.
        """
        if run:
            # Init inputs using the run; first check whether to init everything
            if not ( isinstance(bfield, bool)  and isinstance(efield, bool) and
                     isinstance(plasma, bool)  and isinstance(wall, bool)   and
                     isinstance(neutral, bool) and isinstance(boozer, bool) and
                     isinstance(mhd, bool) ):
                raise
            initall = not any([bfield, efield, plasma, wall, neutral, boozer,
                               mhd])

            # Init either the active run or that which was requested
            run = self.data.active if isinstance(run, bool) else\
                  self.data["q"+run]

            # Find QIDs of the inputs to be initialized
            bfield  = run.bfield.get_qid()  if (initall or bfield)  else False
            efield  = run.efield.get_qid()  if (initall or efield)  else False
            plasma  = run.plasma.get_qid()  if (initall or plasma)  else False
            wall    = run.wall.get_qid()    if (initall or wall)    else False
            neutral = run.neutral.get_qid() if (initall or neutral) else False
            boozer  = run.boozer.get_qid()  if (initall or boozer)  else False
            mhd     = run.mhd.get_qid()     if (initall or mhd)     else False

        else:
            # Init given inputs; first check whether to init everything
            initall = not any([bfield, efield, plasma, wall, neutral, boozer,
                               mhd])

            # Helper function to get the group QID
            def inputqid(key, val):
                if isinstance(val, str):
                    # val is QID already
                    return val
                else:
                    # val is boolean, return active QID if True
                    if val or initall:
                        print((key,val,initall))
                        return getattr(self.data, key).active.get_qid()
                    else:
                        return False

            bfield  = inputqid("bfield",  bfield)
            efield  = inputqid("efield",  efield)
            plasma  = inputqid("plasma",  plasma)
            wall    = inputqid("wall",    wall)
            neutral = inputqid("neutral", neutral)
            boozer  = inputqid("boozer",  boozer)
            mhd     = inputqid("mhd",     mhd)

        self._init(self.data, bfield=bfield, efield=efield, plasma=plasma,
                   wall=wall, neutral=neutral, boozer=boozer, mhd=mhd)

    def input_free(self, bfield=False, efield=False, plasma=False, wall=False,
                   neutral=False, boozer=False, mhd=False):
        """
        Free input used by the Python interface.

        Arguments toggle which input fields are free'd. If called without
        providing any arguments, then all input fields are free'd.

        Args:
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
        self._free(bfield=bfield, efield=efield, plasma=plasma, wall=wall,
                   neutral=neutral, boozer=boozer, mhd=mhd)

    def input_plotrz(self, r, phi, z, t, quantity, clim=(None, None), axes=None,
                     cax=None):
        """
        Plot input quantity on a Rz plane at given coordinates.

        Args:
            r : np.array
                R abscissa where data is evaluated and plotted [m].
            z : np.array
                z abscissa where data is evaluated and plotted [m].
            phi : float
                phi coordinate where data is evaluated [rad].
            t : float
                time coordinate where data is evaluated [s].
            axes : Axes
                plot on these axes instead of creating a new figure.
            clim :
                tuple with minimum and maximum color values.
        """

        out = self.input_eval(r, phi, z, t, quantity, grid=True)
        out = np.transpose(out[:,0,:,0])

        out = np.ma.masked_invalid(out)
        if clim[0] is None:
            clim[0] = np.nanmin(out)
        if clim[1] is None:
            clim[1] = np.nanmax(out)

        mesh = axes.pcolormesh(r, z, out, vmin=clim[0], vmax=clim[1])
        axes.patch.set(hatch='x', edgecolor=[0.9, 0.9, 0.9])

        plt.colorbar(mesh, cax=cax)
        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(r[0], r[-1])
        axes.set_ylim(z[0], z[-1])

    def input_plotseparatrix(self, r, phi, z, t, axes=None):
        """
        Plot plasma separatrix
        """
        out = self.input_eval(r, phi, z, t, "rho", grid=True)

        mesh = axes.contour(r, z, np.transpose(out[:,0,:,0]), [1],
                            colors='black',zorder=1)
