"""Package for processing ASCOT5 data and generating inputs.
"""
import os
import numpy as np
import unyt

from .ascot5io import Ascot5IO
from .ascotpy  import Ascotpy

from .ascotpy.libascot import _LIBASCOT
from .exceptions import *
from .routines.plotting import openfigureifnoaxes, line2d

# Define the unit system ascot uses and add our own unit types
unyt.define_unit("markers", 1*unyt.Dimensionless)
unyt.define_unit("particles", 1*unyt.Dimensionless)
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
                   asigma=False, switch=False):
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
        asigma  : str or bool or dict, optional
            The atomic data input or True to use the one that is active.
        switch  : bool, optional
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
                     isinstance(mhd, bool)     and isinstance(asigma, bool)):
                raise ValueError(
                    "Cannot specify input explicitly when run=True")
            initall = not any([bfield, efield, plasma, wall, neutral, boozer,
                               mhd, asigma])

            # Init either the active run or that which was requested
            run = self.data.active if isinstance(run, bool) else\
                  self.data["q"+run]

            # Find QIDs of the inputs to be initialized
            for inp in ["bfield", "efield", "plasma", "wall", "neutral",
                        "boozer", "mhd", "asigma"]:
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
                               mhd, asigma])

            for inp in ["bfield", "efield", "plasma", "wall", "neutral",
                        "boozer", "mhd", "asigma"]:
                if isinstance(args[inp], dict):
                    # Argument is a dictionary, presumably in the correct format
                    # It is simply passed forward to _init()
                    pass
                elif args[inp] and not inp in self.data:
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
                elif args[inp]:
                    # Argument is True, and data is present
                    args[inp] = getattr(self.data, inp).active.get_qid()
                else:
                    args[inp] = None

        self._init(
            self.data, bfield=args["bfield"], efield=args["efield"],
            plasma=args["plasma"], wall=args["wall"], neutral=args["neutral"],
            boozer=args["boozer"], mhd=args["mhd"], asigma=args["asigma"],
            switch=switch)

    def input_free(self, bfield=False, efield=False, plasma=False, wall=False,
                   neutral=False, boozer=False, mhd=False, asigma=False):
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
        asigma  : bool, optional
            Toggle to free the atomic data input.
        """
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        freeall = not any([bfield, efield, plasma, wall, neutral, boozer,
                           mhd, asigma])
        if freeall:
            self._free(bfield=True, efield=True, plasma=True, wall=True,
                       neutral=True, boozer=True, mhd=True, asigma=True)
        else:
            self._free(bfield=bfield, efield=efield, plasma=plasma, wall=wall,
                       neutral=neutral, boozer=boozer, mhd=mhd, asigma=asigma)

    def preflight_inputspresent(self):
        """Check required inputs are present for this run.
        """

        # Determine what we are simulating.
        runtype = "ascot"

        msg = []
        if runtype == "ascot":
            for inp in ["bfield", "efield", "plasma", "neutral", "options",
                        "marker", "wall", "boozer", "mhd", "asigma"]:
                if not hasattr(self.data, inp):
                    msg += ["Error: %s data is missing" % inp]

        return msg

    def preflight_optionsconsistent(self):
        """Check that options are consistent with inputs and each other.
        """
        msg = []
        opt = self.data.options.active

        #msg0 = opt.validatevalues()
        msg0 = []
        if len(msg0) > 0:
            msg += ["Error: following options had invalid parameters:"]
            msg += msg0

        opt = opt.read()
        if opt["ENABLE_MHD"] == 1 and \
           self.data.bfield.active.get_type() == "B_STS":
            msg += ["Error: cannot enable MHD for stellarators"]

        if opt["SIM_MODE"] in [1,2,3] and \
           self.data.marker.active.get_type() == "mrk_fl":
            msg += ["Error: SIM_MODE invalid for field-line markers"]

        # Make rough estimates on how much memory is consumed by diagnostics
        high_memory_consumption = 1e9 # In bits
        # Orbit = Npoint * Nfields (~10) * Nmrk * 8 bit
        orb_mem = opt["ORBITWRITE_NPOINT"] * 10 \
            * self.data.marker.active.read()["n"] * 8
        # Distributions = NR * Nz * ... * 8 bit
        rzp = opt["DIST_NBIN_R"] * opt["DIST_NBIN_Z"] * opt["DIST_NBIN_PHI"]
        rtp = opt["DIST_NBIN_RHO"] * opt["DIST_NBIN_THETA"] * opt["DIST_NBIN_PHI"]
        p2d = opt["DIST_NBIN_PPA"] * opt["DIST_NBIN_PPE"]
        p3d = opt["DIST_NBIN_PR"] * opt["DIST_NBIN_PZ"] * opt["DIST_NBIN_PPHI"]

        if opt["ENABLE_ORBITWRITE"] == 1 and orb_mem > high_memory_consumption:
            msg += ["Warning: orbit diagnostic memory consumption high (~" +
                    str(int(orb_mem / 1e9)) + "Gb)"]

        if opt["ENABLE_DIST_5D"] == 1 and rzp * p2d * 8 > high_memory_consumption:
            msg += ["Warning: 5D distribution memory consumption high (~" +
                    str(int(rzp * p2d * 8 / 1e9)) + "Gb)"]

        if opt["ENABLE_DIST_6D"] == 1 and rzp * p3d * 8 > high_memory_consumption:
            msg += ["Warning: 6D distribution memory consumption high (~" +
                    str(int(rzp * p3d * 8 / 1e9)) + "Gb)"]

        if opt["ENABLE_DIST_RHO5D"] == 1 and rtp * p2d * 8 > high_memory_consumption:
            msg += ["Warning: rho5D distribution memory consumption high (~" +
                    str(int(rtp * p2d * 8 / 1e9)) + "Gb)"]

        if opt["ENABLE_DIST_RHO6D"] == 1 and rtp * p3d * 8 > high_memory_consumption:
            msg += ["Warning: rho6D distribution memory consumption high (~" +
                    str(int(rtp * p3d * 8 / 1e9)) + "Gb)"]

        return msg

    def preflight_bfieldpsi0(self):
        """Checks whether psi0 given in input is actually extreme value.

        Because psi is interpolated with splines, there might be numerical error
        that causes psi (near the axis) to have more extreme value than psi0
        which is given in input. This leads to imaginary rho and termination
        of the simulation if marker ends up there.

        This check uses Monte Carlo method to 1. Draw phi 2. Evaluate axis (R,z)
        3. Draw random (R,z) coordinates within 10 cm of the axis. 4. Evaluate psi
        at that point and compare to psi0. Process repeats N times. Check passes
        if all evaluations are valid.
        """

        data = self.data.bfield.active.read()
        if "psi0" not in data:
            return []

        psi0 = data["psi0"]
        psi1 = data["psi1"]
        psi0 = -7.1

        N     = 10000
        phi   = np.random.rand(N,) * 360
        theta = np.random.rand(N,) * 2 * np.pi

        axis = self.input_eval(1, phi, 0, 0, "axis")
        z0 = axis["axisz"]
        r0 = axis["axisr"]

        R = 0.1 # 10 cm
        r = R * np.cos(theta) + r0
        z = R * np.sin(theta) + z0
        psi = self.input_eval(r, phi, z, 0, "psi")

        if psi0 < psi1 and any(psi < psi0):
            return ["Error: psi0 = %.2e but we found near axis that psi = %.2e"\
                    % (psi0, np.amin(psi)) ]
        elif psi0 > psi1 and any(psi > psi0):
            return ["Error: psi0 = %.2e but we found near axis that psi = %.2e"\
                    % (psi0, np.amax(psi)) ]
        return []

    @openfigureifnoaxes(projection=None)
    def preflight_plottopview(self, axes=None):
        """Plot top view of the machine showing Ip, Bphi, and markers.

        Assumes bfield is initialized in ascotpy.
        """
        r0, z0 = self.input_eval(1, 0, 0, 0, "axisr", "axisz")
        #z0 = axis["axisz"]
        #r0 = axis["axisr"]

        rmin = r0-r0/2
        rmax = r0+r0/2
        dphi = 10 * np.pi/180 # So that b and j quivers dont overlap

        r   = np.linspace(rmin, rmax, 10)
        phi = np.linspace(0, 360, 18, endpoint=False) * np.pi/180
        phi, r = np.meshgrid(phi, r)
        r   = r.ravel()
        phi = phi.ravel()

        br   = np.squeeze(self.input_eval(r, phi, z0, 0, "br"))
        bphi = np.squeeze(self.input_eval(r, phi, z0, 0, "bphi"))
        jr   = np.squeeze(self.input_eval(r, phi + dphi, z0, 0, "jr"))
        jphi = np.squeeze(self.input_eval(r, phi + dphi, z0, 0, "jphi"))

        x  = np.cos(phi) * r
        y  = np.sin(phi) * r
        bx = np.cos(phi) * br - np.sin(phi) * bphi
        by = np.sin(phi) * br + np.cos(phi) * bphi
        bnorm = np.sqrt(bx**2 + by**2)

        xj = np.cos(phi+dphi) * r
        yj = np.sin(phi+dphi) * r
        jx = np.cos(phi+dphi) * jr - np.sin(phi+dphi) * jphi
        jy = np.sin(phi+dphi) * jr + np.cos(phi+dphi) * jphi
        jnorm = np.sqrt(jx**2 + jy**2)

        axes.quiver(x,y, bx/bnorm, by/bnorm,
                    color="blue", scale=40)
        axes.quiver(xj,yj, jx/jnorm, jy/jnorm,
                    color="red", scale=40)

        axes.set_aspect("equal", adjustable="box")

        marker = self.data.marker.active.read()
        x  = np.cos(marker["phi"] * np.pi/180) * marker["r"]
        y  = np.sin(marker["phi"] * np.pi/180) * marker["r"]

        axes.scatter(x,y, s=1, c="black", zorder=-2)

        axes.legend(("Magnetic field", "Plasma current", "Markers"))

    def input_plotwallcontour(self, phi=0*unyt.deg, axes=None):
        """Plot intersection of the wall and the poloidal plane at the given
        toroidal angle.

        Parameters
        ----------
        phi : float
            Toroidal angle of the poloidal plane.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        """
        ls = self.data.wall.active.getwallcontour(phi=phi)
        line2d(ls[:,:,0], ls[:,:,1], c="black", axesequal=True, axes=axes)
