"""Python interface to libascot.so.
"""
import ctypes
import unyt
import numpy as np
import scipy
import warnings
import wurlitzer # For muting libascot.so

import a5py.physlib as physlib
import a5py.routines.plotting as a5plt
from a5py.routines.plotting import openfigureifnoaxes, plt
from a5py.exceptions import AscotInitException

from .libascot    import LibAscot, _LIBASCOT, PTR_RFOF, _get_struct_class
from .libsimulate import LibSimulate
from .libproviders import LibProviders

if _LIBASCOT:
    from . import ascot2py

class Ascotpy(LibAscot, LibSimulate, LibProviders):
    """Class with methods to initialize and access the data via Python.

    This class will be inherited by Ascot class. Here we hide all the dirty
    implementation of the libascot.so interface in private methods and
    attributes so that the Ascot front is clean for the users.

    Attributes
    ----------
    _sim
        Simulation data struct.
    _nmrk
        Number of markers currently in the marker array.
    _inistate
        Marker input array for interactive simulations.
    _endstate
        Marker output array for interactive simulations.
    _diag_occupied
        Flag indicating if diagnostics output contains data.
    _mute
        Mute output from libascot.so: "yes" - stdout and stderr both muted,
        "no" - output is not muted, "err" - stderr is printed.
    """

    DUMMY_QID = "".encode('UTF-8')
    """Flag to use in _sim.qid_* to indicate that the input is not initialized.
    """

    def __init__(self):
        """Initialize pointers to data that is being passed between Python and C
        via libascot.
        """
        if not _LIBASCOT:
            return

        # Initialize attributes
        self._nmrk = ctypes.c_int32()
        self._inistate = ctypes.POINTER(_get_struct_class("particle_state"))()
        self._endstate = ctypes.POINTER(_get_struct_class("particle_state"))()
        self._diag_occupied = False

        self._sim = _get_struct_class("sim_data")()
        self._mute = "no"

    def _initmpi(self, mpirank, mpisize, mpiroot=0):
        """Initialize MPI data.

        Parameters
        ----------
        mpirank : int
            MPI rank.
        mpisize : int
            MPI size.
        mpiroot : int, optional
            Rank of the root process.

            Only changed when running pseudo-MPI simulation.
        """
        self._sim.mpi_rank = ctypes.c_int32(mpirank)
        self._sim.mpi_size = ctypes.c_int32(mpisize)
        self._sim.mpi_root = ctypes.c_int32(mpiroot)

    def _init(self, data, bfield=None, efield=None, plasma=None,
              wall=None, neutral=None, boozer=None, mhd=None, asigma=None,
              nbi=None, switch=False):
        """Read and initialize input data so it can be accessed
        by libascot.

        The input is initialized by reading the data from the HDF5 file
        if a boolean value or QID is provided as an argument. The whole reading
        and initialization happens on the C-side.

        If the argument is a dictionary containing the input data (i.e. what
        would be provided to a corresponding write_hdf5 function) then the HDF5
        file is completely ignored and the input is initialized from
        the dictionary completely on the Python side, using ctypes, and only
        the init is done in C. See :mod:`a5py.ascotpy.libproviders` for
        details.

        Parameters
        ----------
        data : Ascot5IO
            Ascot5 data on disk which is used in initialization.
        bfield : str or dict
            QID of the magnetic field to be initialized or the data as a
            dictionary.
        efield : str or dict
            QID of the electric field to be initialized or the data as a
            dictionary.
        plasma : str or dict
            QID of the plasma data to be initialized or the data as a
            dictionary.
        wall : str or dict
            QID of the wall data to be initialized or the data as a dictionary.
        neutral : str or dict
            QID of the neutral data to be initialized or the data as a
            dictionary.
        boozer : str or dict
            QID of the boozer to be initialized or the data as a dictionary.
        mhd : str or dict
            QID of the MHD data to be initialized or the data as a dictionary.
        asigma : str or dict
            QID of the atomicdata to be initialized or the data as a dictionary.
        nbi : str or dict
            QID of the NBI data to be initialized or the data as a dictionary.
        switch : bool
            If ``True``, free input that has been
        """
        # Iterate through all inputs and mark those that are initialized
        inputs2read = ctypes.c_int32()
        args = locals() # Contains function arguments and values in a dictionary
        to_be_provided = [] # Inputs to be directly injected (provided)
        for inp in ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd", "asigma", "nbi"]:
            if args[inp] is None:
                # This input is not going to be initialized
                continue
            currentqid = getattr(self._sim, "qid_" + inp)

            # Inject/provide this input
            if isinstance(args[inp],dict):
                if currentqid != Ascotpy.DUMMY_QID:
                    # Other input of same type is initialized. Raise error or
                    # free the old input.
                    if not switch:
                        raise AscotInitException(
                            "Cannot initialize " + inp
                            + ": other input is already initialized.\n"
                            "Either use \"switch=True\" or free("+inp+"=True)")

                    self._free(**{inp : True})
                to_be_provided.append(inp)
                continue

            # Convert QID strings to bytes
            args[inp] = args[inp].encode("UTF-8")
            if args[inp] == currentqid:
                # This input is already initialized
                continue

            if currentqid != Ascotpy.DUMMY_QID:
                # Other input of same type is initialized. Raise error or free
                # the old input.
                if not switch:
                    raise AscotInitException(
                        "Cannot initialize " + inp
                        + ": other input is already initialized.\n"
                        "Either use \"switch=True\" or free(" + inp + "=True)")

                self._free(**{inp : True})

            # Mark this input to be read
            inputs2read.value |= getattr(ascot2py, "hdf5_input_" + inp)

        # Separate loop to mark QIDs of the inputs that are read. (We couldn't
        # do this in the previous loop as that might get interrupted by
        # an exception, in which case sim.qid_* would point to data which is not
        # initialized.
        for inp in ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd", "asigma", "nbi"]:
            if inputs2read.value & getattr(ascot2py, "hdf5_input_" + inp):
                setattr(self._sim, "qid_" + inp, args[inp])

        def hdf5init():
            """Simple wrapper.
            """
            if(inputs2read.value == 0): return
            ascot2py.hdf5_interface_read_input(
                ctypes.byref(self._sim),
                inputs2read,
                None, # Marker array (ignore)
                None  # Number of markers that were read (ignore)
            )

        if self._mute == "no":
            hdf5init()
        else:
            with wurlitzer.pipes() as (out, err):
                hdf5init()
            err = err.read()
            if self._mute == "err" and len(err) > 1: print(err)

        for inp in to_be_provided:
            # Initialize data from dictionaries here
            getattr(self, "_provide_" + inp)(**args[inp])

    def _free(self, bfield=False, efield=False, plasma=False, wall=False,
              neutral=False, boozer=False, mhd=False, asigma=False):
        """Free input data initialized in C-side.
        """
        args = locals() # Contains function arguments and values in a dictionary

        # Iterate through all inputs and free the data if the corresponding
        # argument is True
        for inp in ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd", "asigma"]:
            if args[inp] and \
                getattr(self._sim, "qid_" + inp) != Ascotpy.DUMMY_QID:
                # Deallocate the data allocated in C side
                if inp == "bfield":
                    _LIBASCOT.B_field_free(ctypes.byref(self._sim.B_data))
                elif inp == "efield":
                    _LIBASCOT.E_field_free(ctypes.byref(self._sim.E_data))
                else:
                    data = getattr(self._sim, inp + "_data")
                    getattr(_LIBASCOT, inp + "_free")(ctypes.byref(data))

                # Set QID to dummy value
                setattr(self._sim, "qid_" + inp, Ascotpy.DUMMY_QID)

    def _requireinit(self, *inputs):
        """Raise error if given input parent is not initialized.

        Parameters
        ----------
        *inputs : [str]
            Name(s) of the input parent group(s) e.g. "bfield".

        Raises
        ------
        AscotInitException
            If input parent group has no initialized data.
        """
        for inp in inputs:
            qid = getattr(self._sim, "qid_" + inp)
            if qid == Ascotpy.DUMMY_QID:
                raise AscotInitException(inp + " is not initialized")

    def input_init_rfof(self):
        """
        """
        if self._rfof_initialized:
            raise AscotInitException("RFOF input is already initialized")
        fun = _LIBASCOT.rfof_init_offload
        fun.restype  = None
        fun.argtypes = [PTR_RFOF]
        fun(ctypes.byref(self._sim.rfof_data))
        self._rfof_initialized = True

    def input_free_rfof(self):
        """
        """
        if not self._rfof_initialized:
            raise AscotInitException("RFOF input is not initialized")

        self._rfof_initialized = False

    def input_initialized(self):
        """Get inputs that are currently initialized.

        Returns
        -------
        outputs : dict [str, str]
            Name of the input parent groups that are initialized and QID of the
            initialized input.
        """
        out = {}
        for inp in ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd", "asigma", "nbi"]:
            qid = getattr(self._sim, "qid_" + inp)
            if qid != Ascotpy.DUMMY_QID:
                out[inp] = qid.decode("utf-8")
        return out

    def input_eval_list(self, show=True):
        """Return quantities that can be evaluated from inputs.

        The returned quantities are those that can be evaluated from the inputs
        with `input_eval` method. Units are not listed as that function returns
        quantities with units included.

        The number ion species in plasma input is not fixed so "ni" refers to
        i'th plasma species as listed by `input_getplasmaspecies` (index starts
        from 1).

        Parameters
        ----------
        show : bool, optional
            Print output.

        Returns
        -------
        quantities : dict [str, str]
            Name of each quantity and a short description.
        """

        out = {
            "rho":
            "Square root of normalized poloidal flux",
            "psi":
            "Poloidal flux",
            "br":
            "Magnetic field R component (not including MHD)",
            "bphi":
            "Magnetic field phi component (not including MHD)",
            "bz":
            "Magnetic field z component (not including MHD)",
            "bpol":
            "Magnetic field poloidal component (not including MHD)",
            "bnorm":
            "Magnetic field norm (not including MHD)",
            "brdr":
            "d/dr of magnetic field R component (not including MHD)",
            "brdphi":
            "d/dphi of magnetic field R component (not including MHD)",
            "brdz":
            "d/dz of magnetic field R component (not including MHD)",
            "bphidr":
            "d/dr of magnetic field phi component (not including MHD)",
            "bphidphi":
            "d/dphi of magnetic field phi component (not including MHD)",
            "bphidz":
            "d/dz of magnetic field phi component (not including MHD)",
            "bzdr":
            "d/dr of magnetic field z component (not including MHD)",
            "bzdphi":
            "d/dphi of magnetic field z component (not including MHD)",
            "bzdz":
            "d/dz of magnetic field z component (not including MHD)",
            "divb":
            "Magnetic field divergence (not including MHD)",
            "gradbr":
            "Magnetic field gradient R component (not including MHD)",
            "gradbphi":
            "Magnetic field gradient phi component (not including MHD)",
            "gradbz":
            "Magnetic field gradient z component (not including MHD)",
            "curlbr":
            "Magnetic field curl R component (not including MHD)",
            "curlbphi":
            "Magnetic field curl phi component (not including MHD)",
            "curlbz":
            "Magnetic field curl z component (not including MHD)",
            "jr":
            "Current density R component",
            "jphi":
            "Current density phi component",
            "jz":
            "Current density z component",
            "jnorm":
            "Current density",
            "er":
            "Electric field R component (not including MHD)",
            "ephi":
            "Electric field phi component (not including MHD)",
            "ez":
            "Electric field z component (not including MHD)",
            "ne":
            "Electron density",
            "te":
            "Electron temperature",
            "zeff":
            "Efective charge, i.e., sum (n_i * Z_i^2) / ne",
            "n0":
            "Neutral density",
            "alphaeig":
            "Eigenfunction of the magnetic field MHD perturbation",
            "phieig":
            "Eigenfunction of the electric field MHD perturbation",
            "mhd_br":
            "R component of B-field perturbation due to MHD",
            "mhd_bphi":
            "phi component of B-field perturbation due to MHD",
            "mhd_bz":
            "z component of B-field perturbation due to MHD",
            "mhd_er":
            "R component of E-field perturbation due to MHD",
            "mhd_ephi":
            "phi component of E-field perturbation due to MHD",
            "mhd_ez":
            "z component of E-field perturbation due to MHD",
            "mhd_phi":
            "Electric potential due to MHD",
            "db/b (mhd)":
            """|B_MHD| / |B| where B is evaluated from the magnetic
            field input and both contain all components""",
            "psi (bzr)":
            "Boozer poloidal flux",
            "rho (bzr)":
            "Boozer radial coordinate",
            "theta":
            "Boozer poloidal coordinate",
            "zeta":
            "Boozer toroidal coordinate",
            "dpsidr (bzr)":
            "d/dR of Boozer psi",
            "dpsidphi (bzr)":
            "d/dphi of Boozer psi",
            "dpsidz (bzr)":
            "d/dz of Boozer psi",
            "dthetadr":
            "d/dr of Boozer theta",
            "dthetadphi":
            "d/dphi of Boozer theta",
            "dthetadz":
            "d/dz of Boozer theta",
            "dzetadr":
            "d/dr of Boozer zeta",
            "dzetadphi":
            "d/dphi of Boozer zeta",
            "dzetadz":
            "d/dz of Boozer zeta",
            "qprof":
            "Local value of the q-profile (safety factor)",
            "bjac":
            "Determinant of the Jacobian in Boozer coordinate transformation",
            "bjacxb2":
            "bjac multiplied by B^2 which is constant along flux surfaces",
        }

        #nion, mass, charge, anum, znum = self.input_getplasmaspecies()
        #for i in range(nion):
        #    out["ni" + str(i+1)] = "Ion species (anum, znum) = (%d, %d) density" \
        #        % (anum[i], znum[i])
        #    out["ti" + str(i+1)] = "Ion species (anum, znum) = (%d, %d) temperature" \
        #        % (anum[i], znum[i])

        if show:
            for name, desc in out.items():
                print(name.ljust(15) + " : " + desc)

        return out

    @physlib.parseunits(r="m", phi="deg", z="m", t="s")
    def input_eval(self, r, phi, z, t, *qnt, grid=False):
        """Evaluate input quantities at given coordinates.

        This method uses ASCOT5 C-routines for evaluating inputs exactly as
        they are evaluated during a simulation. See :meth:`input_eval_list` for
        a complete list of available input and derived quantities.

        Parameters
        ----------
        r : array_like
            R coordinates where data is evaluated.
        phi : array_like
            phi coordinates where data is evaluated.
        z : array_like
            z coordinates where data is evaluated.
        t : array_like
            Time coordinates where data is evaluated.
        *qnt : str
            Names of evaluated quantities.
        grid : boolean, optional
            Treat input coordinates as abscissae and return the evaluated
            quantities on a grid instead.

        Returns
        -------
        val : array_like (n,) or (nr,nphi,nz,nt)
            Quantity evaluated in each queried point.

            NaN is returned at those points where the evaluation failed.

            4D matrix is returned when grid=``True``.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        r   = unyt.unyt_array(r).astype(dtype=np.double).ravel()
        phi = unyt.unyt_array(phi).astype(dtype=np.double).ravel()
        z   = unyt.unyt_array(z).astype(dtype=np.double).ravel()
        t   = unyt.unyt_array(t).astype(dtype=np.double).ravel()

        phi = phi.to("rad")

        if grid:
            arrsize = (r.size, phi.size, z.size, t.size)
            r, phi, z, t = np.meshgrid(r, phi, z, t, indexing="ij")
            r   = r.ravel()
            phi = phi.ravel()
            z   = z.ravel()
            t   = t.ravel()
        else:
            # Not a grid so check that dimensions are correct (and make
            # single-valued vectors correct size)
            arrsize = np.amax(np.array([r.size, phi.size, z.size, t.size]))
            errmsg = "Input arrays have inconsistent sizes ({}, {}, {}, {})"
            assert (r.size == 1 or r.size == arrsize) and  \
                (phi.size == 1 or phi.size == arrsize) and \
                (z.size == 1 or z.size == arrsize) and     \
                (t.size == 1 or t.size == arrsize),        \
                errmsg.format(r.size, phi.size, z.size, t.size)

            if r.size == 1:
                r = r[0]*np.ones((arrsize,))
            if phi.size == 1:
                phi = phi[0]*np.ones((arrsize,))
            if z.size == 1:
                z = z[0]*np.ones((arrsize,))
            if t.size == 1:
                t = t[0]*np.ones((arrsize,))

        out = {}
        if any(q in qnt for q in ["rho", "psi", "rhodpsi", "psidr", "psidphi",
                                  "psidz"]):
            out.update(**self._eval_bfield(r, phi, z, t, evalrho=True))
        if any(q in qnt for q in ["br", "bphi", "bz", "brdr", "brdphi", "brdz",
                                  "bphidr", "bphidphi", "bphidz", "bzdr",
                                  "bzdphi", "bzdz", "divb", "bnorm",
                                  "jnorm", "jr", "jphi", "jz", "gradbr",
                                  "gradbphi", "gradbz", "curlbr", "curlbphi",
                                  "curlbz"]):
            out.update(self._eval_bfield(r, phi, z, t, evalb=True))
            out["divb"] = out["br"]/r + out["brdr"] + out["bphidphi"]/r \
                + out["bzdz"]
            out["bnorm"] = np.sqrt(out["br"]**2 + out["bphi"]**2 + out["bz"]**2)
            out["jr"]    = (out["bzdphi"]/r - out["bphidz"]) / unyt.mu_0
            out["jphi"]  = (out["brdz"] - out["bzdr"]) / unyt.mu_0
            out["jz"]    = (out["bphi"]/r + out["bphidr"] - out["brdphi"]/r) \
                            / unyt.mu_0
            out["jnorm"] = np.sqrt(out["jr"]**2 + out["jphi"]**2 + out["jz"]**2)
            out["gradbr"]   = (out["br"]*out["brdr"] + out["bphi"]*out["bphidr"]
                               + out["bz"]*out["bzdr"]) / out["bnorm"]
            out["gradbphi"] = (out["br"]*out["brdphi"]
                               + out["bphi"]*out["bphidphi"]
                               + out["bz"]*out["bzdphi"]) / (out["bnorm"]*r)
            out["gradbz"]   = (out["br"]*out["brdz"] + out["bphi"]*out["bphidz"]
                               + out["bz"]*out["bzdz"]) / out["bnorm"]
            out["curlbr"]   = out["bzdphi"] / r - out["bphidz"]
            out["curlbphi"] = out["brdz"] - out["bzdr"]
            out["curlbz"]   = (out["bphi"] - out["brdphi"]) / r + out["bphidr"]
        if any(q in qnt for q in ["axisr", "axisz", "rminor", "bpol"]):
            out.update(self._eval_bfield(r, phi, z, t, evalaxis=True, evalb=True))
            out["rminor"] = np.sqrt(  ( out["axisr"] - r )**2
                                    + ( out["axisz"] - z )**2 )
            dr = r-out["axisr"]
            dz = z-out["axisz"]
            bpol_sign = np.sign(dr*out["bz"] - dz*out["br"])
            out["bpol"] = bpol_sign*np.sqrt(out["bphi"]**2 + out["bz"]**2)
        if any(q in qnt for q in ["er", "ephi", "ez"]):
            out.update(self._eval_efield(r, phi, z, t))
        if any(q in qnt for q in ["n0"]):
            out.update(self._eval_neutral(r, phi, z, t))
        if any(q in qnt for q in ["psi (bzr)", "theta", "zeta", "dpsidr (bzr)",
                                  "dpsidphi (bzr)", "dpsidz (bzr)", "dthetadr",
                                  "dthetadphi", "dthetadz", "dzetadr",
                                  "dzetadphi", "dzetadz", "rho (bzr)"]):
            out.update(self._eval_boozer(r, phi, z, t))
        if any(q in qnt for q in ["qprof", "bjac", "bjacxb2"]):
            out.update(self._eval_boozer(r, phi, z, t, evalfun=True))
        if any(q in qnt for q in ["alphaeig", "phieig"]):
            out.update(self._eval_mhd(r, phi, z, t, evalpot=True))
        if any(q in qnt for q in ["mhd_br", "mhd_bphi", "mhd_bz", "mhd_er",
                                  "mhd_ephi", "mhd_ez", "mhd_phi"]):
            out.update(self._eval_mhd(r, phi, z, t))
        if any(q in qnt for q in ["db/b (mhd)"]):
            bpert = self._eval_mhd(r, phi, z, t)
            b = self._eval_bfield(r, phi, z, t, evalb=True)
            bpert = np.sqrt(  bpert["mhd_br"]**2
                            + bpert["mhd_bphi"]**2
                            + bpert["mhd_bz"]**2)
            b = np.sqrt(b["br"]**2 + b["bphi"]**2 + b["bz"]**2)
            out["db/b (mhd)"] = bpert/b

        ni = ["ni" + str(i+1) for i in range(99)]
        ti = ["ti" + str(i+1) for i in range(99)]
        if any(q in qnt for q in ["ne", "te", "zeff"] + ni + ti):
            out.update(self._eval_plasma(r, phi, z, t))

        for q in list(out.keys()):
            if q not in qnt:
                del out[q]
            elif grid:
                out[q] = np.reshape(out[q], arrsize)
        if len(qnt) == 1:
            return out[qnt[0]]
        else:
            return [out[q] for q in qnt]

    def input_plotalongaxis(self, qnt, nphi=None, phimin=0, phimax=360,
                            t=0.0*unyt.s, axes=None):
        """Plot qnt on the magnetic axis as a function of phi.


        By default plots qnt in the phigrid of the axis inputs.
        However, if nphi is given, qnt is plotted against
        np.linspace(phimin, phimax, nphi)[:-1]

        Parameters
        ----------
        qnt : str
            Name of the quantity to be plotted along the magnetic axis.
        nphi : int, optional
            number of phi points for plotting
        phimin : float, optional
            Minimum phi value at which qnt is plotted.
        phimax : float, optional
            "Maximum" phivalua for the plot. See the comment above parameters.
        """

        if nphi is None:
            # Phivec corresponds to the points where the axis input is given.
            b = self.data.bfield.active.read()
            phi_vec = np.linspace(b["axis_phimin"][0],b["axis_phimax"][0],
                                  b["axis_nphi"][0]+1)[:-1]*unyt.deg
            qnt_vec = np.empty(b["axis_nphi"][0])
            axisr = b["axisr"]*unyt.m
            axisz = b["axisz"]*unyt.m
        else:
            # Phivev is calculated based on the given resolution, nphi
            phi_vec = np.linspace(phimin, phimax, nphi+1)[:-1]*unyt.deg
            qnt_vec = np.empty(nphi)
            dummyarray = np.ones(nphi)
            out = self._eval_bfield(
                dummyarray, phi_vec.to("rad").v, dummyarray, t, evalaxis=True)
            axisr = out["axisr"]    # these are already in unyt.m
            axisz = out["axisz"]    # these are already in unyt.m
        qnt_vec = self.input_eval(axisr, phi_vec, axisz, t, qnt)
        a5plt.line2d([phi_vec], [qnt_vec], xlabel="phi [deg]",
                     ylabel=f"{qnt:s} [{str(qnt_vec.units):s}]", axes=axes)

    def input_check_bfield_seam(self,
                                qnt,
                                r=None,
                                z=None,
                                axes=None,
                                centremostpercentage=10,
                                diff=True,
                                forceplot=False,
                                n_centremost_elements_for_fit=None,
                                plot_spline_interpolated=False,
                                ):
        """
        Tries to check if an (R,z) plane is missing near the "seam" of the
        bfield input.

        This function is intended to be used to verify that the correct number
        of (R,z) planes is provided when using a stellarator bfield input. To
        recap, psi and Bfield components should be given at the points
        np.linspace(phimin, phimax, nphi+1)[:-1]. (It is possible to
        independently specify phimin, phimax and nphi for both psi and the
        bfield components.)

        If an (R,z) plane might be missing, plots qnt as a function of phi,
        evaluated only at the grid points where the input is given.

        Parameters
        ----------
        qnt : str
            One of the following: {"psi", "br", "bphi", "bz"}
        r : float, optional
            R coordinate where data is evaluated. (The functions uses the R
            value of the input that is closest to this requested one) [m]
        z : float, optional
            z coordinate where data is evaluated. (The functions uses the z
            value of the input that is closest to this requested one) [m]
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        centremostpercentage : float, optional
            Percentage of the grid closest to the requested r and z.
            (Default is 10, i.e. 5 % of the grid on either side)
        diff : bool, optional
            If True, plots np.diff(qnt)/np.diff(phi) instead of qnt. A possible
            kink may be more clearly visible in the diff(qnt) plot.
        forceplot : bool, optional
            If True, plots the requested quantity even if no (R,z) plane is
            deemed missing.
        n_centremost_elements_for_fit : int, optional
            Number of elements from the centre of the quantity array to use for
            least square sum fitting a polynomial. Default is None, in which
            case the whole array, consisting of 'centremostpercentage', is used.
        plot_spline_interpolated : bool, optional
            If True, plots the spline interpolated version of qnt.
        """

        if qnt not in ["psi", "br", "bphi", "bz"]:
            raise ValueError("qnt must be one of the "
                             "following: {\"psi\", \"br\", \"bphi\", \"bz\"}")
        b = self.data.bfield.active.read()

        # Since psi and bfield may have been given in different grids, take
        # explicitly the corresponding grid
        if qnt == "psi":
            name = "psi"
        else:
            name = "b"
        dataphimin = b[name+"_phimin"][0]; dataphimax = b[name+"_phimax"][0]
        datanr, datanphi, datanz = b[qnt].shape
        phi_vec = np.linspace(dataphimin, dataphimax, datanphi+1)[:-1]
        dphi = (dataphimax-dataphimin)/(datanphi)
        datarmin = b[name+"_rmin"][0]; datarmax = b[name+"_rmax"][0]
        datazmin = b[name+"_zmin"][0]; datazmax = b[name+"_zmax"][0]

        # If r or z co-ordinates are given, get qnt at those
        if r is not None:
            rmidindex = int((r-datarmin)/(datarmax-datarmin)*datanr)
        else:
            rmidindex = int(datanr/2)
        # Update r to be the actual r that we are using
        r = np.linspace(datarmin, datarmax, datanr)[rmidindex]

        if z is not None:
            zmidindex = int((z-datazmin)/(datazmax-datazmin)*datanz)
        else:
            zmidindex = int(datanz/2)
        # Update z to be the actual r that we are using
        z = np.linspace(datazmin, datazmax, datanz)[zmidindex]

        qntarray = b[qnt][rmidindex,:,zmidindex]

        # Shift phi_vec values so that the middle of the array correstponds
        # to the "seam" of the input, that is, phi=phimin is put in the middle
        n_elements_rolled = int(datanphi/2)
        phi_vec = phi_vec - n_elements_rolled*dphi
        qntarray = np.roll(qntarray, n_elements_rolled)

        # Select only the centremost points
        frac_points_to_keep = 0.01*centremostpercentage
        start_index = int(datanphi * (1 - frac_points_to_keep) / 2)
        end_index   = datanphi - start_index
        phi_vec = phi_vec[start_index:end_index]
        qntarray = qntarray[start_index:end_index]

        # If diff, plot instead np.diff(qnt)
        if diff:
            qntarray = np.diff(qntarray)   #/np.mean(qntarray)  # normalised diff
            diff_phi = np.diff(phi_vec)
            qntarray = qntarray/diff_phi
            phi_vec = 0.5*(phi_vec[:-1] + phi_vec[1:])
            ylabel = f"diff({qnt:s})/d phi"
        else:
            ylabel = qnt

        # Take only the points closest to the seam
        if n_centremost_elements_for_fit is not None:
            def centremost(arr, n):
                """Return n centremost elements of a 1D numpy array."""
                length = len(arr)
                if n > length:
                    raise ValueError("You requested more elements for the fit" \
                    f" than the input data has inside " \
                    f"centremostpercentage={centremostpercentage:.1f} %. There" \
                    f" are {len(qntarray):d} elements, while you requested" \
                    f" {n:d}. Please, either reduce " \
                    "'n_centremost_elements_for_fit' or increase " \
                    "'centremostpercentage'.")
                start = (length - n) // 2
                end = start + n
                return arr[start:end]
            qntforfit = centremost(qntarray, n_centremost_elements_for_fit)
            phiforfit = centremost(phi_vec, n_centremost_elements_for_fit)
        else:
            qntforfit = qntarray
            phiforfit = phi_vec

        # Fit a polynomial. Third order seems to work well.
        poly_order = 3
        fit = np.polyfit(phiforfit, qntforfit, poly_order)
        fitfunc = np.poly1d(fit)

        # Evaluate the vertical distance between qntarray and the polynomial
        qnt_vertical_distance = np.abs(qntforfit - fitfunc(phiforfit))
        qnt_vertical_distance_avg = np.mean(qnt_vertical_distance)
        qnt_vertical_distance_normalised = qnt_vertical_distance/qnt_vertical_distance_avg

        # If normalised deviation above threshold, the input point may cause a kink
        problemeatic_points_indices = qnt_vertical_distance_normalised > 2.0
        n_problematic_points = np.sum(problemeatic_points_indices)
        if n_problematic_points > 0:
            warnings.warn(f"{n_problematic_points:d} potentially problematic " \
                          f"points in the {qnt:s} input found near the 'seam'" \
                          " in phi direction.")

        # If kink found, plot the figure for the user to double check
        if n_problematic_points or forceplot:
            axes = a5plt.line2d([phi_vec],
                                [qntarray],
                                axes=axes,
                                marker="o",
                                markerfacecolor="none",
                                title=f"R = {r:.2f} m, z = {z:.2f} m",
                                label="input data", skipshow=True)

            # Plot also the spline-interpolated values to assess how well the
            # splines handle the data
            if plot_spline_interpolated:
                phi_high_res = np.linspace(phi_vec[0],phi_vec[-1],100)
                high_res_qnt = self.input_eval(r, phi_high_res, z, 0, qnt)
                if diff:
                    phi_hr_diff = np.diff(phi_high_res)
                    phi_high_res = 0.5*(phi_high_res[:-1] + phi_high_res[1:])
                    high_res_qnt = np.diff(high_res_qnt)
                    high_res_qnt = high_res_qnt/phi_hr_diff
                axes = a5plt.line2d([phi_high_res], [high_res_qnt], axes=axes,
                            label="interpolated", skipshow=True)

            # Plot the polynomial fit to verify that it is sensible
            phi_high_res_centre = np.linspace(phiforfit[0], phiforfit[-1], 50)
            axes = a5plt.line2d([phi_high_res_centre],
                                [fitfunc(phi_high_res_centre)],
                                axes=axes,
                                title=f"R = {r:.2f} m, z = {z:.2f} m",
                                label=f"lsq fit, polynomial order {poly_order:d}",
                                c="k",
                                skipshow=True,
                                )

            # Circle problematic input points with red
            a5plt.scatter2d(phiforfit[problemeatic_points_indices],
                            qntforfit[problemeatic_points_indices],
                            c="r",
                            axes=axes,
                            xlabel="phi [deg]",
                            ylabel=ylabel,
                            marker="o",
                            markersize=20,
                            markerfacecolor="none",
                            label="Potentially problematic points",
                            skipshow=True,
                            )

            a5plt.tight_layout(axes=axes)
            a5plt.legend()
            a5plt.show()


    def input_plotrz(self, r, z, qnt, phi=0*unyt.deg, t=0*unyt.s, rhomax=None,
                     clim=None, cmap=None, axes=None, cax=None,
                     drawcontours=False, contourvalues=None, contourcolors=None,
                     contourlinestyles="-", contourlinewidths=None,
                     labelfontsize=10):
        """Plot input quantity on a (R, z) plane at given coordinates.

        With drawcontours=True, draws also the contourlines.
        To plot quantity on a logarithmic scale (base 10), add "log" in
        the name of the quantity e.g. "log ne".

        Parameters
        ----------
        r : array_like (nr,)
             R abscissa where data is evaluated and plotted.
        z : array_like (nz,)
            z abscissa where data is evaluated and plotted.
        qnt : str
            Name of the plotted quantity.
        phi : float, optional
            Toroidal coordinate where data is evaluated.
        t : float, optional
            Time coordinate where data is evaluated.
        rhomax : float, optional
            Do not plot (R, z) points for which rho > rhomax
        clim : list [float, float], optional
            Minimum and maximum color range values.
        cmap : str, optional
            Colormap.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        drawcontours : boolean, optional
            If true, draws contourlines on top of the heatmap.
        contourvalues : list [float], optional
            Which contour lines to draw. Leave empty for automatic selection
        contourcolors : str or list[str], optional
            Colors of the contour lines. By giving a list, one can match the
            colors with the contourvalues (see above).
        contourlinestyles : str or list[str], optional
            Linestyle of the contour lines -- as the reader might have guessed.
            By giving a list, one can match the styles with the contourvalues (
            see above).
        contourlinewidths : float or list[float], optional
            You can probably guess. By giving a list, one can match the widths
            with the contourvalues (see above).
        labelfontsize : float, optional
            Font size of the labels that indicate contour values.
        """
        logscale = False
        if "log" in qnt:
            qnt = qnt.replace("log", "").strip()
            logscale = True

        try:
            r = r.v
        except AttributeError:
            pass
        try:
            z = z.v
        except AttributeError:
            pass
        out = np.squeeze(self.input_eval(r*unyt.m, phi, z*unyt.m, t, qnt,
                                         grid=True)[:,0,:,0])
        diverging = np.nanmin(out)*np.nanmax(out) < 0 # If data has both -+ vals

        # Handle possible maximum rho value
        if not rhomax==None and rhomax <= 0:
            print(f"Negative rhomax given ({rhomax:.2e}) to input_plotrz(). Ignoring rhomax.")
            rhomax=None
        if not rhomax==None:
            # Throw away values for which rho > rhomax

            # Get rho at these points and set correspondin data to nan
            rho_vec = np.squeeze(self.input_eval(r*unyt.m, phi, z*unyt.m, t,
                                                 "rho", grid=True)[:,0,:,0]).v
            ind_ignore = rho_vec > rhomax
            out[ind_ignore] = np.nan

        a5plt.mesh2d(r, z, out.v, diverging=diverging, logscale=logscale,
                     axesequal=True, xlabel="R [m]", ylabel="z [m]",
                     clabel=qnt + " [" + str(out.units) + "]", clim=clim,
                     cmap=cmap, axes=axes, cax=cax)
        if drawcontours:
            a5plt.contour2d(r, z, out.v, contours=contourvalues, xlabel="R [m]",
                            ylabel="z [m]", axesequal=True,
                            colors=contourcolors,
                            linestyles=contourlinestyles,
                            linewidths=contourlinewidths, axes=axes,
                            contourlabels=True, labelfontsize=labelfontsize)

    @openfigureifnoaxes(projection=None)
    def input_plotrhocontour(self, rho=1.0, phi=0*unyt.deg, t=0*unyt.s,
                             ntheta=180, axes=None, **kwargs):
        """Plot contours of constant rho.

        Calling this without arguments plots the separatrix.

        Parameters
        ----------
        rho : float or array_like, optional
            Value(s) of square root of normalized poloidal flux whose
            contours are plotted.
        phi : float, optional
            Toroidal angle at which contours are plotted.
        t : float, optional
            Time instance when the contours are plotted.
        ntheta : int, optional
            Poloidal resolution of the contours.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        **kwargs
            Arguments passed to :meth:`.matplotlib.pyplot.plot`.
        """
        if not isinstance(rho, list): rho = [rho]
        for rhoval in rho:
            v  = rhoval * np.ones((ntheta,))
            th = np.linspace(0,2*np.pi,ntheta)*unyt.rad
            r, z = self.input_rhotheta2rz(v, th, phi, t)
            axes.plot(r, z, **kwargs)

        axes.set_xlabel("R [m]")
        axes.set_ylabel("z [m]")

    def input_plotphitheta(self, phi, theta, rho, qnt, t=0*unyt.s,
                     clim=None, cmap=None, axes=None, cax=None,
                     drawcontours=False, contourvalues=None, contourcolors=None,
                     contourlinestyles="-", contourlinewidths=None,
                     labelfontsize=10, drawbfielddir=False,
                     bfieldlinecolor="black", bfieldlinedensity=0.5,
                     blinewidth=1.0, blinelength=50.0):
        """Plot input quantity on a (phi, theta) surface at given rho.

        With drawcontours=True, draws also the contourlines.
        To plot quantity on a logarithmic scale (base 10), add "log" in
        the name of the quantity e.g. "log ne".

        Parameters
        ----------
        phi : array_like (nphi,)
            phi abscissa where data is evaluated and plotted.
        theta : array_like (ntheta,)
            theta abscissa where data is evaluated and plotted.
        rho : float
            Rho co-ordinate of the flux surface.
        qnt : str
            Name of the plotted quantity.
        t : float, optional
            Time coordinate where data is evaluated.
        clim : list [float, float], optional
            Minimum and maximum color range values.
        cmap : str, optional
            Colormap.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        drawcontours : boolean, optional
            If true, draws contourlines on top of the heatmap.
        contourvalues : list [float], optional
            Which contour lines to draw. Leave empty for automatic selection
        contourcolors : str or list[str], optional
            Colors of the contour lines. By giving a list, one can match the
            colors with the contourvalues (see above).
        contourlinestyles : str or list[str], optional
            Linestyle of the contour lines -- as the reader might have guessed.
            By giving a list, one can match the styles with the contourvalues (
            see above).
        contourlinewidths : float or list[float], optional
            You can probably guess. By giving a list, one can match the widths
            with the contourvalues (see above).
        labelfontsize : float, optional
            Font size of the labels that indicate contour values.
        drawbfielddir : boolean, optional
            Draw the direction of the bfield on top of the figure as
            streamlines. This makes sense only in a tokamak geometry since the
            poloidal plane is a contant phi plane. If drawfielddir true, one
            must have theta and phi go from 0 to 360.
        bfieldlinecolor : str, optional
            Color of the bfield lines if drawbfielddir=True.
        bfieldlinedensity : float, optional
            Density of the bfield lines.
        blinewidth : float, optional
            Width of the bfield lines.
        blinelength : float, optional
            Length of the bfield lines.
        """
        logscale = False
        if "log" in qnt:
            qnt = qnt.replace("log", "").strip()
            logscale = True

        try:
            phi = phi.v
        except AttributeError:
            pass
        try:
            theta = theta.v
        except AttributeError:
            pass

        # phi and theta inputs are abscissae 1d arrays. We need to make them
        # into a grid and transform them into a pair of 1d arrays containing all
        # grid points where we wnat to plot
        phiphi, thetatheta = np.meshgrid(phi,theta)
        phi_grid1d = phiphi.ravel(); theta_grid1d = thetatheta.ravel()

        # Map the (phi,theta,rho) to (R,z)
        rho_grid1d = np.ones(len(phi_grid1d))*rho
        r_grid1d, z_grid1d = self.input_rhotheta2rz(rho_grid1d,
                                                    theta_grid1d*unyt.deg,
                                                    phi_grid1d*unyt.deg, t)

        # Get the value of the qnt into a long 1d array
        qnt_grid1d = self.input_eval(r_grid1d, phi_grid1d*unyt.deg, z_grid1d,
                                     t, qnt)

        # Go from raveled array to 2d
        qnt_grid2d = qnt_grid1d.reshape(len(phi), len(theta), order="F")

        diverging = np.nanmin(qnt_grid2d)*np.nanmax(qnt_grid2d) < 0 # If data has both -+ vals

        a5plt.mesh2d(phi, theta, qnt_grid2d.v, diverging=diverging, logscale=logscale,
                     axesequal=True, xlabel="phi [deg]", ylabel="theta [deg]",
                     clabel=qnt + " [" + str(qnt_grid2d.units) + "]", clim=clim,
                     cmap=cmap, axes=axes, cax=cax)


        if drawcontours:
            a5plt.contour2d(phi, theta, qnt_grid2d.v, contours=contourvalues,
                            xlabel="phi [deg]", ylabel="theta [deg]",
                            axesequal=True, colors=contourcolors,
                            linestyles=contourlinestyles,
                            linewidths=contourlinewidths, axes=axes,
                            contourlabels=True, labelfontsize=labelfontsize)


        if drawbfielddir:
            # Evaluate the magnitude of the poloidal component
            br, bz, bphi = self.input_eval(r_grid1d, phi_grid1d*unyt.deg,
                                           z_grid1d, t, "br", "bz", "bphi")
            br=br.v; bz=bz.v; bphi=bphi.v
            bpolmag = np.sqrt(br**2 + bz**2)

            # To determine whether the poloidal component rotates clockwise or
            # anti-clockwise, we need to know where the magnetic axis is.
            out = self._eval_bfield(1*unyt.m, phi_grid1d, 1*unyt.m, t,
                                    evalaxis=True)
            axisr = out["axisr"]
            axisz = out["axisz"]

            # We now get the bpol sign from the sign of the cross product
            # "r cross Bpol", where r points from the magnetic axis to the (R,z)
            # location we are considering.
            dr = r_grid1d-axisr     # r distance from the magnetic axis
            dz = z_grid1d-axisz     # z
            bpol_sign = np.sign(dr*bz - dz*br)

            # Get bpol
            bpol = bpol_sign*bpolmag

            # Arrange Bphi and Bpol into 2d arrays
            nphi = len(phi)
            ntheta = len(theta)
            Bphi_2d = bphi.reshape((ntheta, nphi), order='F')
            Btheta_2d = bpol.reshape((ntheta, nphi), order='F')

            # Before making the streamline plot of Bphi and Bpol, we need to
            # normalise them with their respective lengths.

            # Get the poloidal length (circumference) in different phi slices
            R2d = r_grid1d.reshape((ntheta,nphi))
            z2d = z_grid1d.reshape((ntheta,nphi))
            # Distances (each column is some phi value. Rows need to be summed
            # together to get the total arc length on that phi value (column))
            dl_poldir = np.sqrt((R2d[1:,:] - R2d[:-1,:])**2 +
                                (z2d[1:,:] - z2d[:-1,:])**2)
            pol_circumference = np.sum(dl_poldir, axis=0)

            # Get the toroidal length (circumference) in different theta slices
            x_grid1d = np.cos(phi_grid1d*np.pi/180)*r_grid1d
            y_grid1d = np.sin(phi_grid1d*np.pi/180)*r_grid1d
            x2d = x_grid1d.reshape((ntheta,nphi))
            y2d = y_grid1d.reshape((ntheta,nphi))
            dl_tordir = np.sqrt((x2d[:,1:] - x2d[:,:-1])**2 +
                                (y2d[:,1:] - y2d[:,:-1])**2)
            tor_circumference = np.sum(dl_tordir, axis=1)

            # Scale Bphi and Bpol with the respective circumferences
            Bphi_2d_scaled = Bphi_2d/tor_circumference[:,np.newaxis]
            Btheta_2d_scaled = Btheta_2d/pol_circumference[np.newaxis,:]

            axes.streamplot(phi, theta, Bphi_2d_scaled.v, Btheta_2d_scaled.v,
                            color=bfieldlinecolor, density=bfieldlinedensity,
                            linewidth=blinewidth, maxlength=blinelength)


    def input_plotfluxsurface3d(self, rho, phi, theta, qnt=None, t=0*unyt.s,
                     clim=None, cmap=None, axes=None, cax=None, meshalpha = 0.6,
                     tri_lc=None, tri_lw=0.4, opacity=0.8):
        """Plot a flux surface or, given qnt, a quantity on a flux surface in 3D.

        To plot quantity on a logarithmic scale (base 10), add "log" in
        the name of the quantity e.g. "log ne".

        Parameters
        ----------
        rho : float
            Rho co-ordinate of the flux surface.
        phi : array_like (nphi,)
            phi abscissa where data is evaluated and plotted.
        theta : array_like (ntheta,)
            theta abscissa where data is evaluated and plotted.
        qnt : str, optional
            Name of the plotted quantity. E.g. "bnorm"
        t : float, optional
            Time coordinate where data is evaluated.
        clim : list [float, float], optional
            Minimum and maximum color range values.
        cmap : str, optional
            Colormap.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        meshalpha: float, optional
            A parameter used when generating the triangles for the surface. For
            more info, see "alphashape" and "alpha parameter".
        tri_lc : str, optional
            Color of the mesh triangle edges. By default, edges are invisible.
        tri_lw : float, optional
            Width of the mesh triangle edges. To enable edges, set also tri_lc.
        opacity : float, optional
            Opacity of the drawn surface (0: transparent, 1: opaque).
        """
        try:
            phi = phi.v
        except AttributeError:
            pass
        try:
            theta = theta.v
        except AttributeError:
            pass

        # phi and theta inputs are abscissae 1d arrays. We need to make them
        # into a grid and transform them into a pair of 1d arrays containing all
        # grid points where we wnat to plot
        phiphi, thetatheta = np.meshgrid(phi,theta)
        phi_grid1d = phiphi.ravel(); theta_grid1d = thetatheta.ravel()

        # Map the (phi,theta,rho) to (R,z)
        rho_grid1d = np.ones(len(phi_grid1d))*rho
        r_grid1d, z_grid1d = self.input_rhotheta2rz(rho_grid1d,
                                                    theta_grid1d*unyt.deg,
                                                    phi_grid1d*unyt.deg, t)
        # Map (R, phi) to (x,y)
        x_grid1d = np.cos(phi_grid1d*np.pi/180)*r_grid1d
        y_grid1d = np.sin(phi_grid1d*np.pi/180)*r_grid1d

        # Draw surface with or without quantity giving the color
        if not qnt==None:
            # Get the value of the qnt into a long 1d array
            logscale = False
            if "log" in qnt:
                qnt = qnt.replace("log", "").strip()
                logscale = True
            out_grid1d = self.input_eval(r_grid1d, phi_grid1d*unyt.deg, z_grid1d,
                                        t, qnt)
            # If data has both -+ vals
            diverging = np.nanmin(out_grid1d)*np.nanmax(out_grid1d) < 0

            a5plt.surface3d(x_grid1d, y_grid1d, z_grid1d, qnt_grid1d=out_grid1d,
                            diverging=diverging, logscale=logscale,
                            axesequal=True, xlabel="x [m]", ylabel="y [m]",
                            zlabel="z [m]",
                            clabel=qnt + " [" + str(out_grid1d.units) + "]",
                            clim=clim, cmap=cmap,  axes=axes, cax=cax,
                            meshalpha=meshalpha, tri_lc=tri_lc, tri_lw=tri_lw,
                            opacity=opacity)
        else:
            a5plt.surface3d(x_grid1d, y_grid1d, z_grid1d, axesequal=True,
                            xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
                            clim=clim, cmap=cmap,  axes=axes, cax=cax,
                            meshalpha=meshalpha, tri_lc=tri_lc, tri_lw=tri_lw,
                            opacity=opacity)


    def input_plotcontoursurface3d(self, qnt, isovalue, xmin, xmax, ymin, ymax,
                                   zmin, zmax, pointspermeterx=10,
                                   pointspermetery=10, pointspermeterz=10,
                                   t=0*unyt.s, max_rho=1.0, ax=None, c=0.135,
                                   delta=None, tri_lc=None, tri_lw=0.4,
                                   opacity=0.8, reltol=0.01, axscatter=None,
                                   surfcolor="blue"):
        """
        Plot a 3D contour surface of the specified quantity at a given isovalue.

        If there are holes in you surface, try lowering input parameter c or
        increasing reltol.

        Parameters
        ----------
        qnt : str
            Name of the quantity whose contour surface is drawn (e.g., "bnorm").
        isovalue : float
            Value of the quantity at which to generate the isosurface.
        xmin, xmax : float
            Minimum and maximum bounds along the x-axis.
        ymin, ymax : float
            Minimum and maximum bounds along the y-axis.
        zmin, zmax : float
            Minimum and maximum bounds along the z-axis.
        pointspermeterx, pointspermetery, pointspermeterz : int, optional
            Number of grid points along the x, y, and z directions per meter.
            Ideally, use uniform grid spacing in all directions.
        t : float or unyt_quantity, optional
            Time coordinate at which the data is evaluated. Default is 0 seconds.
        max_rho : float, optional
            Maximum `rho` until (surface outside this is not drawn)
        ax : matplotlib.axes.Axes, optional
            Matplotlib 3D axes to plot on. If None, a new figure and axes are
            created.
        c : float, optional
            Scaling parameter for the grid spacing or contour thresholding.
            Default is 0.135.
        delta : float or None, optional
            Characteristic grid spacing. Ideally, use uniform grid of points.
        tri_lc : str or None, optional
            Color of the mesh triangle edges. If None, edges are invisible.
        tri_lw : float, optional
            Width of the mesh triangle edges. Default is 0.4.
        opacity : float, optional
            Opacity of the drawn surface. (0: transparent, 1: opaque)
        reltol : float, optional
            Relative tolerance for isosurface extraction. Default is 0.01.
        axscatter : 3D axes, optional
            If given, plots the grid points of the surface on these axes. This
            can be useful especially when verifying that the the plotted surface
            is sensible.
        """

        nx = max(int(pointspermeterx*(xmax-xmin)), 1)
        ny = max(int(pointspermetery*(ymax-ymin)), 1)
        nz = max(int(pointspermeterz*(zmax-zmin)), 1)

        # Construct a grid of points
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        z = np.linspace(zmin, zmax, nz)
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        x1d = xx.ravel(); y1d = yy.ravel(); z1d = zz.ravel()
        points = np.column_stack((x1d, y1d, z1d))

        # R and phi co-ordinates needed because we call input_eval()
        r1d = np.sqrt(x1d**2 + y1d**2)
        phi1d = np.mod(np.arctan2(y1d, x1d), 2 * np.pi)

        # Evaluate qnt and rho in all grid points
        qnt_grid1d, rho1d = self.input_eval(r1d * unyt.m, phi1d * unyt.rad, z1d * unyt.m, t, qnt, "rho")

        # If not given, calculate delta, the average grid spacing
        if delta==None:
            delta = (xmax-xmin)/nx

        # This parameter specifies the inverse of the size of triangles that are
        # kept. Consequently, it should be inversely proportional to the grid
        # spacing. You don't want to throw away too small triangles because the
        # surface then starts to have holes. On the other hand, you want to
        # throw away the largest triangles because otherwise there will be
        # artificial connections/triangles in places where there should be no
        # surface at all. The severity of this problem greatly depends on the
        # surface. The parameter c has been defined for fine tuning. You may try
        # to adjust it manually to get better results.
        meshalpha=c*1/delta

        #This is the tolerance when discarding points outside the contour
        #surface. If the grid spacing is larger you must have a larger tolerance
        # to avoid holes. Feel free to adjust the reltol if you think you can
        # find a more suitable value.
        tol = reltol*delta

        # Throw away points that are not on the contour surface. Also throw away
        # points outside max_rho
        mask = (
            (rho1d < max_rho)
            & (qnt_grid1d >= (1.0 - tol) * isovalue)
            & (qnt_grid1d <= (1.0 + tol) * isovalue)
        )
        surface_points = points[mask]

        # Plot the surface
        a5plt.surface3d(surface_points[:,0], surface_points[:,1],
                        surface_points[:,2], axesequal=True, xlabel="x [m]",
                        ylabel="y [m]", zlabel="z [m]", axes=ax,
                        meshalpha=meshalpha, tri_lc=tri_lc, tri_lw=tri_lw,
                        opacity=opacity, surfacecolor=surfcolor)

        # Scatter plot for sanity check (you can verify if there are
        # sufficiently many points)
        if not axscatter==None:
            a5plt.scatter3d(surface_points[:,0], surface_points[:,1],
                            surface_points[:,2], xlabel="x [m]", ylabel="y [m]",
                            zlabel="z [m]", markersize=1, axesequal=True, axes=axscatter)

    def get_plasmaquantities(self):
        """Return species present in plasma input.
        """
        spec = self.get_plasmaspecies()
        quantities = []
        for i in range(1,spec["nspecies"]):
            quantities.append("ni" + str(i))
            quantities.append("ti" + str(i))

        return quantities

    def input_rhovolume(self, nrho=10, ntheta=360, nphi=10, method="mc",
                        tol=1e-1, t=0*unyt.s, return_area=False,
                        return_coords=False, minrho=0, maxrho=1, minphi=0, maxphi=360, mintheta=0, maxtheta=360):
        """Evaluate flux surface volumes.

        Parameters
        ----------
        nrho : int, optional
            Number of radial grid edges between [0, 1].
        ntheta : int, optional
            Number of poloidal grid edges.
        nphi : int, optional
            Number of toroidal grid edges.
        method : {"prism", "mc"}, optional
            Method used to evaluate the volumes.

            "prism" divides the volume into structured 3D grid consisting of
            triangular prims where the triangular plane is poloidal and one
            of the vertices is at magnetic axis while the other two are along
            a given flux surface. The volume is the estimated as a sum of these
            prisms.

            "mc" divides the volume into toroidal sectors of same size, and on
            each sector fins the separatrix bounding box. Bounding box is
            assumed to be fixed along the sector and markers are randomly
            sampled within it. The volumes are then estimated via the Monte
            Carlo method.
        tol : float, optional
            If relative difference in volume in subsequent iterations falls
            below this value, the algorithm finishes.
        t : float, optional
            Time slice when the volumes are computed.
        return_area : bool, optional
            Return also the area on (R, z) plane.
        return_coords : bool, optional
            Return also the (R, phi, z) coordinates for grid center points.
        minrho : float, optional
            minimum value of rho
        maxrho : float, optional
            maximum value of rho
        minphi : float, optional
            minimum value of phi (in degrees)
        maxphi : float, optional
            maximum value of phi (in degrees)
        mintheta : float, optional
            minimum value of theta (in degrees)
        maxtheta : float, optional
            maximum value of theta (in degrees)

        Returns
        -------
        volume : array_like (nrho-1, ntheta-1, nphi-1)
            Volume of each bin.
        area : array_like (nrho-1, ntheta-1, nphi-1)
            Area of the poloidal cross section of each bin if
            ``return_area=True``.
        r : array_like (nrho-1, ntheta-1, nphi-1)
            R coordinates of bin centers if ``return_coords=True``.
        p : array_like (nrho-1, ntheta-1, nphi-1)
            phi coordinates of bin centers if ``return_coords=True``.
        z : array_like (nrho-1, ntheta-1, nphi-1)
            z coordinates of bin centers if ``return_coords=True``.
        """
        if nrho < 2 or ntheta < 2 or nphi < 2:
            raise ValueError(
                "Minimum number of edges is 2 but nrho=%d, ntheta=%d, nphi%d"
                + " were given", (nrho, ntheta, nphi))

        # Grid edges and area and volume of the grid cells
        rho    = np.linspace(minrho, maxrho,   nrho) * unyt.dimensionless
        phi    = np.linspace(minphi, maxphi/180*np.pi, nphi) * unyt.rad
        theta  = np.linspace(mintheta, maxtheta/180*np.pi, ntheta) * unyt.rad
        area   = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.m**2
        volume = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.m**3

        # Compute grid center points in (R,phi,z)
        rhoc, thc, p = np.meshgrid(
            ( rho[1:]   + rho[:-1]   ) / 2,
            ( theta[1:] + theta[:-1] ) / 2,
            ( phi[1:]   + phi[:-1]   ) / 2, indexing="ij")
        r, z = self.input_rhotheta2rz(
            rhoc.ravel(), thc.ravel(), p.ravel(), t)
        del thc, rhoc
        r = np.reshape(r, p.shape)
        z = np.reshape(z, p.shape)

        def integrate_prism(rho, theta, phimin, phimax):
            """Integrate volume by dividing it into prisms.
            """
            # Integration variable theta
            th = np.linspace(0, 2*np.pi, 1000)

            # Get axis location
            out = self._eval_bfield(
                1*unyt.m, (phimin + phimax) / 2, 1*unyt.m, t, evalaxis=True)
            r0 = out["axisr"]
            z0 = out["axisz"]

            # Get contour points for each rho value
            rhoc, thc, phic = np.meshgrid(
                rho[1:], th, ( phimin + phimax ) / 2, indexing="ij" )
            rc, zc = self.input_rhotheta2rz(
                rhoc.ravel(), thc.ravel()*unyt.rad, phic.ravel(), t)
            rc = rc.reshape(rhoc.shape)[:,:,0]
            zc = zc.reshape(rhoc.shape)[:,:,0]

            # Calculate the triangle areas and triangle centroid
            area = 0.5 * np.abs(  r0        * ( zc[:,:-1] - zc[:,1:]  )
                                + rc[:,:-1] * ( zc[:,1:]  - z0        )
                                + rc[:,1:]  * ( z0        - zc[:,:-1] ) )
            h       = (phimax - phimin).v
            centerr = (r0 + rc[:,:-1] + rc[:,1:]) / 3
            volume  = h * centerr * area

            # Now divide these to theta sectors that were requested
            V = np.zeros((rho.size-1, theta.size-1)) * unyt.m**3
            A = np.zeros((rho.size-1, theta.size-1)) * unyt.m**2
            th = th[:-1] # Drop last index as it is same as first
            for it in range(theta.size-1):
                idx = np.logical_and(th >= theta[it], th < theta[it+1])
                V[:,it] = np.sum(volume[:, idx], axis=1)
                A[:,it] = np.sum(area[:, idx], axis=1)

            V[1:,:] -= V[:-1,:]
            A[1:,:] -= A[:-1,:]
            return V, A

        def integrate_mc(rho, theta, phimin, phimax):
            """Integrate volume by dividing it into prisms.
            """
            # Determine the bounding box
            th = np.linspace(0, 2*np.pi, 1000)
            rhoc, thc, phic = np.meshgrid(1.0*unyt.dimensionless, th,
                                          unyt.unyt_array([phimin, phimax]),
                                          indexing="ij" )
            rc, zc = self.input_rhotheta2rz(
                rhoc.ravel(), thc.ravel()*unyt.rad, phic.ravel(), t)
            bbox = np.ones((4,)) * unyt.m
            bbox[0] = np.nanmin(rc)
            bbox[1] = np.nanmax(rc)
            bbox[2] = np.nanmin(zc)
            bbox[3] = np.nanmax(zc)

            # Get axis location
            out = self._eval_bfield(
                1*unyt.m, (phimin + phimax) / 2, 1*unyt.m, t, evalaxis=True)
            r0 = out["axisr"]
            z0 = out["axisz"]

            # Number of markers drawn per iteration, total number of markers
            # drawn and area and volume of the test region
            ndraw  = 10000
            ntot   = 0
            atotal = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
            vtotal = 0.5 * (bbox[1]**2 - bbox[0]**2) \
                * ( bbox[3] - bbox[2] ) * ( phimax - phimin ).v
            deltarho = 1.0 / (rho.size-1)
            deltath  = 2*np.pi / (theta.size-1)

            v0     = np.zeros((rho.size-1, theta.size-1)) * unyt.m**3
            area   = np.zeros((rho.size-1, theta.size-1)) * unyt.m**2
            volume = np.zeros((rho.size-1, theta.size-1)) * unyt.m**3
            points_in_bins = np.zeros((rho.size-1, theta.size-1))
            while np.sum(v0) == 0 or \
                  np.max(np.abs( volume[v0>0] - v0[v0>0] ) / v0[v0>0])  > tol:
                v0[:,:] = volume
                pi = np.random.uniform(phimin, phimax,   ndraw) * unyt.rad
                ri = np.random.uniform(bbox[0], bbox[1], ndraw) * unyt.m
                zi = np.random.uniform(bbox[2], bbox[3], ndraw) * unyt.m
                th = np.arctan2(zi-z0,ri-r0) + np.pi

                rhoi    = self.input_eval(ri, pi, zi, t, "rho").v
                ind     = rhoi <= 1 # Reject markers outside separatrix
                i_rho   = np.floor(rhoi[ind] / deltarho).astype(int)
                i_theta = np.floor(th[ind] / deltath).astype(int)

                np.add.at(points_in_bins, (i_rho, i_theta), 1)
                ntot += ndraw

                area[:,:]   = atotal * points_in_bins / ntot
                volume[:,:] = vtotal * points_in_bins / ntot

            return volume, area

        if method == "prism":
            integrate = integrate_prism
        elif method == "mc":
            integrate = integrate_mc

        for ip in range(nphi-1):
            volume[:,:,ip], area[:,:,ip] = \
                integrate(rho, theta, phi[ip], phi[ip+1])
            dphi = (phi[ip+1] - phi[ip]) / 2
            v, _   = integrate(rho, theta, phi[ip], phi[ip] + dphi)
            v1, _  = integrate(rho, theta, phi[ip] + dphi, phi[ip] + dphi*2)
            v += v1
            while(np.abs(np.sum(v) - np.sum(volume[:,:,ip])) > tol):
                volume[:,:,ip] = v
                dphi /= 2
                v[:] = 0 * unyt.m**3
                for i in range(int((phi[ip+1] - phi[ip])/dphi)):
                    v1, _ = integrate(rho, theta, phi[ip] + dphi*i,
                                      phi[ip] + dphi*(i+1))
                    v += v1
            volume[:,:,ip] = v

        if np.sum(volume == 0):
            warnings.warn(
                "Volume calculation resulted in one or more elements with zero"
              + " volume. This might happen when psi0 differs from the actual"
              + " value.")

        if return_area and return_coords:
            return volume, area, r, p, z
        elif return_area:
            return volume, area
        elif return_coords:
            return volume, r, p, z
        else:
            return volume

    def input_eval_safetyfactor(self, rho, nth=10000, return_ftrap=False):
        """Evaluate safety factor and associated quantities.

        This function most likely works for tokamaks only.

        Parameters
        ----------
        rho : array_like, (n,)
            Radial positions where the safety factor is evaluated.

            Note that evaluating it near the axis or exactly at the separatrix
            may yield mess.
        return_bminmax : bool, optional
            Return also minimum and maximum values of B along the flux surface.

        Returns
        -------
        qprof : array_like, (n,)
            Safety factor.
        Iprof : array_like, (n,)
            Toroidal current term which, when multiplied with 2pi/mu0, gives the
            enclosed toroidal current.
        gprof : array_like, (n,)
            This is just R * Bphi which is constant since Bphi ~ 1/R.
        ftrap : array_like, (n,)
            Trapped particle fraction if `return_ftrap` is true.
        """
        ftrap = np.ones(rho.shape)
        qprof = np.zeros(rho.shape)
        Iprof = np.zeros(rho.shape) * unyt.m*unyt.T
        gprof = np.zeros(rho.shape) * unyt.m*unyt.T
        thgrid = np.linspace(0, 2*np.pi, nth)
        for i in range(rho.size):
            # Interpolate the contour points on the fixed theta grid
            r, z = self.input_rhotheta2rz(
                rho[i]*np.ones(thgrid.shape), thgrid*unyt.rad,
                np.zeros(thgrid.shape)*unyt.rad, 0*unyt.s)

            # Magnetic field along the contour (psi can be used to check that
            # the contour was set properly). Drop the last element in r and z
            # as it is the same as first.
            br, bphi, bz, psi = self.input_eval(
                r[:-1], 0*unyt.rad, z[:-1], 0*unyt.s,
                "br", "bphi", "bz", "psi")

            bpol  = np.sqrt(br**2 + bz**2)
            bnorm = np.sqrt(br**2 + bphi**2 + bz**2)
            ds    = (np.diff(r) * br + np.diff(z) * bz) / bpol # darc dot e^_pol
            r = r[:-1]; z = z[:-1]

            # The toroidal current term (multiplying this with mu0/2 pi gets
            # enclosed toroidal current)
            Iprof[i] = np.sum( ds * bpol ) / ( 2*np.pi )

            # g = R*Bphi. Since Bphi ~ 1/R, this is a constant
            gprof[i] = r[0] * bphi[0]

            # The (global) safety factor q(psi)
            qprof[i] = np.sum( ds * gprof[i] / ( r**2 * bpol ) ) / ( 2*np.pi )

            # Trapped particle fraction
            dth = np.diff(thgrid)
            bmax = np.max(bnorm)
            b_favg = np.sum( bnorm * dth / bpol ) / np.sum( dth / bpol )
            dtau = 1.0 / 99
            for tau in np.linspace(0,1,100):
                favg = ( np.sum( np.sqrt(1 - tau * bnorm / bmax) * dth / bpol )
                        / np.sum( dth / bpol ) )
                ftrap[i] -= (3.0/4) * (b_favg / bmax)**2 * tau * dtau / favg

        if return_ftrap:
            return qprof, Iprof, gprof, ftrap
        return qprof, Iprof, gprof

    def _input_eval_orbitboundaries(self, mugrid, ptorgrid, ekin):
        """Plot boundaries of different orbit topologies in (mu, ekin)
        phase-space when energy is fixed.

        WIP

        Parameters
        ----------
        mugrid : array_like
            Grid for the magnetic moment abscissa.
        ptorgrid : array_like
            Grid for the the canonical toroidal angular momentum abscissa.
        ekin : float
            Particle energy.
        """
        pass

    def input_eval_ripple(
            self, rgrid, zgrid, rlarmor, rho0=1.0*unyt.dimensionless,
            theta0=0*unyt.rad, t=0*unyt.s, nphi=362, plot=True, axes1=None,
            axes2=None):
        """Evaluate various ripple quantities.

        The variation of toroidal field strength as a function of toroidal angle
        is called the toroidal field (TF) ripple. The TF ripple can lead to
        ripple-induced losses provided that i) the particle is trapped in
        a ripple well in which case it is promptly lost or ii) the ripple causes
        particle motion to become stochastic which leads to radial diffusion.

        This function evaluates the quantities relevant for ripple transport.
        Note that ``deltacrit`` and ``ripplewell`` are approximate quantities
        based on analytical estimates.

        Parameters
        ----------
        rgrid : array_like
            R points where ``delta``, ``deltacrit``, and ``ripplewell`` are
            evaluated.
        zgrid : array_like
            z points where ``delta``, ``deltacrit``, and ``ripplewell`` are
            evaluated.
        rlarmor : float or tuple (str, int, float)
            Larmor radius of the test particle or tuple with
            (species, charge, energy).
        rho0 : float, optional
            Radial coordinate where ``amplitude`` is evaluated.
        theta0 : float, optional
            Poloidal coordinate where ``amplitude`` is evaluated.
        t : float, optional
            Time instance for time-dependent fields.
        nphi : int, optional
            Number of toroidal grid points in ripple evaluation.

            There is probably never a need to change this number.
        plot : bool, optional
            If True, the evaluated quantities are also plotted.
        axes1 : , optional
            Axes for plotting ``amplitude`` if ``plot`` = True.
        axes2 : , optional
            Axes for plotting ``delta``, ``deltacrit``, and ``ripplewell`` if
            ``plot`` = True.

        Returns
        -------
        phigrid : array_like, (nphi)
            Toroidal grid [-2, 362] where ``amplitude`` is given.

            This grid is not limited to [0,360] as usual so that one can also
            use it to spot any discontinuities in field if the input was poorly
            constructed.
        amplitude : array_like, (nphi,)
            Variation in toroidal field strength, i.e.
            Bphi(phi) - Bphi_axisymmetric.
        delta : array_like, (nr,nz)
            Ripple magnitude defined as  (Bmax-Bmin) / (Bmax+Bmin), where Bmax
            and Bmin are toroidal extrema of Bphi.
        deltacrit : array_like, (nr,nz)
            Critical ripple magnitude divided by ``delta``.

            Particles are subject to stochastic ripple diffusion when
            ``deltacrit`` < 1.
        ripplewell : array_like, (nr,nz)
            Ripple well parameter.

            Particles are bound to be ripple-trapped when ``ripplewell`` < 1.
            However, this does not necessarily mean that the particle is lost as
            it may escape the well when ripple-trapping causes it to drift
            outwards. As particles that are ripple-trapped drift vertically,
            a good indication that a trapped particle will be lost if the ripple
            well region covers everything vertically below/above that point.
        """
        phi = np.linspace(-2, 362, nphi) * unyt.deg
        bphi, rho, drhodpsi = self.input_eval(
            rgrid, phi, zgrid, 0*unyt.s, "bphi", "rho", "rhodpsi", grid=True)
        bphi = np.squeeze(np.abs(bphi))
        bmax = np.nanmax(bphi, axis=1)
        bmin = np.nanmin(bphi, axis=1)
        delta = (bmax - bmin) / (bmax + bmin)
        rho   = np.squeeze(rho[:,0,:,0])
        drhodpsi = np.squeeze(drhodpsi[:,0,:,0])

        if not isinstance(rlarmor, float):
            prt   = physlib.species.species(rlarmor[0], charge=rlarmor[1])
            ekin  = rlarmor[2] * unyt.eV
            bnorm = np.mean(bphi)
            rlarmor = physlib.gyrolength(
                prt["mass"], prt["charge"], ekin, 0.0, bnorm)

        r0, z0 = self.input_rhotheta2rz(rho0*np.ones((nphi,)),
                                        theta0*np.ones((nphi,)), phi, 0*unyt.s)
        amplitude = self.input_eval(r0, phi, z0, 0*unyt.s, "bphi")
        meanbphi  = np.mean(amplitude)
        amplitude -= meanbphi

        # Use fft to find ripple periodicity
        ff      = np.fft.fftfreq(len(phi), (phi[1]-phi[0]))
        Fyy     = abs(np.fft.fft(amplitude))
        Nripple = np.round(abs(ff[np.argmax(Fyy[1:])+1]) * (phi[-1] - phi[0]))

        rhogrid = np.linspace(0, 1, 100)
        q, _, _ = self.input_eval_safetyfactor(rhogrid)
        q = scipy.interpolate.CubicSpline(rhogrid, q, extrapolate=False)

        raxis, zaxis = self.input_rhotheta2rz(
            0.0, 0*unyt.rad, 0*unyt.rad, 0*unyt.s)

        rg, zg = np.meshgrid(rgrid, zgrid, indexing="ij")
        rminor = np.sqrt( (rg - raxis)**2 + (zg - zaxis)**2 )
        sinth  = (zg - zaxis) / rminor
        eps    = rminor / raxis
        qfac   = np.abs(q(rho))
        dqdpsi = np.abs(q(rho,1) * drhodpsi)

        deltacrit = np.power( eps / (np.pi * Nripple * qfac), 3.0/2 ) / \
            ( rlarmor * dqdpsi * delta )
        ripplewell = rg * np.abs(sinth) / ( Nripple * qfac * delta )

        if plot:
            newfig = False
            if axes1 is None or axes2 is None:
                newfig = True
                fig = plt.figure(figsize=(10,10))
                if axes1 is None and axes2 is None:
                    axes1 = fig.add_subplot(2,1,1)
                    axes2 = fig.add_subplot(2,1,2)
                else:
                    raise ValueError("Please provide both axes for plotting")
            axes1.plot(phi, meanbphi + amplitude)
            axes2.contourf(rgrid, zgrid, deltacrit.T.v, [0.0, 1.0],
                           colors="C0", alpha=0.5)
            axes2.contourf(rgrid, zgrid, ripplewell.T.v, [0.0, 1.0],
                           colors="C3",alpha=0.5)
            cs = axes2.contour(rgrid, zgrid, delta.T*100,
                               [0.1, 0.3, 0.5, 1.0, 5.0, 10.0], colors="black")

            def fmt(x):
                s = f"{x:.1f}"
                if s.endswith("0"):
                    s = f"{x:.0f}"
                return rf"{s} %"
            axes2.clabel(cs, cs.levels,fmt=fmt)

            axes2.set_aspect("equal", adjustable="box")

            axes1.set_xlabel("Toroidal angle [deg]")
            axes1.set_ylabel("Toroidal field [T]")

            axes2.set_xlabel("R [m]")
            axes2.set_ylabel("z [m]")

            if newfig: plt.show()

        return phi, amplitude, delta, deltacrit, ripplewell

    @openfigureifnoaxes()
    def input_eval_orbitresonance(
            self, rhogrid, xigrid, egrid, species, plot=False, n=1, p=0,
            omega=0.0, axes=None, cax=None):
        """Evaluate the resonance condition between particles and perturbations.

        Orbit resonance refers to the phenomenon where the frequency of
        a particle's orbit in a magnetic field matches the frequency of
        a perturbation, such as a wave or an external perturbation. This
        resonance can lead to increased particle transport.

        The resonance condition for orbit resonance is given by:

        f_{tor} / f_{pol} = m / n

        where:
        - f_{tor} is the toroidal orbit frequency of the particle,
        - f_{pol} is the poloidal orbit frequency of the particle,
        - n is the toroidal mode number of the perturbation.
        - m is an arbitrary integer.

        This function can evaluate and plot the resonance in 1D or 2D space.
        Note that one of the ('rhogrid', 'xigrid', 'egrid') must be a scalar.

        Parameters
        ----------
        rhogrid : array_like, (nrho)
            Value of the radial coordinate of radial grid for which
            the resonance is evaluated.
        xigrid : float or array_like, (nxi)
            Value of the pitch of pitch grid for which the resonance is
            evaluated.
        egrid : float or array_like, (ne)
            Value of the energy of energy grid for which the resonance is
            evaluated.
        species : str
            Name of the particle species.
        plot : bool, optional
            If True, the results are not returned but visualized instead.
        resonance : float or array_like, optional
            Resonance(s) to be indicated on the plot e.g. m/n=2/1.
        omega : float, optional
            Perturbation mode frequency.

            The default value assumes stationary perturbation.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.

        Returns
        -------
        torfreq : array_like (nrho, nxi) or (nrho, ne) or (nxi, ne)
            Evaluated toroidal orbit frequencies if 'plot' = False.
        polfreq : array_like (nrho, nxi) or (nrho, ne) or (nxi, ne)
            Evaluated poloidal orbit frequencies if 'plot' = False.
        qeff : array_like (nrho, nxi) or (nrho, ne) or (nxi, ne)
            Evaluated effective qfactor = torfreq / polfreq if 'plot' = False.
        """
        rhogrid = np.atleast_1d(rhogrid)
        xigrid  = np.atleast_1d(xigrid)
        egrid   = np.atleast_1d(egrid)

        if rhogrid.size > 1 and xigrid.size > 1 and egrid.size > 1:
            raise ValueError("One of the input grids must be a scalar")

        mrk = self.data.create_input(
            "marker resonance", species=species, rhogrid=rhogrid, xigrid=xigrid,
            egrid=egrid, dryrun=True)
        opt = self.data.create_input(
            "options singleorbit", ntor=12, npol=6, dryrun=True)

        self.simulation_initmarkers(**mrk)
        self.simulation_initoptions(**opt)
        vrun = self.simulation_run(printsummary=False)

        t_phi0, z_phi0, pitch_phi0, ids_phi0 = vrun.getorbit(
            'mileage', 'z', 'pitch', 'ids', pncrid=0)
        t_omp,  z_omp,  pitch_omp, ids_omp, phi_omp  = vrun.getorbit(
            'mileage', 'z', 'pitch', 'ids', 'phi', pncrid=1)
        self.simulation_free(markers=True, diagnostics=True)

        if rhogrid.size == 1:
            n1 = xigrid.size
            n2 = egrid.size
        elif xigrid.size == 1:
            n1 = rhogrid.size
            n2 = egrid.size
        elif egrid.size == 1:
            n1 = rhogrid.size
            n2 = xigrid.size

        nmrk = mrk['n']
        passing = np.zeros((n1, n2))
        torfreq = np.zeros((n1, n2)) + np.nan
        polfreq = np.zeros((n1, n2)) + np.nan
        for i in range(nmrk):
            ixi  = int(i / n1)
            irho = i % n1
            idx_omp  = ids_omp  == i+1
            idx_phi0 = ids_phi0 == i+1

            # Is this particle passing or trapped? (Does omp crossings have
            # different sign of pitch)
            passing[irho, ixi] = np.all(pitch_omp[idx_omp] < 0) or \
                                 np.all(pitch_omp[idx_omp] > 0)

            if passing[irho, ixi]:
                # Passing
                torfreq[irho, ixi] = 1.0 / np.mean(np.diff(t_phi0[idx_phi0]))
                polfreq[irho, ixi] = 1.0 / np.mean(np.diff(t_omp[idx_omp]))
            else:
                #Trapped
                deltat = np.diff(t_omp[idx_omp][::2])
                polfreq[irho, ixi] = 1.0 / np.mean(deltat)
                deltaphi = np.abs(np.diff(phi_omp[idx_omp][::2]))
                torfreq[irho, ixi] = 1.0 / np.mean(deltat * 360 / deltaphi)

        if not plot:
            return torfreq, polfreq, torfreq/polfreq

        # The rest is plotting
        omega = np.abs(torfreq * n - p * polfreq)
        def plot2d(x, y, title, xlabel, ylabel):
            im   = axes.pcolormesh(x, y, np.log10(1.0/omega).T)
            axes.set_title(title)
            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

            cbar = plt.colorbar(im, ax=axes, orientation='horizontal')
            cbar.set_label(
                r"$log_{10} | \omega_{tor} n - p \omega_{pol} -\omega |$")

        axes.set_xlabel(r"$\rho$")
        x = None; y = None; xlabel = None; ylabel = None
        if rhogrid.size > 1 and xigrid.size > 1:
            plot2d(rhogrid, xigrid, r"$E=%g$ MeV" % (egrid/1e6),
                   r"$\rho$", r"$v_\parallel/v$")
        elif rhogrid.size > 1 and egrid.size > 1:
            plot2d(rhogrid, egrid, r"$v_\parallel/v=%f$" % (xigrid),
                   r"$\rho$", r"$E_{kin}$ [eV]")
        elif xigrid.size > 1 and egrid.size > 1:
            plot2d(xigrid, egrid, r"$\rho=%f$" % (rhogrid),
                   r"$v_\parallel/v$", r"$E$ [MeV]")
        elif rhogrid.size > 1:
            axes.plot(rhogrid, torfreq / polfreq)
            axes.set_ylabel(r"$f_{tor} / f_{pol}$")
