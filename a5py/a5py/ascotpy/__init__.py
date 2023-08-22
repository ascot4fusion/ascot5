"""Python interface to libascot.so.
"""
import ctypes
import unyt
import numpy as np

import a5py.routines.plotting as a5plt
from a5py.physlib import parseunits
from a5py.routines.plotting import openfigureifnoaxes
from a5py.exceptions import AscotInitException

from .libascot    import LibAscot, _LIBASCOT
from .libsimulate import LibSimulate
from .libproviders import LibProviders

if _LIBASCOT:
    from . import ascot2py

class Ascotpy(LibAscot, LibSimulate, LibProviders):
    """Class with methods to initialize and access the data via Python.

    This class will be inherited by Ascot class. Here we hide all the dirty
    implementation of the libascot.so interface in private methods and
    attributes so that the Ascot front is clean for the users.

    MPI parallelism has not been tested.

    Attributes
    ----------
    _sim
        Simulation offload data struct.
    _offload_ready
        Flag indicating if inputs are packed.
    _nmrk
        Number of markers currently in the marker array.
    _inistate
        Marker input array for interactive simulations.
    _endstate
        Marker output array for interactive simulations.
    _diag_occupied
        Flag indicating if diagnostics array contains data.
    _diag_offload_array
        Diagnostics data offload array.
    _offload_array
        Offload array containing the input data when this instance is packed.
    _int_offload_array
        Offload array for integers containing the input data when this instance
        is packed.
    _bfield_offload_array
        Offload data for the magnetic field input.
    _efield_offload_array
        Offload data for the electric field input.
    _plasma_offload_array
        Offload data for the plasma input.
    _neutral_offload_array
        Offload data for the neutral input.
    _wall_offload_array
        Offload data for the wall (float) input.
    _wall_int_offload_array
        Offload data for the wall (int) input.
    _boozer_offload_array
        Offload data for the Boozer input.
    _mhd_offload_array
        Offload data for the MHD input.
    _mpi_root
        Rank of the root MPI process.
    _mpi_rank
        Rank of this MPI process.
    _mpi_size
        Number of MPI processes.
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
        self._offload_ready    = False
        self._nmrk             = ctypes.c_int32()
        self._diag_occupied    = False
        self._offload_data      = ascot2py.struct_c__SA_offload_package()
        self._offload_array     = ctypes.POINTER(ctypes.c_double)()
        self._int_offload_array = ctypes.POINTER(ctypes.c_int   )()
        self._inistate = ctypes.POINTER(ascot2py.struct_c__SA_particle_state)()
        self._endstate = ctypes.POINTER(ascot2py.struct_c__SA_particle_state)()

        self._sim = ascot2py.struct_c__SA_sim_offload_data()
        self._bfield_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self._efield_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self._plasma_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self._neutral_offload_array = ctypes.POINTER(ctypes.c_double)()
        self._wall_offload_array    = ctypes.POINTER(ctypes.c_double)()
        self._boozer_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self._mhd_offload_array     = ctypes.POINTER(ctypes.c_double)()
        self._asigma_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self._diag_offload_array    = ctypes.POINTER(ctypes.c_double)()

        self._wall_int_offload_array = ctypes.POINTER(ctypes.c_int)()

        # MPI init needs to be called (only) once so that we can run simulations
        argc       = ctypes.c_int32()
        argc.value = 0
        argv       = ctypes.POINTER(ctypes.c_char)()
        self._mpi_rank = ctypes.c_int32()
        self._mpi_root = ctypes.c_int32()
        self._mpi_size = ctypes.c_int32()
        ascot2py.mpi_interface_init(
            argc, ctypes.byref(argv), ctypes.byref(self._sim),
            ctypes.byref(self._mpi_rank), ctypes.byref(self._mpi_size),
            ctypes.byref(self._mpi_root))

    def _init(self, data, bfield=None, efield=None, plasma=None,
              wall=None, neutral=None, boozer=None, mhd=None,
              switch=False):
        """Read, offload, and initialize input data so it can be accessed
        by libascot.

        Assumes data is present in the hdf5 file and the provided QIDs are
        valid. If the parameter is given as a dict in same format as it is
        written, it is used directly instead of reading the from the file.

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
        switch : bool
            If ``True``, free input that has been
        """
        if self._offload_ready:
            raise AscotInitException("This instance has been packed")

        # Iterate through all inputs and mark those that are initialized
        inputs2read = ctypes.c_int32()
        args = locals() # Contains function arguments and values in a dictionary

        # List here dependencies to be directly injected (provided)
        to_be_provided = []
        for inp in ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd"]:
            if args[inp] is None:
                # This input is not initialized
                continue

            # Inject/provide this input
            if isinstance(args[inp],dict):
                to_be_provided.append(inp)
                continue

            # Convert QID strings to bytes
            args[inp] = args[inp].encode("UTF-8")

            currentqid = getattr(self._sim, "qid_" + inp)
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
                    "mhd"]:
            if inputs2read.value & getattr(ascot2py, "hdf5_input_" + inp):
                setattr(self._sim, "qid_" + inp, args[inp])

        ascot2py.hdf5_interface_read_input(
            ctypes.byref(self._sim),
            inputs2read,
            ctypes.byref(self._bfield_offload_array),
            ctypes.byref(self._efield_offload_array),
            ctypes.byref(self._plasma_offload_array),
            ctypes.byref(self._neutral_offload_array),
            ctypes.byref(self._wall_offload_array),
            ctypes.byref(self._wall_int_offload_array),
            ctypes.byref(self._boozer_offload_array),
            ctypes.byref(self._mhd_offload_array),
            ctypes.byref(self._asigma_offload_array),
            None, # Marker array (ignore)
            None  # Number of markers that were read (ignore)
            )

        for inp in to_be_provided:
            if inp == "wall":
                if 'x1x2x3' in args[inp]:
                    self.provide_wall_3d(  args[inp]['x1x2x3'],  args[inp]['y1y2y3'],  args[inp]['z1z2z3']  )
                elif 'r' in args[inp]:
                    self.provide_wall_2d(  args[inp]['r'],  args[inp]['z']  )
                else:
                    raise AscotInitException("Unsupported dict for input '{}' passed for injection.".format(inp))                    
            elif inp == "bfield":
                if 'b_rmin' in args[inp]:
                    self.provide_BSTS( b_rmin=args[inp]['b_rmin'], b_rmax=args[inp]['b_rmax'], b_nr=args[inp]['b_nr'],
                                       b_zmin=args[inp]['b_zmin'], b_zmax=args[inp]['b_zmax'], b_nz=args[inp]['b_nz'],
                                       b_phimin=args[inp]['b_phimin'], b_phimax=args[inp]['b_phimax'], b_nphi=args[inp]['b_nphi'],
                                       psi0=args[inp]['psi0'], psi1=args[inp]['psi1'],
                                       br=args[inp]['br'], bphi=args[inp]['bphi'], bz=args[inp]['bz'], psi=args[inp]['psi'],
                                       axis_phimin=args[inp]['axis_phimin'], axis_phimax=args[inp]['axis_phimax'],
                                       axis_nphi=args[inp]['axis_nphi'],
                                       axisr=args[inp]['axisr'], axisz=args[inp]['axisz'],
                                       psi_rmin=args[inp]['psi_rmin'],     psi_rmax=args[inp]['psi_rmax'], psi_nr=args[inp]['psi_nr'],
                                       psi_zmin=args[inp]['psi_zmin'],     psi_zmax=args[inp]['psi_zmax'], psi_nz=args[inp]['psi_nz'],
                                       psi_phimin=args[inp]['psi_phimin'], psi_phimax=args[inp]['psi_phimax'], psi_nphi=args[inp]['psi_nphi'] )
                else:
                    raise AscotInitException("Unsupported dict for input '{}' passed for injection.".format(inp))
            else:
                raise AscotInitException("Unsupported input to inject: '{}'".format(inp))

    def _free(self, bfield=False, efield=False, plasma=False, wall=False,
              neutral=False, boozer=False, mhd=False):
        """Free input data initialized in C-side.
        """
        if self._offload_ready:
            raise AscotInitException("This instance has been packed")

        args = locals() # Contains function arguments and values in a dictionary

        # Iterate through all inputs and free the data if the corresponding
        # argument is True
        for inp in ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd"]:
            if args[inp] and \
               getattr(self._sim, "qid_" + inp) != Ascotpy.DUMMY_QID:

                # Deallocate the data allocated in C side
                array = getattr(self, "_" + inp + "_offload_array")
                ascot2py.libascot_deallocate(array)
                if inp == "wall":
                    array = getattr(self, "_" + inp + "_int_offload_array")
                    ascot2py.libascot_deallocate(array)

                # Set the offload array length to zero
                if inp == "bfield":
                    self._sim.B_offload_data.offload_array_length = 0
                elif inp == "efield":
                    self._sim.E_offload_data.offload_array_length = 0
                else:
                    data = getattr(self._sim, inp + "_offload_data")
                    data.offload_array_length = 0
                    if inp == "wall":
                        data.int_offload_array_length = 0

                # Set QID to dummy value
                setattr(self._sim, "qid_" + inp, Ascotpy.DUMMY_QID)

    def _pack(self):
        """Pack offload arrays as one making this instance ready for simulation.

        Note that inputs cannot be changed or freed before calling _unpack. Make
        sure all required data is initialized before packing.
        """
        if self._offload_ready:
            raise AscotInitException("This instance is already packed")

        # This call internally frees individual offload arrays and initializes
        # the common ones.
        ascot2py.pack_offload_array(
            ctypes.byref(self._sim), ctypes.byref(self._offload_data),
            self._bfield_offload_array, self._efield_offload_array,
            self._plasma_offload_array, self._neutral_offload_array,
            self._wall_offload_array,   self._wall_int_offload_array,
            self._boozer_offload_array, self._mhd_offload_array,
            self._asigma_offload_array,
            ctypes.byref(self._offload_array),
            ctypes.byref(self._int_offload_array))
        self._offload_ready = True

        # Set pointers to correct locations (based on the order arrays are
        # packed) so that we can continue using evaluation routines.
        def advance(prevptr, increment):
            ptr = ctypes.addressof(prevptr) + 2 * increment
            ptr = ctypes.cast(ctypes.c_void_p(ptr),
                              ctypes.POINTER(ctypes.c_double))
            return ptr

        self._bfield_offload_array = self._offload_array
        self._efield_offload_array = advance(
            self._bfield_offload_array.contents,
            self._sim.B_offload_data.offload_array_length)
        self._plasma_offload_array = advance(
            self._efield_offload_array.contents,
            self._sim.E_offload_data.offload_array_length)
        self._neutral_offload_array = advance(
            self._plasma_offload_array.contents,
            self._sim.plasma_offload_data.offload_array_length)
        self._wall_offload_array = advance(
            self._neutral_offload_array.contents,
            self._sim.neutral_offload_data.offload_array_length)
        self._boozer_offload_array = advance(
            self._wall_offload_array.contents,
            self._sim.wall_offload_data.offload_array_length)
        self._mhd_offload_array = advance(
            self._boozer_offload_array.contents,
            self._sim.boozer_offload_data.offload_array_length)

        self._wall_int_offload_array = self._int_offload_array


    def _unpack(self, bfield=True, efield=True, plasma=True, wall=True,
               neutral=True, boozer=True, mhd=True):
        """Unpack simulation arrays, i.e. free offload array and re-read data.

        After unpacking the inputs can be changed or freed again but simulations
        cannot be performed.
        """
        if not self._offload_ready:
            raise AscotInitException("This instance hasn't been packed")

        ascot2py.libascot_deallocate(self._offload_array)
        ascot2py.libascot_deallocate(self._int_offload_array)
        self._offload_data.offload_array_length     = 0
        self._offload_data.int_offload_array_length = 0
        self._offload_data.unpack_pos               = 0
        self._offload_data.int_unpack_pos           = 0

        inputs2read = {}
        def readornot(name, read, qid):
            if read:
                inputs2read[name] = qid.decode("utf-8")
            return Ascotpy.DUMMY_QID

        self._sim.qid_bfield  = readornot("bfield",  bfield,
                                          self._sim.qid_bfield)
        self._sim.qid_efield  = readornot("efield",  efield,
                                          self._sim.qid_efield)
        self._sim.qid_plasma  = readornot("plasma",  plasma,
                                          self._sim.qid_plasma)
        self._sim.qid_neutral = readornot("neutral", neutral,
                                          self._sim.qid_neutral)
        self._sim.qid_wall    = readornot("wall",    wall,
                                          self._sim.qid_wall)
        self._sim.qid_boozer  = readornot("boozer",  boozer,
                                          self._sim.qid_boozer)
        self._sim.qid_mhd     = readornot("mhd",     mhd,
                                          self._sim.qid_mhd)

        self._offload_ready = False
        self.input_init(**inputs2read)

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
                    "mhd"]:
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
            "bnorm":
            "Magnetic field norm (not including MHD)",
            "brdr":
            "d/dR of magnetic field R component (not including MHD)",
            "brdphi":
            "d/dphi of magnetic field R component (not including MHD)",
            "brdz":
            "d/dz of magnetic field R component (not including MHD)",
            "bphidr":
            "d/dR of magnetic field phi component (not including MHD)",
            "bphidphi":
            "d/dphi of magnetic field phi component (not including MHD)",
            "bphidz":
            "d/dz of magnetic field phi component (not including MHD)",
            "bzdr":
            "d/dR of magnetic field z component (not including MHD)",
            "bzdphi":
            "d/dphi of magnetic field z component (not including MHD)",
            "bzdz":
            "d/dz of magnetic field z component (not including MHD)",
            "divergence":
            "Magnetic field divergence (not including MHD)",
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

        nion = self.input_getplasmaspecies()[0]
        for i in range(1, nion+1):
            out["ni" + str(i)] = "Ion species (anum, znum) = () density"
            out["ti" + str(i)] = "Ion species (anum, znum) = () temperature"

        if show:
            for name, desc in out.items():
                print(name.ljust(15) + " : " + desc)

        return out

    @parseunits(r="m", phi="deg", z="m", t="s")
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
            r, phi, z, t = np.meshgrid(r, phi, z, t)
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
        if any(q in qnt for q in ["rho", "psi"]):
            out.update(**self._eval_bfield(r, phi, z, t, evalrho=True))
        if any(q in qnt for q in ["br", "bphi", "bz", "brdr", "brdphi", "brdz",
                                  "bphidr", "bphidphi", "bphidz", "bzdr",
                                  "bzdphi", "bzdz", "divergence", "bnorm",
                                  "jnorm", "jr", "jphi", "jz"]):
            out.update(self._eval_bfield(r, phi, z, t, evalb=True))
            out["divergence"] = out["br"]/r + out["brdr"] + out["bphidphi"]/r \
                + out["bzdz"]
            out["bnorm"] = np.sqrt(out["br"]**2 + out["bphi"]**2 + out["bz"]**2)
            out["jr"]    = (out["bzdphi"]/r - out["bphidz"]) / unyt.mu_0
            out["jphi"]  = (out["brdz"] - out["bzdr"]) / unyt.mu_0
            out["jz"]    = (out["bphi"]/r + out["bphidr"] - out["brdphi"]/r) \
                            / unyt.mu_0
            out["jnorm"] = np.sqrt(out["jr"]**2 + out["jphi"]**2 + out["jz"]**2)
        if any(q in qnt for q in ["axisr", "axisz"]):
            out.update(self._eval_bfield(r, phi, z, t, evalaxis=True))
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
        if any(q in qnt for q in ["ne", "te"] + ni + ti):
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

    def input_plotrz(self, r, z, qnt, phi=0*unyt.deg, t=0*unyt.s,
                     clim=[None, None], cmap=None, axes=None, cax=None):
        """Plot input quantity on a (R, z) plane at given coordinates.

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
        clim : list [float, float], optional
            Minimum and maximum color range values.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        """
        log = False
        if "log" in qnt:
            qnt = qnt.replace("log", "").strip()
            log = True

        out = np.squeeze(self.input_eval(r, phi, z, t, qnt, grid=True)[:,0,:,0])
        diverging = np.nanmin(out)*np.nanmax(out) < 0 # If data has both -+ vals
        try:
            r = r.v
        except AttributeError:
            pass
        try:
            z = z.v
        except AttributeError:
            pass
        a5plt.mesh2d(r, z, out.v, diverging=diverging, log=log,
                     axesequal=True, xlabel="R [m]", ylabel="z [m]",
                     clabel=qnt + " [" + str(out.units) + "]", clim=clim,
                     cmap=cmap, axes=axes, cax=cax)

    def input_plotseparatrix(self, phi=0*unyt.deg, t=0*unyt.s,
                             rlim=[0.1*unyt.m, 20*unyt.m],
                             zlim=[-10*unyt.m, 10*unyt.m], axes=None):
        """Plot separatrix on (R,z) plane.
        """
        r = np.linspace(rlim[0], rlim[1], 100)
        z = np.linspace(zlim[0], zlim[1], 100)
        v = np.squeeze(self.input_eval(r, phi, z, t, "rho", grid=True)[:,0,:,0])
        a5plt.contour2d(r, z, v, [1], xlabel="R [m]", ylabel="z [m]",
                        colors="black", linestyles="solid", linewidths=1,
                        axesequal=True, axes=axes)

    def get_plasmaquantities(self):
        """Return species present in plasma input.
        """
        spec = self.get_plasmaspecies()
        for i in range(1,spec["nspecies"]):
            quantities.append("ni" + str(i))
            quantities.append("ti" + str(i))

        return quantities

    def input_rhovolume(self, nrho=10, method="int", tol=1e-2, ntheta=360,
                        nphi=360, t=0*unyt.s, nperiod=1,
                        return_area=False, return_coords=False):
        """Evaluate flux surface volumes.

        Parameters
        ----------
        nrho : int, optional
            Number of radial grid edges between [0, 1].

        method : {"int", "mc"}, optional
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
            below this value, the Monte Carlo algorithm in "mc" finishes.
        ntheta : int, optional
            Number of poloidal grid edges.
        nphi : int, optional
            Number of toroidal grid edges.
        t : float, optional
            Time slice when the volumes are computed.
        nperiod : int, optional
            How many periods toroidal volume is assumed to have.
        return_area : bool, optional
            Return also the area on (R, z) plane.
        return_coords : bool, optional
            Return also the (R, phi, z) coordinates for grid center points.

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
        rho   = np.linspace(0, 1,   nrho) * unyt.dimensionless
        phi   = np.linspace(0, 2*np.pi / nperiod, nphi) * unyt.rad
        theta = np.linspace(0, 2*np.pi, ntheta) * unyt.rad

        r = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.m
        p = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.rad
        z = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.m

        area   = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.m**2
        volume = np.zeros((nrho-1, ntheta-1, nphi-1)) * unyt.m**3

        if method == "prism":
            for ip in range(nphi-1):
                # Get axis location
                out = self._eval_bfield(1*unyt.m, phi[ip], 1*unyt.m, t,
                                        evalaxis=True)
                r0 = out["axisr"]
                z0 = out["axisz"]
                cr = np.zeros((nrho-1, ntheta-1)) * unyt.m

                # Contour coordinates (excluding rho=0 which is on axis)
                rc0, zc0 = self.input_rhotheta2rz(rho[1:], theta[0], phi[ip], t)
                for it in range(ntheta-1):
                    # (R,z) values corresponding to rho and current theta & phi
                    r[:,it,ip], z[:,it,ip] = self.input_rhotheta2rz(
                        (rho[1:] + rho[:-1]) / 2, (theta[it+1] + theta[it]) / 2,
                        (phi[ip+1] + phi[ip]) / 2, t)
                    p[:,:,ip] = (phi[ip+1] + phi[ip]) / 2

                    rc, zc = self.input_rhotheta2rz(
                        rho[1:], theta[it+1], phi[ip], t)
                    # Calculate the triangle areas and triangle centroid
                    cr[:,it] = (r0 + rc + rc0) / 3
                    area[:,it,ip] = 0.5 * np.abs(  r0*(zc0 - zc) + rc0*(zc - z0)
                                                 + rc*(z0 - zc0) )

                    # Save points so we can calculate the area at the next slice
                    rc0 = rc
                    zc0 = zc

                h = (phi[ip+1] - phi[ip]).v
                volume[:,:,ip] = h * cr * area[:,:,ip]

                # Calculated areas and volumes are for up to given rho value but
                # we ant them to be between two rho values
                area[1:,:,ip]   -= area[:-1,:,ip]
                volume[1:,:,ip] -= volume[:-1,:,ip]

        elif method == "mc":
            for ip in range(nphi-1):
                out = self._eval_bfield(0*unyt.m, phi[ip], 0*unyt.m, t,
                                        evalaxis=True)
                r0 = out["axisr"]
                z0 = out["axisz"]
                bbox = np.array([r0.v, r0.v, z0.v, z0.v]) * unyt.m

                for it in range(ntheta-1):
                    # (R,z) values corresponding to rho and current theta & phi
                    r[:,it,ip], z[:,it,ip] = self.input_rhotheta2rz(
                        (rho[1:] + rho[:-1]) / 2, (theta[it+1] + theta[it]) / 2,
                        (phi[ip+1] + phi[ip]) / 2, t)
                    p[:,:,ip] = (phi[ip+1] + phi[ip]) / 2

                    # Determine bounding box for the separatrix
                    rc, zc = self.input_rhotheta2rz(
                        rho[-1], theta[it], phi[ip], t)
                    bbox[0] = np.nanmin([np.nanmin(rc), bbox[0]])
                    bbox[1] = np.nanmax([np.nanmax(rc), bbox[1]])
                    bbox[2] = np.nanmin([np.nanmin(zc), bbox[2]])
                    bbox[3] = np.nanmax([np.nanmax(zc), bbox[3]])

                # Number of markers drawn per iteration, total number of markers
                # drawn and area and volume of the test region
                ndraw  = 10000
                ntot   = 0
                atotal = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
                vtotal = 0.5 * (bbox[1]**2 - bbox[0]**2) \
                    * ( bbox[3] - bbox[2] ) * ( phi[ip+1] - phi[ip] ).v
                deltarho = 1.0 / (nrho-1)
                deltath  = 2*np.pi / (ntheta-1)

                v0 = volume[:,:,ip].copy()
                points_in_bins = np.zeros((nrho-1, ntheta-1))
                while np.sum(v0) == 0 or \
                      np.nanmax(np.abs( volume[:,:,ip] - v0 ) / v0)  > tol:

                    v0[:,:] = volume[:,:,ip]
                    pi = np.random.uniform(phi[ip], phi[ip+1], ndraw) * unyt.rad
                    ri = np.random.uniform(bbox[0], bbox[1], ndraw) * unyt.m
                    zi = np.random.uniform(bbox[2], bbox[3], ndraw) * unyt.m
                    theta = np.arctan2(zi-z0,ri-r0) + np.pi

                    rho   = self.input_eval(ri, pi, zi, t, "rho").v
                    ind   = rho <= 1 # Reject markers outside separatrix
                    i_rho = np.floor(rho[ind] / deltarho).astype(int)
                    i_theta = np.floor(theta[ind] / deltath).astype(int)

                    np.add.at(points_in_bins, (i_rho, i_theta), 1)
                    ntot += ndraw

                    area[:,:,ip] = atotal * points_in_bins / ntot
                    volume[:,:,ip] = vtotal * points_in_bins / ntot

        if return_area and return_coords:
            return volume, area, r, p, z
        elif return_area:
            return volume, area
        elif return_coords:
            return volume, r, p, z
        else:
            return volume

    def _evaluateripple(self, R, z, t, nphi):
        """TBD
        """
        nin = R.size
        phigrid = np.linspace(0, 2*np.pi, nphi+1)[:-1]
        R = np.meshgrid(R, phigrid, indexing="ij")[0]
        z = np.meshgrid(z, phigrid, indexing="ij")[0]
        t, phigrid = np.meshgrid(t, phigrid, indexing="ij")
        out = np.abs(self.eval_bfield(R.ravel(), phigrid.ravel(), z.ravel(),
                                      t.ravel(), evalb=True)["bphi"])
        out = out.reshape(nin, nphi)
        bmax = np.nanmax(out, axis=1)
        bmin = np.nanmin(out, axis=1)
        return (bmax - bmin) / (bmax + bmin)


    def _evaluateripplewell(self, R, z, t, nphi):
        """TBD
        """
        nin = R.size
        phigrid = np.linspace(0, 2*np.pi, nphi+1)[:-1]
        R = np.meshgrid(R, phigrid, indexing="ij")[0]
        z = np.meshgrid(z, phigrid, indexing="ij")[0]
        t, phigrid = np.meshgrid(t, phigrid, indexing="ij")
        out = self.eval_bfield(R.ravel(), phigrid.ravel(), z.ravel(),
                               t.ravel(), evalb=True)

        bnorm = np.sqrt( out["br"]*out["br"] + out["bphi"]*out["bphi"]
                         + out["bz"]*out["bz"] )
        bhat_r   = out["br"]   / bnorm
        bhat_phi = out["bphi"] / bnorm
        bhat_z   = out["bz"]   / bnorm

        bbar = (bhat_r * out["bphidr"] +  bhat_z * out["bphidz"]).reshape(nin, nphi)
        btil = (bhat_phi * out["bphidphi"] / R.ravel()).reshape(nin, nphi)

        return np.nanmean(np.abs(bbar), axis=1) \
            / np.nanmax(np.abs(btil), axis=1)

    def _plotripple(self, Rgrid, zgrid, time, nphi, axes=None):
        """TBD
        """
        R, z, t = np.meshgrid(Rgrid, zgrid, time, indexing="ij")
        rip = self.evaluateripple(R, z, t, nphi).reshape(Rgrid.size, zgrid.size)

        newfig = axes is None
        if newfig:
            plt.figure()
            axes = plt.gca()

        CS = axes.contour(Rgrid, zgrid, 100 * rip.transpose(),
                          [0.01, 0.05, 0.1, 0.5, 1, 5])
        axes.clabel(CS)

        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(Rgrid[0], Rgrid[-1])
        axes.set_ylim(zgrid[0], zgrid[-1])

        if newfig:
            plt.show(block=False)

    def _plotripplewell(self, Rgrid, zgrid, time, nphi, axes=None,
                       clevel=[-1,0,1], clabel=True, **kwargs):
        """TBD
        """
        R, z, t = np.meshgrid(Rgrid, zgrid, time, indexing="ij")
        rip = self.evaluateripplewell(R, z, t, nphi).reshape(Rgrid.size, zgrid.size)

        newfig = axes is None
        if newfig:
            plt.figure()
            axes = plt.gca()

        CS = axes.contour(Rgrid, zgrid, np.log10(rip.transpose()),
                          clevel, **kwargs)
        if clabel:
            axes.clabel(CS)

        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(Rgrid[0], Rgrid[-1])
        axes.set_ylim(zgrid[0], zgrid[-1])

        if newfig:
            plt.show(block=False)
