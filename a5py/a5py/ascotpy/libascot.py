"""
Python side of the interactive Python interface.

This module defines LibAscot class whose methods can be used in Python to call
Ascot5 functions (written in C) directly. The callable functions are defined in
library module libascot.c which must be compiled first with make libascot. This
module acts as a wrapper for those functions. More advanced functionality should
be implemented in other modules.

File: libascot.py
"""
import ctypes
import sys
import os
import warnings
import numpy as np

from scipy import interpolate

from a5py.ascotpy import ascotpy2

from ctypes.util import find_library
from numpy.ctypeslib import ndpointer
from a5py.ascot5io.ascot5 import Ascot

class AscotpyInitException(Exception):
    """Exception raised when Ascotpy object could not be initialized."""
    pass

class LibAscot:
    """
    An object representing a running ascot5 process.
    """

    DUMMY_QID = "".encode('UTF-8')

    def __init__(self, h5fn=None, libpath="libascot.so"):
        """
        Initialize and start Ascot5 process using given HDF5 file as an input.

        This function starts the process but does not read or initialize any
        data yet. All flags are set to False and function arguments are defined
        explicitly (because C uses strong typing whereas Python does not).

        Args:
            libpath : str, optional <br>
                Path to libascot.so library file. Default is to assume it
                is in the LD_LIBRARY_PATH.
            h5fn : str, optional <br>
                Path to HDF5 from which inputs are read. Default is "ascot.h5"
                in same folder the script is executed.
        """

        # Open library
        try:
            self.libascot = ctypes.CDLL(libpath)
        except Exception as e:
            msg = "\nCould not locate libascot.so. Try " \
                + "export LD_LIBRARY_PATH=/spam/ascot5 " \
                + "or "                                  \
                + "export LD_LIBRARY_PATH=/spam/ascot5:$LD_LIBRARY_PATH"
            raise AscotpyInitException(msg)

        # Initialize attributes
        self._offload_ready    = False
        self._nmrk             = ctypes.c_int32()
        self._diag_occupied    = False
        self.offload_data      = ascotpy2.struct_c__SA_offload_package()
        self.offload_array     = ctypes.POINTER(ctypes.c_double)()
        self.int_offload_array = ctypes.POINTER(ctypes.c_int   )()
        self.markers = ctypes.POINTER(ascotpy2.struct_c__SA_particle_state)()

        self.sim = ascotpy2.struct_c__SA_sim_offload_data()
        self.bfield_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self.efield_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self.plasma_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self.neutral_offload_array = ctypes.POINTER(ctypes.c_double)()
        self.wall_offload_array    = ctypes.POINTER(ctypes.c_double)()
        self.boozer_offload_array  = ctypes.POINTER(ctypes.c_double)()
        self.mhd_offload_array     = ctypes.POINTER(ctypes.c_double)()
        self.diag_offload_array    = ctypes.POINTER(ctypes.c_double)()

        self.wall_int_offload_array = ctypes.POINTER(ctypes.c_int)()

        # Check that the HDF5 file exists
        self.reload(h5fn)

        # MPI init needs to be called (only) once so that we can run simulations
        argc       = ctypes.c_int32()
        argc.value = 0
        argv       = ctypes.POINTER(ctypes.c_char)()
        self._mpi_rank = ctypes.c_int32()
        self._mpi_root = ctypes.c_int32()
        self._mpi_size = ctypes.c_int32()
        ascotpy2.mpi_interface_init(
            argc, ctypes.byref(argv), ctypes.byref(self.sim),
            ctypes.byref(self._mpi_rank), ctypes.byref(self._mpi_size),
            ctypes.byref(self._mpi_root))

        # Declare functions found in libascot. We have to do this manually
        # (instead of using clang2py) since we want to use numpy arrays
        # for arguments. This can be done using this pointer.
        real_p = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
        sim_p  = ctypes.POINTER(ascotpy2.struct_c__SA_sim_offload_data)
        cdbl_p = ctypes.POINTER(ctypes.c_double)

        # B field functions.
        fun = self.libascot.libascot_B_field_eval_B_dB
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p]

        fun = self.libascot.libascot_B_field_eval_rho
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p]

        fun = self.libascot.libascot_B_field_get_axis
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p]

        fun = self.libascot.libascot_B_field_eval_rhovals
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p,
                        ctypes.c_int,    ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        real_p, real_p, real_p]

        fun = self.libascot.libascot_E_field_eval_E
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p]

        fun = self.libascot.libascot_plasma_get_n_species
        fun.restype  = ctypes.c_int
        fun.argtypes = [sim_p, cdbl_p]

        fun = self.libascot.libascot_plasma_get_species_mass_and_charge
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p, real_p, real_p]

        fun = self.libascot.libascot_plasma_eval_background
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p]

        fun = self.libascot.libascot_neutral_eval_density
        fun.restype  = None
        fun.argtypes = [sim_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p]

        fun = self.libascot.libascot_boozer_eval_psithetazeta
        fun.restype  = ctypes.c_int
        fun.argtypes = [sim_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p]

        fun = self.libascot.libascot_boozer_eval_fun
        fun.restype  = ctypes.c_int
        fun.argtypes = [sim_p, cdbl_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p]

        fun = self.libascot.libascot_mhd_eval
        fun.restype  = ctypes.c_int
        fun.argtypes = [sim_p, cdbl_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p]

        fun = self.libascot.libascot_mhd_eval_perturbation
        fun.restype  = ctypes.c_int
        fun.argtypes = [sim_p, cdbl_p, cdbl_p, cdbl_p,
                        ctypes.c_int, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p]

        fun = self.libascot.libascot_eval_collcoefs
        fun.restype  = ctypes.c_int
        fun.argtypes = [sim_p, cdbl_p, cdbl_p, cdbl_p,
                        ctypes.c_int, real_p, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p]


    def reload(self, h5fn):
        """
        Change HDF5 file and free resources from old one.

        Args:
            h5fn : str <br>
                Name of the new HDF5 file.
        """
        self.free(bfield=True, efield=True, plasma=True, neutral=True,
                  boozer=True, mhd=True)

        if h5fn == None:
            self.h5fn = None
            self.hdf5 = Ascot()
        else:
            self.h5fn = h5fn if type(h5fn) is bytes else h5fn.encode('UTF-8')

            try:
                self.hdf5 = Ascot(str(self.h5fn,'UTF-8'))
            except:
                self.h5fn = None
                self.hdf5 = None
                raise

            self.sim.hdf5_in = self.h5fn


    def get_filepath(self):
        """
        Get path to the filename or None if no file is open.
        """
        return self.h5fn


    def init(self, bfield=False, efield=False, plasma=False, wall=False,
             neutral=False, boozer=False, mhd=False, ignorewarnings=False):
        if self._offload_ready:
            raise AscotpyInitException("This instance has been packed")

        a5 = self.hdf5
        inputs2read = ctypes.c_int32()
        def init_field(field, qid, array, currentqid, byterep,
                       intarray=None):
            """
            Initialize field's offload array which has given qid.
            """
            if isinstance(qid, str):
                # Convert given qid to bytes
                qid = qid.encode('UTF-8')
            elif qid == True:
                # qid == True means we use the active field
                qid = a5[field].active.get_qid().encode('UTF-8')
            else:
                # This field is not read
                return currentqid

            if qid == currentqid:
                if not ignorewarnings:
                    warnings.warn(field + " already initialized.", Warning)
                return currentqid
            else:
                if qid == LibAscot.DUMMY_QID:
                    ascotpy2.libascot_deallocate(array)
                    if intarray is not None:
                        ascotpy2.libascot_deallocate(intarray)
                inputs2read.value = inputs2read.value | byterep
                return qid

        self.sim.qid_bfield = init_field(
            "bfield", bfield, self.bfield_offload_array, self.sim.qid_bfield,
            ascotpy2.hdf5_input_bfield)

        self.sim.qid_efield = init_field(
            "efield", efield, self.efield_offload_array, self.sim.qid_efield,
            ascotpy2.hdf5_input_efield)

        self.sim.qid_plasma = init_field(
            "plasma", plasma, self.plasma_offload_array, self.sim.qid_plasma,
            ascotpy2.hdf5_input_plasma)

        self.sim.qid_neutral = init_field(
            "neutral", neutral, self.neutral_offload_array, self.sim.qid_neutral,
            ascotpy2.hdf5_input_neutral)

        self.sim.qid_wall = init_field(
            "wall", wall, self.wall_offload_array, self.sim.qid_wall,
            ascotpy2.hdf5_input_wall, self.wall_int_offload_array)

        self.sim.qid_boozer = init_field(
            "boozer", boozer, self.boozer_offload_array, self.sim.qid_boozer,
            ascotpy2.hdf5_input_boozer)

        self.sim.qid_mhd = init_field(
            "mhd", mhd, self.mhd_offload_array, self.sim.qid_mhd,
            ascotpy2.hdf5_input_mhd)

        ascotpy2.hdf5_interface_read_input(
            ctypes.byref(self.sim),
            inputs2read,
            ctypes.byref(self.bfield_offload_array),
            ctypes.byref(self.efield_offload_array),
            ctypes.byref(self.plasma_offload_array),
            ctypes.byref(self.neutral_offload_array),
            ctypes.byref(self.wall_offload_array),
            ctypes.byref(self.wall_int_offload_array),
            ctypes.byref(self.boozer_offload_array),
            ctypes.byref(self.mhd_offload_array),
            None, # Marker array
            None  # Number of markers that were read
            )


    def free(self, bfield=False, efield=False, plasma=False, wall=False,
             neutral=False, boozer=False, mhd=False):
        if self._offload_ready:
            raise AscotpyInitException("This instance has been packed")

        if bfield:
            self.sim.qid_bfield = LibAscot.DUMMY_QID
            self.sim.B_offload_data.offload_array_length = 0
            ascotpy2.libascot_deallocate(self.bfield_offload_array)
        if efield:
            self.sim.qid_efield = LibAscot.DUMMY_QID
            self.sim.E_offload_data.offload_array_length = 0
            ascotpy2.libascot_deallocate(self.efield_offload_array)
        if plasma:
            self.sim.qid_plasma = LibAscot.DUMMY_QID
            self.sim.plasma_offload_data.offload_array_length = 0
            ascotpy2.libascot_deallocate(self.plasma_offload_array)
        if wall:
            self.sim.qid_wall = LibAscot.DUMMY_QID
            self.sim.wall_offload_data.offload_array_length = 0
            self.sim.wall_offload_data.int_offload_array_length = 0
            ascotpy2.libascot_deallocate(self.wall_offload_array)
            ascotpy2.libascot_deallocate(self.wall_int_offload_array)
        if neutral:
            self.sim.qid_neutral = LibAscot.DUMMY_QID
            self.sim.neutral_offload_data.offload_array_length = 0
            ascotpy2.libascot_deallocate(self.neutral_offload_array)
        if boozer:
            self.sim.qid_boozer = LibAscot.DUMMY_QID
            self.sim.boozer_offload_data.offload_array_length = 0
            ascotpy2.libascot_deallocate(self.boozer_offload_array)
        if mhd:
            self.sim.qid_mhd = LibAscot.DUMMY_QID
            self.sim.mhd_offload_data.offload_array_length = 0
            ascotpy2.libascot_deallocate(self.mhd_offload_array)


    def init_from_run(self, run, bfield=False, efield=False, plasma=False,
                      wall=False, neutral=False, boozer=False, mhd=False):
        """
        Initialize data that was used in the given run.
        """
        if isinstance(run, str):
            run = self.hdf5["q"+run]

        if bfield:
            bfield = run.bfield.get_qid()
        if efield:
            efield = run.efield.get_qid()
        if plasma:
            plasma = run.plasma.get_qid()
        if wall:
            wall = run.wall.get_qid()
        if neutral:
            neutral = run.neutral.get_qid()
        if boozer:
            boozer = run.boozer.get_qid()
        if mhd:
            mhd = run.mhd.get_qid()

        self.init(bfield=bfield, efield=efield, plasma=plasma, wall=wall,
                  neutral=neutral, boozer=boozer, mhd=mhd, ignorewarnings=True)


    def pack(self):
        """
        Pack offload arrays as one, making this instance ready for simulation.

        Note that inputs cannot be changed or freed before calling unpack. Make
        sure all required data is initialized before packing.
        """
        if self._offload_ready:
            raise AscotpyInitException("This instance is already packed")

        # This call internally frees individual offload arrays and initializes
        # the common ones.
        ascotpy2.pack_offload_array(
            ctypes.byref(self.sim), ctypes.byref(self.offload_data),
            self.bfield_offload_array, self.efield_offload_array,
            self.plasma_offload_array, self.neutral_offload_array,
            self.wall_offload_array, self.wall_int_offload_array,
            self.boozer_offload_array, self.mhd_offload_array,
            ctypes.byref(self.offload_array),
            ctypes.byref(self.int_offload_array))
        self._offload_ready = True

        # Set pointers to correct locations (based on the order arrays are
        # packed) so that we can continue using evaluation routines.
        def advance(prevptr, increment):
            ptr = ctypes.addressof(prevptr) + 2 * increment
            ptr = ctypes.cast(ctypes.c_void_p(ptr),
                              ctypes.POINTER(ctypes.c_double))
            return ptr

        self.bfield_offload_array = self.offload_array
        self.efield_offload_array = advance(
            self.bfield_offload_array.contents,
            self.sim.B_offload_data.offload_array_length)
        self.plasma_offload_array = advance(
            self.efield_offload_array.contents,
            self.sim.E_offload_data.offload_array_length)
        self.neutral_offload_array = advance(
            self.plasma_offload_array.contents,
            self.sim.plasma_offload_data.offload_array_length)
        self.wall_offload_array = advance(
            self.neutral_offload_array.contents,
            self.sim.neutral_offload_data.offload_array_length)
        self.boozer_offload_array = advance(
            self.wall_offload_array.contents,
            self.sim.wall_offload_data.offload_array_length)
        self.mhd_offload_array = advance(
            self.boozer_offload_array.contents,
            self.sim.boozer_offload_data.offload_array_length)

        self.wall_int_offload_array = self.int_offload_array


    def unpack(self, bfield=True, efield=True, plasma=True, wall=True,
               neutral=True, boozer=True, mhd=True):
        """
        Unpack simulation arrays, i.e. free offload array and re-read data.

        After unpacking the inputs can be changed or freed again but simulations
        cannot be performed.
        """
        if not self._offload_ready:
            raise AscotpyInitException("This instance hasn't been packed")

        ascotpy2.libascot_deallocate(self.offload_array)
        ascotpy2.libascot_deallocate(self.int_offload_array)
        self.offload_data.offload_array_length     = 0
        self.offload_data.int_offload_array_length = 0
        self.offload_data.unpack_pos               = 0
        self.offload_data.int_unpack_pos           = 0

        inputs2read = {}
        def readornot(name, read, qid):
            if read:
                inputs2read[name] = qid.decode("utf-8")
            return LibAscot.DUMMY_QID

        self.sim.qid_bfield  = readornot("bfield",  bfield,
                                         self.sim.qid_bfield)
        self.sim.qid_efield  = readornot("efield",  efield,
                                         self.sim.qid_efield)
        self.sim.qid_plasma  = readornot("plasma",  plasma,
                                         self.sim.qid_plasma)
        self.sim.qid_neutral = readornot("neutral", neutral,
                                         self.sim.qid_neutral)
        self.sim.qid_wall    = readornot("wall",    wall,
                                         self.sim.qid_wall)
        self.sim.qid_boozer  = readornot("boozer",  boozer,
                                         self.sim.qid_boozer)
        self.sim.qid_mhd     = readornot("mhd",     mhd,
                                         self.sim.qid_mhd)

        self._offload_ready = False
        self.init(**inputs2read)


    def eval_bfield(self, R, phi, z, t, evalb=False, evalrho=False,
                    evalaxis=False):
        """
        Evaluate magnetic field quantities at given coordinates.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].
            evalb : bool, optional <br>
                Evaluate magnetic field vector and derivatives.
            evalrho : bool, optional <br>
                Evaluate poloidal flux and normalized poloidal flux.
            evalaxis : bool, optional <br>
                Evaluate magnetic axis.

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
            "Magnetic field not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}

        if evalb:
            out["br"]       = np.zeros(R.shape, dtype="f8") + np.nan
            out["bphi"]     = np.zeros(R.shape, dtype="f8") + np.nan
            out["bz"]       = np.zeros(R.shape, dtype="f8") + np.nan
            out["brdr"]     = np.zeros(R.shape, dtype="f8") + np.nan
            out["brdphi"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["brdz"]     = np.zeros(R.shape, dtype="f8") + np.nan
            out["bphidr"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["bphidphi"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["bphidz"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["bzdr"]     = np.zeros(R.shape, dtype="f8") + np.nan
            out["bzdphi"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["bzdz"]     = np.zeros(R.shape, dtype="f8") + np.nan

            self.libascot.libascot_B_field_eval_B_dB(
                ctypes.byref(self.sim), self.bfield_offload_array,
                Neval, R, phi, z, t, out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"], out["bphidr"],
                out["bphidphi"], out["bphidz"], out["bzdr"], out["bzdphi"],
                out["bzdz"])

        if evalrho:
            out["rho"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["psi"] = np.zeros(R.shape, dtype="f8") + np.nan
            self.libascot.libascot_B_field_eval_rho(
                ctypes.byref(self.sim), self.bfield_offload_array,
                Neval, R, phi, z, t, out["rho"], out["psi"])

        if evalaxis:
            out["axisr"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["axisz"] = np.zeros(R.shape, dtype="f8") + np.nan
            self.libascot.libascot_B_field_get_axis(
                ctypes.byref(self.sim), self.bfield_offload_array,
                Neval, phi, out["axisr"], out["axisz"])

        return out


    def eval_efield(self, R, phi, z, t):
        """
        Evaluate electric field quantities at given coordinates.

        Note that magnetic field has to be initialized as well.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
            "Magnetic field not initialized"
        assert self.sim.qid_efield != LibAscot.DUMMY_QID, \
            "Electric field not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["er"]   = np.zeros(R.shape, dtype="f8") + np.nan
        out["ephi"] = np.zeros(R.shape, dtype="f8") + np.nan
        out["ez"]   = np.zeros(R.shape, dtype="f8") + np.nan

        self.libascot.libascot_E_field_eval_E(
            ctypes.byref(self.sim), self.efield_offload_array,
            Neval, R, phi, z, t, out["er"], out["ephi"], out["ez"])

        return out


    def get_plasmaspecies(self):
        """
        Get plasma species information.

        Returns;
            Dictionary containing nspecies, and anum, znum, charge, and mass for
            each species.
        """
        assert self.sim.qid_plasma != LibAscot.DUMMY_QID, \
            "Plasma not initialized"

        out = {}
        out["nspecies"] = \
            self.libascot.libascot_plasma_get_n_species(
                ctypes.byref(self.sim), self.plasma_offload_array)

        out["mass"]     = np.zeros((out["nspecies"],), dtype="f8")
        out["charge"]   = np.zeros((out["nspecies"],), dtype="f8")
        self.libascot.libascot_plasma_get_species_mass_and_charge(
            ctypes.byref(self.sim), self.plasma_offload_array,
            out["mass"], out["charge"])

        return out


    def eval_plasma(self, R, phi, z, t):
        """
        Evaluate plasma quantities at given coordinates.

        Note that magnetic field has to be initialized as well.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
            "Magnetic field not initialized"
        assert self.sim.qid_plasma != LibAscot.DUMMY_QID, \
            "Plasma not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size

        # Allocate enough space for electrons and all ion species.
        nspecies = self.libascot.libascot_plasma_get_n_species()
        rawdens = np.zeros((Neval*nspecies,), dtype="f8") + np.nan
        rawtemp = np.zeros((Neval*nspecies,), dtype="f8") + np.nan

        self.libascot.libascot_plasma_eval_background(
            ctypes.byref(self.sim), self.bfield_offload_array,
            self.plasma_offload_array,
            Neval, R, phi, z, t, rawdens, rawtemp)

        out = {}
        out["ne"] = rawdens[0:Neval]
        out["te"] = rawtemp[0:Neval]
        for i in range(1, nspecies):
            out["ni"+str(i)] = rawdens[(Neval)*i:(Neval)*(i+1)]
            out["ti"+str(i)] = rawtemp[(Neval)*i:(Neval)*(i+1)]

        return out


    def eval_neutral(self, R, phi, z, t):
        """
        Evaluate neutral quantities at given coordinates.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.sim.qid_neutral != LibAscot.DUMMY_QID, \
            "Neutral data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["density"] = np.zeros(R.shape, dtype="f8") + np.nan

        self.libascot.libascot_neutral_eval_density(
            ctypes.byref(self.sim), self.neutral_offload_array,
            Neval, R, phi, z, t, out["density"])

        return out


    def eval_boozer(self, R, phi, z, t, evalfun=False):
        """
        Evaluate boozer quantities at given coordinates.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].
            evalfun : boolean, optional <br>
                evaluate the boozer functions instead of the coordinates.

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        if evalfun:
            assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
                "Magnetic field not initialized"
        assert self.sim.qid_boozer != LibAscot.DUMMY_QID, \
            "Boozer data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}

        if evalfun:
            out["qprof"]      = np.zeros(R.shape, dtype="f8") + np.nan
            out["jacobian"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["jacobianb2"] = np.zeros(R.shape, dtype="f8") + np.nan
            self.libascot.libascot_boozer_eval_fun(
                ctypes.byref(self.sim), self.bfield_offload_array,
                self.boozer_offload_array,
                Neval, R, phi, z, t, out["qprof"], out["jacobian"],
                out["jacobianb2"])
        else:
            out["psi"]        = np.zeros(R.shape, dtype="f8") + np.nan
            out["theta"]      = np.zeros(R.shape, dtype="f8") + np.nan
            out["zeta"]       = np.zeros(R.shape, dtype="f8") + np.nan
            out["dpsidr"]     = np.zeros(R.shape, dtype="f8") + np.nan
            out["dpsidphi"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dpsidz"]     = np.zeros(R.shape, dtype="f8") + np.nan
            out["dthetadr"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dthetadphi"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["dthetadz"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dzetadr"]    = np.zeros(R.shape, dtype="f8") + np.nan
            out["dzetadphi"]  = np.zeros(R.shape, dtype="f8") + np.nan
            out["dzetadz"]    = np.zeros(R.shape, dtype="f8") + np.nan
            out["rho"]        = np.zeros(R.shape, dtype="f8") + np.nan
            self.libascot.libascot_boozer_eval_psithetazeta(
                ctypes.byref(self.sim), self.boozer_offload_array,
                Neval, R, phi, z, t, out["psi"], out["theta"], out["zeta"],
                out["dpsidr"], out["dpsidphi"], out["dpsidz"],
                out["dthetadr"], out["dthetadphi"], out["dthetadz"],
                out["dzetadr"], out["dzetadphi"], out["dzetadz"], out["rho"])

        return out


    def eval_mhd(self, R, phi, z, t, evalpot=False):
        """
        Evaluate mhd perturbation EM components at given coordinates.

        Args:
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                time coordinates where data is evaluated [s].
            evalpot : boolean, optional <br>
                flag to evaluate also alpha and phi.

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
            "Magnetic field not initialized"
        assert self.sim.qid_boozer != LibAscot.DUMMY_QID, \
            "Boozer data not initialized"
        assert self.sim.qid_mhd != LibAscot.DUMMY_QID, \
            "MHD data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["mhd_br"]   = np.zeros(R.shape, dtype="f8") + np.nan
        out["mhd_bphi"] = np.zeros(R.shape, dtype="f8") + np.nan
        out["mhd_bz"]   = np.zeros(R.shape, dtype="f8") + np.nan
        out["mhd_er"]   = np.zeros(R.shape, dtype="f8") + np.nan
        out["mhd_ephi"] = np.zeros(R.shape, dtype="f8") + np.nan
        out["mhd_ez"]   = np.zeros(R.shape, dtype="f8") + np.nan
        out["mhd_phi"]  = np.zeros(R.shape, dtype="f8") + np.nan

        self.libascot.libascot_mhd_eval_perturbation(
            ctypes.byref(self.sim), self.bfield_offload_array,
            self.boozer_offload_array, self.mhd_offload_array,
            Neval, R, phi, z, t, out["mhd_br"], out["mhd_bphi"], out["mhd_bz"],
            out["mhd_er"], out["mhd_ephi"], out["mhd_ez"], out["mhd_phi"])

        if evalpot:
            out["alpha"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["phi"]   = np.zeros(R.shape, dtype="f8") + np.nan

            out["dadt"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dadr"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dadphi"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["dadz"]   = np.zeros(R.shape, dtype="f8") + np.nan

            out["dphidt"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dphidr"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["dphidphi"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["dphidz"]   = np.zeros(R.shape, dtype="f8") + np.nan

            self.libascot.libascot_mhd_eval(
                Neval, R, phi, z, t, out["alpha"], out["dadr"], out["dadphi"],
                out["dadz"], out["dadt"], out["phi"], out["dphidr"],
                out["dphidphi"], out["dphidz"], out["dphidt"])

        return out


    def eval_collcoefs(self, ma, qa, R, phi, z, t, va):
        """
        Evaluate Coulomb collision coefficients for given particle and coords.

        Collision coefficients are evaluated by interpolating the plasma
        parameters on given coordinates, and using those and given test particle
        mass and charge to evaluate the coefficients as a function of test
        particle velocity.

        Args:
            ma : float <br>
                Test particle mass [kg].
            qa : float <br>
                Test particle charge [C].
            R : array_like <br>
                R coordinates where data is evaluated [m].
            phi : array_like <br>
                phi coordinates where data is evaluated [rad].
            z : array_like <br>
                z coordinates where data is evaluated [m].
            t : array_like <br>
                Time coordinates where data is evaluated [s].
            va : array_like <br>
                Test particle velocities [m/s].

        Returns:
            Dictionary with collision coefficients of shape
            (R.size, n_species, va.size).
        """
        assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
            "Magnetic field not initialized"
        assert self.sim.qid_plasma != LibAscot.DUMMY_QID, \
            "Plasma not initialized"


        ma  = float(ma)
        qa  = float(qa)
        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")
        va  = np.asarray(va).ravel().astype(dtype="f8")
        Neval = va.size

        n_species = \
            self.libascot.libascot_plasma_get_n_species(
                ctypes.byref(self.sim), self.plasma_offload_array)

        out = {}
        out["F"]      = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["Dpara"]  = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["Dperp"]  = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["K"]      = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["nu"]     = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["Q"]      = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["dQ"]     = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["dDpara"] = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["clog"]   = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["mu0"]    = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["mu1"]    = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["dmu0"]   = np.zeros((R.size, n_species, va.size), dtype="f8")

        for i in range(R.size):
            F      = np.zeros((n_species, va.size), dtype="f8")
            Dpara  = np.zeros((n_species, va.size), dtype="f8")
            Dperp  = np.zeros((n_species, va.size), dtype="f8")
            K      = np.zeros((n_species, va.size), dtype="f8")
            nu     = np.zeros((n_species, va.size), dtype="f8")
            Q      = np.zeros((n_species, va.size), dtype="f8")
            dQ     = np.zeros((n_species, va.size), dtype="f8")
            dDpara = np.zeros((n_species, va.size), dtype="f8")
            clog   = np.zeros((n_species, va.size), dtype="f8")
            mu0    = np.zeros((n_species, va.size), dtype="f8")
            mu1    = np.zeros((n_species, va.size), dtype="f8")
            dmu0   = np.zeros((n_species, va.size), dtype="f8")
            self.libascot.libascot_eval_collcoefs(
                ctypes.byref(self.sim), self.bfield_offload_array,
                self.plasma_offload_array,
                Neval, va, R[i], phi[i], z[i], t[i], ma, qa, F, Dpara, Dperp,
                K, nu, Q, dQ, dDpara, clog, mu0, mu1, dmu0)
            out["F"][i,:,:]      = F[:,:]
            out["Dpara"][i,:,:]  = Dpara[:,:]
            out["Dperp"][i,:,:]  = Dperp[:,:]
            out["K"][i,:,:]      = K[:,:]
            out["nu"][i,:,:]     = nu[:,:]
            out["Q"][i,:,:]      = Q[:,:]
            out["dQ"][i,:,:]     = dQ[:,:]
            out["dDpara"][i,:,:] = dDpara[:,:]
            out["clog"][i,:,:]   = clog[:,:]
            out["mu0"][i,:,:]    = mu0[:,:]
            out["mu1"][i,:,:]    = mu1[:,:]
            out["dmu0"][i,:,:]   = dmu0[:,:]

        return out


    def get_rhotheta_rz(self, rhovals, theta, phi, time):
        assert self.sim.qid_bfield != LibAscot.DUMMY_QID, \
            "Magnetic field not initialized"

        rhovals = np.asarray(rhovals).ravel().astype(dtype="f8")

        ngrid = 100
        r   = np.zeros((ngrid,), dtype="f8")
        z   = np.zeros((ngrid,), dtype="f8")
        rho = np.zeros((ngrid,), dtype="f8")

        self.libascot.libascot_B_field_eval_rhovals(
            ctypes.byref(self.sim), self.bfield_offload_array,
            ngrid, np.min(rhovals), np.max(rhovals), theta, phi, time,
            r, z, rho)

        rho[0]  = rhovals[0]
        rho[-1] = rhovals[-1]

        fr = interpolate.interp1d(rho, r, fill_value="extrapolate")
        fz = interpolate.interp1d(rho, z, fill_value="extrapolate")

        return (fr(rhovals), fz(rhovals))


def test():
    """
    For testing purposes.
    """
    ascot = LibAscot(h5fn="ascot.h5")
    ascot.init(bfield=True, efield=True, plasma=True, wall=True,
               neutral=True)

    R   = np.array([6.2,   7, 8])
    phi = np.array([  0,   0, 0])
    z   = np.array([0.0, 0.2, 0.2])
    t   = np.array([0.0])

    bvals       = ascot.eval_bfield(R, phi, z, t, evalb=True, evalpsi=True,
                                            evalrho=True, evalaxis=True)
    evals       = ascot.eval_efield(R, phi, z, t)
    plasmavals  = ascot.eval_plasma(R, phi, z, t)
    neutralvals = ascot.eval_neutral(R, phi, z, t)

    print(bvals)

    ascot.free(bfield=True, efield=True, plasma=True, wall=True,
               neutral=True)


if __name__ == '__main__':
    test()
