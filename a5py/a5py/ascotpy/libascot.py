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
import unyt

from scipy import interpolate

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
        self.bfield_initialized  = None
        self.efield_initialized  = None
        self.plasma_initialized  = None
        self.neutral_initialized = None
        self.wall_initialized    = None
        self.boozer_initialized  = None
        self.mhd_initialized     = None
        self.asigma_initialized  = None

        # Declare functions found in libascot

        real_p = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")

        # Init and free functions.
        try:
            fun = self.libascot.libascot_init
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                            ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                            ctypes.c_char_p, ctypes.c_char_p]
        except AttributeError:
            warnings.warn("libascot_init not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_free
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            ctypes.c_int]
        except AttributeError:
            warnings.warn("libascot_free not found", Warning)
            pass

        # B field functions.
        try:
            fun = self.libascot.libascot_B_field_eval_B_dB
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_B_field_eval_B_dB not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_B_field_eval_rho
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_B_field_eval_rho not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_B_field_get_axis
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_B_field_get_axis not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_B_field_eval_rhovals
            fun.restype  = None
            fun.argtypes = [ctypes.c_int,    ctypes.c_double, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_B_field_eval_rhovals not found", Warning)
            pass

        # E field functions.
        try:
            fun = self.libascot.libascot_E_field_eval_E
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_E_field_eval_E not found", Warning)
            pass

        # Plasma functions.
        try:
            fun = self.libascot.libascot_plasma_get_n_species
            fun.restype  = ctypes.c_int
        except AttributeError:
            warnings.warn("libascot_plasma_get_n_species not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_plasma_get_species_mass_and_charge
            fun.restype  = None
            fun.argtypes = [real_p, real_p]
        except AttributeError:
            warnings.warn(
                "libascot_plasma_get_species_mass_and_charge not found",
                Warning)
            pass

        try:
            fun = self.libascot.libascot_plasma_eval_background
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_plasma_eval_background not found", Warning)
            pass

        # Neutral functions.
        try:
            fun = self.libascot.libascot_neutral_eval_density
            fun.restype  = None
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p]
        except AttributeError:
            warnings.warn("libascot_neutral_eval_density not found", Warning)
            pass

        # Boozer and MHD functions.
        try:
            fun = self.libascot.libascot_boozer_eval_psithetazeta
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p,
                            real_p]
        except AttributeError:
            warnings.warn("libascot_boozer_eval_psithetazeta not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_boozer_eval_fun
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_boozer_eval_fun not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_mhd_eval
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_mhd_eval not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_mhd_eval_perturbation
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p,
                            real_p]
        except AttributeError:
            warnings.warn("libascot_mhd_eval_perturbation not found", Warning)
            pass

        # Collision coefficients.
        try:
            fun = self.libascot.libascot_eval_collcoefs
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_int, real_p, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_eval_collcoefs not found", Warning)
            pass

        # Check that the HDF5 file exists
        self.reload(h5fn)


    def reload(self, h5fn):
        """
        Change HDF5 file and free resources from old one.

        Args:
            h5fn : str <br>
                Name of the new HDF5 file.
        """
        self.free(bfield  = self.bfield_initialized,
                  efield  = self.efield_initialized,
                  plasma  = self.plasma_initialized,
                  wall    = self.wall_initialized,
                  neutral = self.neutral_initialized,
                  mhd     = self.mhd_initialized,
                  boozer  = self.boozer_initialized,
                  asigma=self.asigma_initialized)

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


    def get_filepath(self):
        """
        Get path to the filename or None if no file is open.
        """
        return self.h5fn


    def init(self, bfield=False, efield=False, plasma=False, wall=False,
             neutral=False, boozer=False, mhd=False, asigma=False,
             ignorewarnings=False):
        """
        Initialize input data.

        If argument is False, that data is not read. If True, then active input
        is read. Optionally, a QID of the input to be read can be given.

        Args:
            bfield : bool or str, optional <br>
                Flag for initializing magnetic field.
            efield : bool or str, optional <br>
                Flag for initializing electric field.
            plasma : bool or str optional <br>
                Flag for initializing plasma data.
            wall : bool or str optional <br>
                Flag for initializing wall data.
            neutral : bool or str optional <br>
                Flag for initializing neutral data.
            boozer : bool, optional <br>
                Flag for initializing boozer data.
            mhd : bool, optional <br>
                Flag for initializing mhd data.
            asigma : bool or str optional <br>
                Flag for initializing atomic sigma data.
            ignorewarnings : bool, optional <br>
                If True, warnings are not displayed when trying to initialize
                data that is already initialized

        Raises:
            RuntimeError if initialization failed.
        """
        a5 = self.hdf5

        if isinstance(bfield, str):
            qid    = bfield.encode('UTF-8')
            bfield = True
        elif bfield:
            qid = a5.bfield.active.get_qid().encode('UTF-8')

        if bfield and self.bfield_initialized == qid:
            if not ignorewarnings:
                warnings.warn("Magnetic field already initialized.", Warning)
        elif bfield:
            self.free(bfield=True)
            if self.libascot.libascot_init(self.h5fn, qid, None, None, None,
                                           None, None, None, None):
                raise RuntimeError("Failed to initialize magnetic field")

            self.bfield_initialized = qid


        if isinstance(efield, str):
            qid    = efield.encode('UTF-8')
            efield = True
        elif efield:
            qid = a5.efield.active.get_qid().encode('UTF-8')

        if efield and self.efield_initialized == qid:
            if not ignorewarnings:
                warnings.warn("Electric field already initialized.", Warning)
        elif efield:
            self.free(efield=True)
            if self.libascot.libascot_init(self.h5fn, None, qid, None, None,
                                           None, None, None, None) :
                raise RuntimeError("Failed to initialize electric field")

            self.efield_initialized = qid


        if isinstance(plasma, str):
            qid    = plasma.encode('UTF-8')
            plasma = True
        elif plasma:
            qid = a5.plasma.active.get_qid().encode('UTF-8')

        if plasma and self.plasma_initialized == qid:
            if not ignorewarnings:
                warnings.warn("Plasma already initialized.", Warning)
        elif plasma:
            self.free(plasma=True)
            if self.libascot.libascot_init(self.h5fn, None, None, qid, None,
                                           None, None, None, None) :
                raise RuntimeError("Failed to initialize plasma")

            self.plasma_initialized = qid


        if isinstance(wall, str):
            qid  = wall.encode('UTF-8')
            wall = True
        elif wall:
            qid = a5.wall.active.get_qid().encode('UTF-8')

        if wall and self.wall_initialized == qid:
            if not ignorewarnings:
                warnings.warn("Wall already initialized.", Warning)
        elif wall:
            self.free(wall=True)
            if self.libascot.libascot_init(self.h5fn, None, None, None, qid,
                                           None, None, None, None) :
                raise RuntimeError("Failed to initialize wall")

            self.wall_initialized = qid


        if isinstance(neutral, str):
            qid     = neutral.encode('UTF-8')
            neutral = True
        elif neutral:
            qid = a5.neutral.active.get_qid().encode('UTF-8')

        if neutral and self.neutral_initialized == qid:
            if not ignorewarnings:
                warnings.warn("Neutral data already initialized.", Warning)
        elif neutral:
            self.free(neutral=True)
            if self.libascot.libascot_init(self.h5fn, None, None, None, None,
                                           qid, None, None, None) :
                raise RuntimeError("Failed to initialize neutral data")

            self.neutral_initialized = qid

        if isinstance(boozer, str):
            qid    = boozer.encode('UTF-8')
            boozer = True
        elif boozer:
            qid = a5.boozer.active.get_qid().encode('UTF-8')

        if boozer and self.boozer_initialized == qid:
            if not ignorewarnings:
                warnings.warn("Boozer data already initialized.", Warning)
        elif boozer:
            self.free(boozer=True)
            if self.libascot.libascot_init(self.h5fn, None, None, None, None,
                                           None, qid, None) :
                raise RuntimeError("Failed to initialize boozer data")

            self.boozer_initialized = qid

        if isinstance(mhd, str):
            qid = mhd.encode('UTF-8')
            mhd = True
        elif mhd:
            qid = a5.mhd.active.get_qid().encode('UTF-8')

        if mhd and self.mhd_initialized == qid:
            if not ignorewarnings:
                warnings.warn("MHD data already initialized.", Warning)
        elif mhd:
            self.free(mhd=True)
            if self.libascot.libascot_init(self.h5fn, None, None, None, None,
                                           None, None, qid) :
                raise RuntimeError("Failed to initialize MHD data")

            self.mhd_initialized = qid


        if isinstance(asigma, str):
            qid    = asigma.encode('UTF-8')
            asigma = True
        elif asigma:
            qid = a5.asigma.active.get_qid().encode('UTF-8')

        if asigma and self.asigma_initialized:
            warnings.warn("Atomic sigma data already initialized.", Warning)
        elif asigma:
            if self.libascot.libascot_init(self.h5fn, None, None, None, None,
                                           None, qid) :
                raise RuntimeError("Failed to initialize atomic sigma data")

            self.asigma_initialized = True


    def free(self, bfield=False, efield=False, plasma=False, wall=False,
             neutral=False, boozer=False, mhd=False, asigma=False):
        """
        Free input data.

        Args:
            bfield : bool, optional <br>
                Flag for freeing magnetic field.
            efield : bool, optional <br>
                Flag for freeing electric field.
            plasma : bool, optional <br>
                Flag for freeing plasma data.
            wall : bool, optional <br>
                Flag for freeing wall data.
            neutral : bool, optional <br>
                Flag for freeing neutral data.
            boozer : bool, optional <br>
                Flag for freeing boozer data.
            mhd : bool, optional <br>
                Flag for freeing mhd data.
            asigma : bool, optional <br>
                Flag for freeing atomic sigma data.
        """
        if bfield and self.bfield_initialized is not None:
            self.libascot.libascot_free(1, 0, 0, 0, 0, 0, 0)
            self.bfield_initialized = None

        if efield and self.efield_initialized is not None:
            self.libascot.libascot_free(0, 1, 0, 0, 0, 0, 0)
            self.efield_initialized = None

        if plasma and self.plasma_initialized is not None:
            self.libascot.libascot_free(0, 0, 1, 0, 0, 0, 0)
            self.plasma_initialized = None

        if wall and self.wall_initialized is not None:
            self.libascot.libascot_free(0, 0, 0, 1, 0, 0, 0)
            self.wall_initialized = None

        if neutral and self.neutral_initialized is not None:
            self.libascot.libascot_free(0, 0, 0, 0, 1, 0, 0)
            self.neutral_initialized = None

        if boozer and self.boozer_initialized is not None:
            self.libascot.libascot_free(0, 0, 0, 0, 0, 1, 0)
            self.boozer_initialized = None

        if mhd and self.mhd_initialized is not None:
            self.libascot.libascot_free(0, 0, 0, 0, 0, 0, 1)
            self.mhd_initialized = None

        if asigma and self.asigma_initialized:
            self.libascot.libascot_free(0, 0, 0, 0, 0, 0, 0, 1)
            self.asigma_initialized = None


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
        assert self.bfield_initialized is not None, \
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
                Neval, R, phi, z, t,
                out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"],
                out["bphidr"], out["bphidphi"], out["bphidz"],
                out["bzdr"], out["bzdphi"], out["bzdz"])

        if evalrho:
            out["rho"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["psi"] = np.zeros(R.shape, dtype="f8") + np.nan
            self.libascot.libascot_B_field_eval_rho(Neval, R, phi, z, t,
                                                    out["rho"], out["psi"])

        if evalaxis:
            out["axisr"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["axisz"] = np.zeros(R.shape, dtype="f8") + np.nan
            self.libascot.libascot_B_field_get_axis(Neval, phi, out["axisr"],
                                                    out["axisz"])

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
        assert self.bfield_initialized is not None, \
            "Magnetic field not initialized"
        assert self.efield_initialized is not None, \
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
            Neval, R, phi, z, t, out["er"], out["ephi"], out["ez"])

        return out


    def get_plasmaspecies(self):
        """
        Get plasma species information.

        Returns;
            Dictionary containing nspecies, and anum, znum, charge, and mass for
            each species.
        """
        assert self.plasma_initialized is not None, \
            "Plasma not initialized"

        out = {}
        out["nspecies"] = self.libascot.libascot_plasma_get_n_species()
        out["mass"]     = np.zeros((out["nspecies"],), dtype="f8")
        out["charge"]   = np.zeros((out["nspecies"],), dtype="f8")
        self.libascot.libascot_plasma_get_species_mass_and_charge(
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
        assert self.bfield_initialized is not None, \
            "Magnetic field not initialized"
        assert self.plasma_initialized is not None, \
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
        assert self.neutral_initialized is not None, \
            "Neutral data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["density"] = np.zeros(R.shape, dtype="f8") + np.nan

        self.libascot.libascot_neutral_eval_density(Neval, R, phi, z, t,
                                                    out["density"])

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
        assert self.boozer_initialized is not None, \
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
        assert self.bfield_initialized is not None, \
            "Magnetic field not initialized"
        assert self.boozer_initialized is not None, \
            "Boozer data not initialized"
        assert self.mhd_initialized is not None, \
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
        assert self.bfield_initialized is not None, \
            "Magnetic field not initialized"
        assert self.plasma_initialized is not None, \
            "Plasma not initialized"


        ma  = float(ma)
        qa  = float(qa)
        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")
        va  = np.asarray(va).ravel().astype(dtype="f8")
        Neval = va.size

        n_species = self.libascot.libascot_plasma_get_n_species()

        out = {}
        temp = np.zeros((R.size, n_species, va.size), dtype="f8")
        out["F"]      = np.copy(temp) * unyt.m / unyt.s**2
        out["Dpara"]  = np.copy(temp) * unyt.m**2 / unyt.s**3
        out["Dperp"]  = np.copy(temp) * unyt.m**2 / unyt.s**3
        out["K"]      = np.copy(temp) * unyt.m / unyt.s**2
        out["nu"]     = np.copy(temp) / unyt.s
        out["Q"]      = np.copy(temp) * unyt.m / unyt.s**2
        out["dQ"]     = np.copy(temp) / unyt.s
        out["dDpara"] = np.copy(temp) * unyt.m / unyt.s**2
        out["clog"]   = np.copy(temp)
        out["mu0"]    = np.copy(temp)
        out["mu1"]    = np.copy(temp)
        out["dmu0"]   = np.copy(temp)

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
            self.libascot.libascot_eval_collcoefs(Neval, va, R[i], phi[i], z[i],
                                                  t[i], ma, qa, F, Dpara, Dperp,
                                                  K, nu, Q, dQ, dDpara, clog,
                                                  mu0, mu1, dmu0)
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
        assert self.bfield_initialized is not None, \
            "Magnetic field not initialized"

        rhovals = np.asarray(rhovals).ravel().astype(dtype="f8")

        ngrid = 100
        r   = np.zeros((ngrid,), dtype="f8")
        z   = np.zeros((ngrid,), dtype="f8")
        rho = np.zeros((ngrid,), dtype="f8")

        self.libascot.libascot_B_field_eval_rhovals(
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
