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

from ctypes.util import find_library
from numpy.ctypeslib import ndpointer
from a5py.ascot5io.ascot5 import Ascot

class LibAscot:
    """
    An object representing a running ascot5 process.
    """

    def __init__(self, h5fn="ascot.h5", libpath="libascot.so"):
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
            raise type(e)(str(e) + msg).with_traceback(sys.exc_info()[2])

        # Check that the HDF5 file exists
        self.h5fn = h5fn.encode('UTF-8')
        Ascot(self.h5fn)

        # Initialize attributes
        self.bfield_initialized  = False
        self.efield_initialized  = False
        self.plasma_initialized  = False
        self.neutral_initialized = False
        self.wall_initialized    = False
        self.boozer_initialized  = False
        self.mhd_initialized     = False

        # Declare functions found in libascot

        real_p = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")

        # Init and free functions.
        try:
            fun = self.libascot.libascot_init
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            ctypes.c_int, ctypes.c_int]
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
                            real_p, real_p, real_p, real_p, real_p, real_p]
        except AttributeError:
            warnings.warn("libascot_boozer_eval_psithetazeta not found", Warning)
            pass

        try:
            fun = self.libascot.libascot_mhd_eval_perturbation
            fun.restype  = ctypes.c_int
            fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p,
                            real_p, real_p, real_p, real_p, real_p, real_p]
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


    def reload(self, h5fn):
        """
        Change HDF5 file and free resources from old one.

        Args:
            h5fn : str <br>
                Name of the new HDF5 file.
        """
        self.free(bfield=self.bfield_initialized,
                  efield=self.efield_initialized,
                  plasma=self.plasma_initialized,
                  wall=self.wall_initialized,
                  neutral=self.neutral_initialized)
        self.h5fn = h5fn.encode('UTF-8')
        Ascot(self.h5fn)


    def init(self, bfield=False, efield=False, plasma=False, wall=False,
             neutral=False, boozer=False, mhd=False):
        """
        Initialize input data.

        Args:
            bfield : bool, optional <br>
                Flag for initializing magnetic field.
            efield : bool, optional <br>
                Flag for initializing electric field.
            plasma : bool, optional <br>
                Flag for initializing plasma data.
            wall : bool, optional <br>
                Flag for initializing wall data.
            neutral : bool, optional <br>
                Flag for initializing neutral data.
            boozer : bool, optional <br>
                Flag for initializing boozer data.
            mhd : bool, optional <br>
                Flag for initializing mhd data.

        Raises:
            RuntimeError if initialization failed.
        """
        if bfield and self.bfield_initialized:
            warnings.warn("Magnetic field already initialized.", Warning)
        elif bfield:
            if self.libascot.libascot_init(self.h5fn, 1, 0, 0, 0, 0, 0, 0) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.bfield_initialized = True

        if efield and self.efield_initialized:
            warnings.warn("Electric field already initialized.", Warning)
        elif efield:
            if self.libascot.libascot_init(self.h5fn, 0, 1, 0, 0, 0, 0, 0) :
                raise RuntimeError("Failed to initialize electric field")

            self.efield_initialized = True

        if plasma and self.plasma_initialized:
            warnings.warn("Plasma already initialized.", Warning)
        elif plasma:
            if self.libascot.libascot_init(self.h5fn, 0, 0, 1, 0, 0, 0, 0) :
                raise RuntimeError("Failed to initialize plasma")

            self.plasma_initialized = True

        if wall and self.wall_initialized:
            warnings.warn("Wall already initialized.", Warning)
        elif wall:
            if self.libascot.libascot_init(self.h5fn, 0, 0, 0, 1, 0, 0, 0) :
                raise RuntimeError("Failed to initialize wall")

            self.wall_initialized = True

        if neutral and self.neutral_initialized:
            warnings.warn("Neutral data already initialized.", Warning)
        if neutral:
            if self.libascot.libascot_init(self.h5fn, 0, 0, 0, 0, 1, 0, 0) :
                raise RuntimeError("Failed to initialize neutral data")

            self.neutral_initialized = True

        if boozer and self.boozer_initialized:
            warnings.warn("Boozer data already initialized.", Warning)
        if boozer:
            if self.libascot.libascot_init(self.h5fn, 0, 0, 0, 0, 0, 1, 0) :
                raise RuntimeError("Failed to initialize boozer data")

            self.boozer_initialized = True

        if mhd and self.mhd_initialized:
            warnings.warn("MHD data already initialized.", Warning)
        if mhd:
            if self.libascot.libascot_init(self.h5fn, 0, 0, 0, 0, 0, 0, 1) :
                raise RuntimeError("Failed to initialize MHD data")

            self.mhd_initialized = True


    def free(self, bfield=False, efield=False, plasma=False, wall=False,
             neutral=False, boozer=False, mhd=False):
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
        """
        if bfield and self.bfield_initialized:
            self.libascot.libascot_free(1, 0, 0, 0, 0, 0, 0)
            self.bfield_initialized = False

        if efield and self.efield_initialized:
            self.libascot.libascot_free(0, 1, 0, 0, 0, 0, 0)
            self.efield_initialized = False

        if plasma and self.plasma_initialized:
            self.libascot.libascot_free(0, 0, 1, 0, 0, 0, 0)
            self.plasma_initialized = False

        if wall and self.wall_initialized:
            self.libascot.libascot_free(0, 0, 0, 1, 0, 0, 0)
            self.wall_initialized = False

        if neutral and self.neutral_initialized:
            self.libascot.libascot_free(0, 0, 0, 0, 1, 0, 0)
            self.neutral_initialized = False

        if boozer and self.boozer_initialized:
            self.libascot.libascot_free(0, 0, 0, 0, 0, 1, 0)
            self.boozer_initialized = False

        if mhd and self.mhd_initialized:
            self.libascot.libascot_free(0, 0, 0, 0, 0, 0, 1)
            self.mhd_initialized = False


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
        assert self.bfield_initialized, "Magnetic field not initialized"

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
        assert self.bfield_initialized, "Magnetic field not initialized"
        assert self.efield_initialized, "Electric field not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["er"]   = np.zeros(R.shape, dtype="f8")
        out["ephi"] = np.zeros(R.shape, dtype="f8")
        out["ez"]   = np.zeros(R.shape, dtype="f8")

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
        assert self.plasma_initialized, "Plasma not initialized"

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
        assert self.bfield_initialized, "Magnetic field not initialized"
        assert self.plasma_initialized, "Plasma not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size

        # Allocate enough space for electrons and all ion species.
        nspecies = self.libascot.libascot_plasma_get_n_species()
        rawdens = np.zeros((Neval*nspecies,), dtype="f8")
        rawtemp = np.zeros((Neval*nspecies,), dtype="f8")

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
        Evaluate plasma quantities at given coordinates.

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
        assert self.neutral_initialized, "Neutral data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["density"] = np.zeros(R.shape, dtype="f8")

        self.libascot.libascot_neutral_eval_density(Neval, R, phi, z, t,
                                                    out["density"])

        return out


    def eval_boozer(self, R, phi, z, t):
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

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.boozer_initialized, "Boozer data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["psi"]        = np.zeros(R.shape, dtype="f8")
        out["theta"]      = np.zeros(R.shape, dtype="f8")
        out["zeta"]       = np.zeros(R.shape, dtype="f8")
        out["dpsidr"]     = np.zeros(R.shape, dtype="f8")
        out["dpsidphi"]   = np.zeros(R.shape, dtype="f8")
        out["dpsidz"]     = np.zeros(R.shape, dtype="f8")
        out["dthetadr"]   = np.zeros(R.shape, dtype="f8")
        out["dthetadphi"] = np.zeros(R.shape, dtype="f8")
        out["dthetadz"]   = np.zeros(R.shape, dtype="f8")
        out["dzetadr"]    = np.zeros(R.shape, dtype="f8")
        out["dzetadphi"]  = np.zeros(R.shape, dtype="f8")
        out["dzetadz"]    = np.zeros(R.shape, dtype="f8")

        self.libascot.libascot_boozer_eval_psithetazeta(
            Neval, R, phi, z, t, out["psi"], out["theta"], out["zeta"],
            out["dpsidr"], out["dpsidphi"], out["dpsidz"],
            out["dthetadr"], out["dthetadphi"], out["dthetadz"],
            out["dzetadr"], out["dzetadphi"], out["dzetadz"])

        return out


    def eval_mhd(self, R, phi, z, t):
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

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert self.bfield_initialized, "Magnetic field not initialized"
        assert self.boozer_initialized, "Boozer data not initialized"
        assert self.mhd_initialized,    "MHD data not initialized"

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["br"]   = np.zeros(R.shape, dtype="f8")
        out["bphi"] = np.zeros(R.shape, dtype="f8")
        out["bz"]   = np.zeros(R.shape, dtype="f8")
        out["er"]   = np.zeros(R.shape, dtype="f8")
        out["ephi"] = np.zeros(R.shape, dtype="f8")
        out["ez"]   = np.zeros(R.shape, dtype="f8")

        self.libascot.libascot_mhd_eval_perturbation(Neval, R, phi, z, t,
                                                     out["br"],   out["bphi"],
                                                     out["bz"],   out["er"],
                                                     out["ephi"], out["ez"])

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
        assert self.bfield_initialized, "Magnetic field not initialized"
        assert self.plasma_initialized, "Plasma not initialized"


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
        assert self.bfield_initialized, "Magnetic field not initialized"

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
