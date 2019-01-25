"""
Python side of the interactive Python interface.

This module defines Ascotpy class whose methods can be used in Python to call
Ascot5 functions (written in C) directly. The callable functions are defined in
library module ascotpy.c which must be compiled first with make ascotpy. This
module acts as a wrapper for those functions.

File: ascotpymod.py
"""
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

class Ascotpy:
    """
    An object representing a running ascot5 process.
    """

    def __init__(self, libpath, h5fn):
        """
        Initialize and start Ascot5 process using given HDF5 file as an input.

        This function starts the process but does not read or initialize any
        data yet. All flags are set to False and function arguments are defined
        explicitly (because C uses strong typing whereas Python does not).
        """
        self.ascotlib = ctypes.CDLL(libpath)
        self.h5fn = h5fn.encode('UTF-8')

        self.bfield_initialized  = False
        self.efield_initialized  = False
        self.plasma_initialized  = False
        self.neutral_initialized = False
        self.wall_initialized    = False

        real_p = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")

        # Init and free functions.
        fun = self.ascotlib.ascotpy_init
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
                        ctypes.c_int, ctypes.c_int, ctypes.c_int]

        fun = self.ascotlib.ascotpy_free
        fun.restype  = None
        fun.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                        ctypes.c_int, ctypes.c_int]

        # B field functions.
        fun = self.ascotlib.ascotpy_B_field_eval_B_dB
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p]

        fun = self.ascotlib.ascotpy_B_field_eval_psi
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p]

        fun = self.ascotlib.ascotpy_B_field_eval_rho
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p]

        fun = self.ascotlib.ascotpy_B_field_get_axis
        fun.restype  = None
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p]

        # E field functions.
        fun = self.ascotlib.ascotpy_E_field_eval_E
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p,
                        real_p, real_p, real_p]

        # Plasma functions.
        fun = self.ascotlib.ascotpy_plasma_get_n_species
        fun.restype  = ctypes.c_int

        fun = self.ascotlib.ascotpy_plasma_get_species_mass_and_charge
        fun.restype  = None
        fun.argtypes = [real_p, real_p]

        fun = self.ascotlib.ascotpy_plasma_eval_background
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p, real_p]

        # Neutral functions.
        fun = self.ascotlib.ascotpy_neutral_eval_density
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p, real_p]

    def ascotpy_reload(self, h5fn):
        """
        Change HDF5 file and free resources from old one.
        """
        self.ascotpy_free(bfield=self.bfield_initialized,
                          efield=self.efield_initialized,
                          plasma=self.plasma_initialized,
                          wall=self.wall_initialized,
                          neutral=self.neutral_initialized)
        self.h5fn = h5fn.encode('UTF-8')


    def ascotpy_init(self, bfield=False, efield=False, plasma=False, wall=False,
                     neutral=False):
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

        Raises:
            RuntimeError if initialization failed.
        """
        if bfield:
            if self.ascotlib.ascotpy_init(self.h5fn, 1, 0, 0, 0, 0) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.bfield_initialized  = True

        if efield:
            if self.ascotlib.ascotpy_init(self.h5fn, 0, 1, 0, 0, 0) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.efield_initialized  = True

        if plasma:
            if self.ascotlib.ascotpy_init(self.h5fn, 0, 0, 1, 0, 0) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.plasma_initialized  = True

        if wall:
            if self.ascotlib.ascotpy_init(self.h5fn, 0, 0, 0, 1, 0) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.wall_initialized  = True

        if neutral:
            if self.ascotlib.ascotpy_init(self.h5fn, 0, 0, 0, 0, 1) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.neutral_initialized  = True


    def ascotpy_free(self, bfield=False, efield=False, plasma=False, wall=False,
                     neutral=False):
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
        """
        if bfield:
            self.ascotlib.ascotpy_free(1, 0, 0, 0, 0)
            self.bfield_initialized  = False

        if efield:
            self.ascotlib.ascotpy_free(0, 1, 0, 0, 0)
            self.bfield_initialized  = False

        if plasma:
            self.ascotlib.ascotpy_free(0, 0, 1, 0, 0)
            self.bfield_initialized  = False

        if wall:
            self.ascotlib.ascotpy_free(0, 0, 0, 1, 0)
            self.bfield_initialized  = False

        if neutral:
            self.ascotlib.ascotpy_free(0, 0, 0, 0, 1)
            self.bfield_initialized  = False


    def ascotpy_eval_bfield(self, R, phi, z, evalb=False, evalpsi=False,
                            evalrho=False, evalaxis=False):
        """
        Evaluate magnetic field quantities at given coordinates.

        Args:
            R : array_like <br>
                Vector of R coordinates where field is evaluated [m].
            phi : array_like <br>
                Vector of phi coordinates where field is evaluated [rad].
            z : array_like <br>
                Vector of z coordinates where field is evaluated [m].
            evalb : bool, optional <br>
                Evaluate magnetic field vector and derivatives.
            evalpsi : bool, optional <br>
                Evaluate poloidal flux.
            evalrho : bool, optional <br>
                Evaluate normalized poloidal flux.
            evalaxis : bool, optional <br>
                Evaluate magnetic axis.

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert(self.bfield_initialized)

        R   = R.astype(dtype="f8")
        phi = phi.astype(dtype="f8")
        z   = z.astype(dtype="f8")

        Neval = R.size
        out = {}

        if evalb:
            out["br"]       = np.zeros(R.shape, dtype="f8")
            out["bphi"]     = np.zeros(R.shape, dtype="f8")
            out["bz"]       = np.zeros(R.shape, dtype="f8")
            out["brdr"]     = np.zeros(R.shape, dtype="f8")
            out["brdphi"]   = np.zeros(R.shape, dtype="f8")
            out["brdz"]     = np.zeros(R.shape, dtype="f8")
            out["bphidr"]   = np.zeros(R.shape, dtype="f8")
            out["bphidphi"] = np.zeros(R.shape, dtype="f8")
            out["bphidz"]   = np.zeros(R.shape, dtype="f8")
            out["bzdr"]     = np.zeros(R.shape, dtype="f8")
            out["bzdphi"]   = np.zeros(R.shape, dtype="f8")
            out["bzdz"]     = np.zeros(R.shape, dtype="f8")

            self.ascotlib.ascotpy_B_field_eval_B_dB(
                Neval, R, phi, z,
                out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"],
                out["bphidr"], out["bphidphi"], out["bphidz"],
                out["bzdr"], out["bzdphi"], out["bzdz"])

        if evalpsi:
            out["psi"] = np.zeros(R.shape, dtype="f8")
            self.ascotlib.ascotpy_B_field_eval_psi(Neval, R, phi, z,out["psi"])

        if evalrho:
            out["rho"] = np.zeros(R.shape, dtype="f8")
            self.ascotlib.ascotpy_B_field_eval_rho(Neval, R, phi, z,out["rho"])

        if evalaxis:
            out["axisr"] = np.zeros(R.shape, dtype="f8")
            out["axisz"] = np.zeros(R.shape, dtype="f8")
            self.ascotlib.ascotpy_B_field_get_axis(Neval, phi, out["axisr"],
                                                   out["axisz"])

        return out


    def ascotpy_eval_efield(self, R, phi, z):
        """
        Evaluate electric field quantities at given coordinates.

        Note that magnetic field has to be initialized as well.

        Args:
            R : array_like <br>
                Vector of R coordinates where field is evaluated [m].
            phi : array_like <br>
                Vector of phi coordinates where field is evaluated [rad].
            z : array_like <br>
                Vector of z coordinates where field is evaluated [m].

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert(self.bfield_initialized and self.efield_initialized)

        R   = R.astype(dtype="f8")
        phi = phi.astype(dtype="f8")
        z   = z.astype(dtype="f8")

        Neval = R.size
        out = {}
        out["er"]   = np.zeros(R.shape, dtype="f8")
        out["ephi"] = np.zeros(R.shape, dtype="f8")
        out["ez"]   = np.zeros(R.shape, dtype="f8")

        self.ascotlib.ascotpy_E_field_eval_E(
                Neval, R, phi, z, out["er"], out["ephi"], out["ez"])

        return out


    def ascotpy_eval_plasma(self, R, phi, z):
        """
        Evaluate plasma quantities at given coordinates.

        Note that magnetic field has to be initialized as well.

        Args:
            R : array_like <br>
                Vector of R coordinates where data is evaluated [m].
            phi : array_like <br>
                Vector of phi coordinates where data is evaluated [rad].
            z : array_like <br>
                Vector of z coordinates where data is evaluated [m].

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert(self.bfield_initialized and self.plasma_initialized)

        R   = R.astype(dtype="f8")
        phi = phi.astype(dtype="f8")
        z   = z.astype(dtype="f8")

        # First get background species info.
        out = {}
        out["n_species"] = self.ascotlib.ascotpy_plasma_get_n_species()
        out["mass"]      = np.zeros((out["n_species"],), dtype="f8")
        out["charge"]    = np.zeros((out["n_species"],), dtype="f8")
        self.ascotlib.ascotpy_plasma_get_species_mass_and_charge(
            out["mass"], out["charge"])

        Neval = R.size

        # Allocate enough space for electrons and all ion species.
        rawdens = np.zeros((Neval*(out["n_species"]),), dtype="f8")
        rawtemp = np.zeros((Neval*(out["n_species"]),), dtype="f8")

        self.ascotlib.ascotpy_plasma_eval_background(
                Neval, R, phi, z, rawdens, rawtemp)

        out["ne"] = rawdens[0:Neval]
        out["Te"] = rawtemp[0:Neval]
        for i in range(1, out["n_species"]):
            out["ni"+str(i)] = rawdens[(Neval)*i:(Neval)*(i+1)]
            out["Ti"+str(i)] = rawtemp[(Neval)*i:(Neval)*(i+1)]

        return out

    def ascotpy_eval_neutral(self, R, phi, z):
        """
        Evaluate plasma quantities at given coordinates.

        Args:
            R : array_like <br>
                Vector of R coordinates where data is evaluated [m].
            phi : array_like <br>
                Vector of phi coordinates where data is evaluated [rad].
            z : array_like <br>
                Vector of z coordinates where data is evaluated [m].

        Returns:
            Dictionary containing evaluated quantities.

        Raises:
            AssertionError if this is called data uninitialized.
            RuntimeError if evaluation failed.
        """
        assert(self.neutral_initialized)

        R   = R.astype(dtype="f8")
        phi = phi.astype(dtype="f8")
        z   = z.astype(dtype="f8")

        Neval = R.size
        out = {}
        out["n0"]   = np.zeros(R.shape, dtype="f8")

        self.ascotlib.ascotpy_neutral_eval_density(Neval, R, phi, z, out["n0"])

        return out

if __name__ == '__main__':
    # For testing purposes.
    import os
    ascot = Ascotpy(os.path.abspath("ascotpy.so"), "ascot.h5")
    ascot.ascotpy_init(bfield=True, efield=True, plasma=True, wall=True,
                       neutral=True)

    R   = np.array([6.2,   7, 8])
    phi = np.array([  0,   0, 0])
    z   = np.array([0.0, 0.2, 0.2])

    bvals       = ascot.ascotpy_eval_bfield(R, phi, z, evalb=True, evalpsi=True,
                                            evalrho=True, evalaxis=True)
    evals       = ascot.ascotpy_eval_efield(R, phi, z)
    plasmavals  = ascot.ascotpy_eval_plasma(R, phi, z)
    neutralvals = ascot.ascotpy_eval_neutral(R, phi, z)

    print(bvals)

    ascot.ascotpy_free(bfield=True, efield=True, plasma=True, wall=True,
                       neutral=True)
