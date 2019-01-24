"""
Python side of the interactive Python interface.

This module defines Ascotpy class whose methods can be used in Python to call
Ascot5 functions (written in C) directly. The callable functions are defined in
library module ascotpy.c which must be compiled first with make ascotpy. This
module acts as a wrapper for those functions.

File: ascotpy.py
"""
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

class Ascotpy:
    """
    An object representing a running ascot5 process.
    """

    def __init__(self, libpath, ascotfn):
        """
        Initialize and start Ascot5 process using given HDF5 file as an input.

        This function starts the process but does not read or initialize any
        data yet. All flags are set to False and function arguments are defined
        explicitly (because C uses strong typing whereas Python does not).
        """
        self.ascotlib = ctypes.CDLL(libpath)
        self.ascotfn = ascotfn.encode('UTF-8')

        self.bfield_initialized  = False
        self.efield_initialized  = False
        self.plasma_initialized  = False
        self.neutral_initialized = False
        self.wall_initialized    = False

        real_p = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")

        fun = self.ascotlib.ascotpy_init
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
                        ctypes.c_int, ctypes.c_int, ctypes.c_int]

        fun = self.ascotlib.ascotpy_free
        fun.restype  = None
        fun.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                        ctypes.c_int, ctypes.c_int]

        fun = self.ascotlib.ascotpy_bfield_eval_B_dB
        fun.restype  = ctypes.c_int
        fun.argtypes = [ctypes.c_int, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p,
                        real_p, real_p, real_p, real_p, real_p, real_p]

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
            if self.ascotlib.ascotpy_init(self.ascotfn, 1, 0, 0, 0, 0) :
                raise RuntimeError("Failed to initialize magnetic field")

            self.bfield_initialized  = True

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
            AssertionError if this is called magnetic field uninitialized.
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

            self.ascotlib.ascotpy_bfield_eval_B_dB(
                Neval, R, phi, z,
                out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"],
                out["bphidr"], out["bphidphi"], out["bphidz"],
                out["bzdr"], out["bzdphi"], out["bzdz"])

        return out

if __name__ == '__main__':
    # For testing purposes.
    import os
    ascot = Ascotpy(os.path.abspath("ascotpy.so"), "ascot.h5")
    ascot.ascotpy_init(bfield=True)

    R   = np.array([6.2,   7, 8])
    phi = np.array([  0,   0, 0])
    z   = np.array([0.0, 0.2, 0.2])

    bvals = ascot.ascotpy_eval_bfield(R, phi, z, evalb=True)

    ascot.ascotpy_free(bfield=True)
