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

from numpy.ctypeslib import ndpointer

try:
    from . import ascot2py
    PTR_REAL = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
    PTR_SIM  = ctypes.POINTER(ascot2py.struct_c__SA_sim_offload_data)
    PTR_ARR  = ctypes.POINTER(ctypes.c_double)
    PTR_MRK  = ctypes.POINTER(ascot2py.struct_c__SA_particle_state)
    STRUCT_OFFLOAD_PACKAGE = ascot2py.struct_c__SA_offload_package
    STRUCT_OFFLOAD_DATA    = ascot2py.struct_c__SA_sim_offload_data
    LIBASCOT = ascot2py._libraries['libascot.so']
    LIBRARY_AVAILABLE = True
except OSError as error:
    LIBRARY_AVAILABLE = False
    msg = """
    Warning: Failed to import libascot.so. Some functionalities of Ascot\n
    are not available. Verify that libascot.so has been compiled, it can be\n
    found in LD_LIBRARY_PATH, and dependencies (e.g. HDF5) are available.\n
    """
    print(error.strerror)
    warnings.warn(msg)

class LibAscot:
    """Python wrapper of libascot.so.
    """

    def _eval_bfield(self, R, phi, z, t, evalb=False, evalrho=False,
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
        #assert self._sim.qid_bfield != LibAscot.DUMMY_QID, \
        #    "Magnetic field not initialized"

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

            fun = LIBASCOT.libascot_B_field_eval_B_dB
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                Neval, R, phi, z, t, out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"], out["bphidr"],
                out["bphidphi"], out["bphidz"], out["bzdr"], out["bzdphi"],
                out["bzdz"])

        if evalrho:
            out["rho"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["psi"] = np.zeros(R.shape, dtype="f8") + np.nan

            fun = LIBASCOT.libascot_B_field_eval_rho
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                Neval, R, phi, z, t, out["rho"], out["psi"])

        if evalaxis:
            out["axisr"] = np.zeros(R.shape, dtype="f8") + np.nan
            out["axisz"] = np.zeros(R.shape, dtype="f8") + np.nan

            fun = LIBASCOT.libascot_B_field_get_axis
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                Neval, phi, out["axisr"], out["axisz"])

        return out

    def _eval_efield(self, R, phi, z, t):
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

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["er"]   = np.zeros(R.shape, dtype="f8") + np.nan
        out["ephi"] = np.zeros(R.shape, dtype="f8") + np.nan
        out["ez"]   = np.zeros(R.shape, dtype="f8") + np.nan

        fun = LIBASCOT.libascot_E_field_eval_E
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim), self._efield_offload_array,
            self._bfield_offload_array,
            Neval, R, phi, z, t, out["er"], out["ephi"], out["ez"])

        return out


    def _get_plasmaspecies(self):
        """
        Get plasma species information.

        Returns;
            Dictionary containing nspecies, and anum, znum, charge, and mass for
            each species.
        """

        fun = LIBASCOT.libascot_plasma_get_n_species
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR]

        out = {}
        out["nspecies"] = fun(
            ctypes.byref(self._sim), self._plasma_offload_array)

        out["mass"]     = np.zeros((out["nspecies"],), dtype="f8")
        out["charge"]   = np.zeros((out["nspecies"],), dtype="f8")

        fun = LIBASCOT.libascot_plasma_get_species_mass_and_charge
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim), self._plasma_offload_array,
            out["mass"], out["charge"])

        return out


    def _eval_plasma(self, R, phi, z, t):
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

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size

        # Allocate enough space for electrons and all ion species.
        nspecies = LIBASCOT.libascot_plasma_get_n_species()
        rawdens = np.zeros((Neval*nspecies,), dtype="f8") + np.nan
        rawtemp = np.zeros((Neval*nspecies,), dtype="f8") + np.nan

        fun = LIBASCOT.libascot_plasma_eval_background
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim), self._bfield_offload_array,
            self._plasma_offload_array,
            Neval, R, phi, z, t, rawdens, rawtemp)

        out = {}
        out["ne"] = rawdens[0:Neval]
        out["te"] = rawtemp[0:Neval]
        for i in range(1, nspecies):
            out["ni"+str(i)] = rawdens[(Neval)*i:(Neval)*(i+1)]
            out["ti"+str(i)] = rawtemp[(Neval)*i:(Neval)*(i+1)]

        return out


    def _eval_neutral(self, R, phi, z, t):
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

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}
        out["n0"] = np.zeros(R.shape, dtype="f8") + np.nan

        fun = LIBASCOT.libascot_neutral_eval_density
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL]

        fun(ctypes.byref(self._sim), self._neutral_offload_array,
            Neval, R, phi, z, t, out["n0"])

        return out


    def _eval_boozer(self, R, phi, z, t, evalfun=False):
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

        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")

        Neval = R.size
        out = {}

        if evalfun:
            out["qprof"]   = np.zeros(R.shape, dtype="f8") + np.nan
            out["bjac"]    = np.zeros(R.shape, dtype="f8") + np.nan
            out["bjacxb2"] = np.zeros(R.shape, dtype="f8") + np.nan

            fun = LIBASCOT.libascot_boozer_eval_fun
            fun.restype  = ctypes.c_int
            fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                self._boozer_offload_array,
                Neval, R, phi, z, t, out["qprof"], out["bjac"],
                out["bjacxb2"])
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

            fun = LIBASCOT.libascot_boozer_eval_psithetazeta
            fun.restype  = ctypes.c_int
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._boozer_offload_array,
                Neval, R, phi, z, t, out["psi"], out["theta"], out["zeta"],
                out["dpsidr"], out["dpsidphi"], out["dpsidz"],
                out["dthetadr"], out["dthetadphi"], out["dthetadz"],
                out["dzetadr"], out["dzetadphi"], out["dzetadz"], out["rho"])

        return out


    def _eval_mhd(self, R, phi, z, t, evalpot=False):
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

        fun = LIBASCOT.libascot_mhd_eval_perturbation
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim), self._bfield_offload_array,
            self._boozer_offload_array, self._mhd_offload_array,
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

            fun = LIBASCOT.libascot_mhd_eval
            fun.restype  = ctypes.c_int
            fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL]

            fun(ctypes.byref(self._sim), self._boozer_offload_array,
                self._mhd_offload_array, Neval, R, phi, z, t, out["alpha"],
                out["dadr"], out["dadphi"], out["dadz"], out["dadt"],
                out["phi"], out["dphidr"], out["dphidphi"], out["dphidz"],
                out["dphidt"])

        return out


    def input_eval_collcoefs(self, ma, qa, R, phi, z, t, va):
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


        ma  = float(ma)
        qa  = float(qa)
        R   = np.asarray(R).ravel().astype(dtype="f8")
        phi = np.asarray(phi).ravel().astype(dtype="f8")
        z   = np.asarray(z).ravel().astype(dtype="f8")
        t   = np.asarray(t).ravel().astype(dtype="f8")
        va  = np.asarray(va).ravel().astype(dtype="f8")
        Neval = va.size

        n_species = \
            LIBASCOT.libascot_plasma_get_n_species(
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

        fun = LIBASCOT.libascot_eval_collcoefs
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]
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
            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                self._plasma_offload_array,
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


    def input_rhotheta2rz(self, rhovals, theta, phi, time):

        rhovals = np.asarray(rhovals).ravel().astype(dtype="f8")

        ngrid = 100
        r   = np.zeros((ngrid,), dtype="f8")
        z   = np.zeros((ngrid,), dtype="f8")
        rho = np.zeros((ngrid,), dtype="f8")

        fun = LIBASCOT.libascot_B_field_eval_rhovals
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR,
                        ctypes.c_int,    ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        PTR_REAL, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim), self._bfield_offload_array,
            ngrid, np.min(rhovals), np.max(rhovals), theta, phi, time,
            r, z, rho)

        rho[0]  = rhovals[0]
        rho[-1] = rhovals[-1]

        fr = interpolate.interp1d(rho, r, fill_value="extrapolate")
        fz = interpolate.interp1d(rho, z, fill_value="extrapolate")

        return (fr(rhovals), fz(rhovals))
