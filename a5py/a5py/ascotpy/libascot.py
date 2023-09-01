"""Python side of the interactive Python interface.

This module defines LibAscot class whose methods can be used in Python to call
Ascot5 functions (written in C) directly. The callable functions are defined in
library module libascot.c which must be compiled first with make libascot. This
module acts as a wrapper for those functions. More advanced functionality should
be implemented in other modules.
"""
import ctypes
import sys
import os
import warnings
import numpy as np
import unyt

from scipy import interpolate
from numpy.ctypeslib import ndpointer

from a5py.physlib import parseunits

try:
    from . import ascot2py
    PTR_REAL = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")
    PTR_LONG = ndpointer(ctypes.c_long, flags="C_CONTIGUOUS")
    PTR_SIM  = ctypes.POINTER(ascot2py.struct_c__SA_sim_offload_data)
    PTR_ARR  = ctypes.POINTER(ctypes.c_double)
    _LIBASCOT = ascot2py._libraries['libascot.so']
except OSError as error:
    _LIBASCOT = None
    msg = """
    Warning: Failed to import libascot.so. Some functionalities of Ascot
    are not available. Verify that libascot.so has been compiled, it can be
    found in LD_LIBRARY_PATH, and dependencies (e.g. HDF5) are available.
    """
    warnings.warn(msg)

class LibAscot:
    """Python wrapper of libascot.so.
    """

    def _eval_bfield(self, r, phi, z, t, evalb=False, evalrho=False,
                     evalaxis=False):
        """Evaluate magnetic field quantities at given coordinates.

        Parameters
        ----------
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].
        evalb : bool, optional
            Evaluate B_field_eval_B_dB.
        evalrho : bool, optional
            Evaluate B_field_eval_rho.
        evalaxis : bool, optional
            Evaluate B_field_get_axis.

        Returns
        -------
        out : dict [str, array_like (n,)]
            Evaluated quantities in a dictionary.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield")
        Neval = r.size
        out = {}

        if evalb:
            Tperm = unyt.T / unyt.m
            out["br"]       = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.T
            out["bphi"]     = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.T
            out["bz"]       = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.T
            out["brdr"]     = (np.zeros(r.shape, dtype="f8") + np.nan) * Tperm
            out["brdphi"]   = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.T
            out["brdz"]     = (np.zeros(r.shape, dtype="f8") + np.nan) * Tperm
            out["bphidr"]   = (np.zeros(r.shape, dtype="f8") + np.nan) * Tperm
            out["bphidphi"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.T
            out["bphidz"]   = (np.zeros(r.shape, dtype="f8") + np.nan) * Tperm
            out["bzdr"]     = (np.zeros(r.shape, dtype="f8") + np.nan) * Tperm
            out["bzdphi"]   = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.T
            out["bzdz"]     = (np.zeros(r.shape, dtype="f8") + np.nan) * Tperm

            fun = _LIBASCOT.libascot_B_field_eval_B_dB
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                Neval, r, phi, z, t, out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"], out["bphidr"],
                out["bphidphi"], out["bphidz"], out["bzdr"], out["bzdphi"],
                out["bzdz"])

        if evalrho:
            out["psi"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.Wb
            out["rho"] = (np.zeros(r.shape, dtype="f8") + np.nan) \
                * unyt.dimensionless

            fun = _LIBASCOT.libascot_B_field_eval_rho
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                Neval, r, phi, z, t, out["rho"], out["psi"])

        if evalaxis:
            out["axisr"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.m
            out["axisz"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.m

            fun = _LIBASCOT.libascot_B_field_get_axis
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                Neval, phi, out["axisr"], out["axisz"])

        return out

    def _eval_efield(self, r, phi, z, t):
        """Evaluate electric field quantities at given coordinates.

        Parameters
        ----------
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].

        Returns
        -------
        out : dict [str, array_like (n,)]
            Evaluated quantities in a dictionary.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield", "efield")
        Neval = r.size
        out = {}
        Vperm = unyt.V / unyt.m
        out["er"]   = (np.zeros(r.shape, dtype="f8") + np.nan) * Vperm
        out["ephi"] = (np.zeros(r.shape, dtype="f8") + np.nan) * Vperm
        out["ez"]   = (np.zeros(r.shape, dtype="f8") + np.nan) * Vperm

        fun = _LIBASCOT.libascot_E_field_eval_E
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim), self._efield_offload_array,
            self._bfield_offload_array,
            Neval, r, phi, z, t, out["er"], out["ephi"], out["ez"])

        return out

    def _eval_plasma(self, r, phi, z, t):
        """Evaluate plasma quantities at given coordinates.

        Parameters
        ----------
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].

        Returns
        -------
        out : dict [str, array_like (n,)]
            Evaluated quantities in a dictionary.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield", "plasma")
        Neval = r.size

        # Allocate enough space for electrons and all ion species.
        m3 = unyt.m**3
        eV = unyt.eV
        nspecies  = self.input_getplasmaspecies()[0] + 1
        rawdens = (np.zeros((Neval*nspecies,), dtype="f8") + np.nan) / m3
        rawtemp = (np.zeros((Neval*nspecies,), dtype="f8") + np.nan) * eV

        fun = _LIBASCOT.libascot_plasma_eval_background
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim), self._bfield_offload_array,
            self._plasma_offload_array,
            Neval, r, phi, z, t, rawdens, rawtemp)

        out = {}
        out["ne"] = rawdens[0:Neval]
        out["te"] = rawtemp[0:Neval]
        for i in range(1, nspecies):
            out["ni"+str(i)] = rawdens[(Neval)*i:(Neval)*(i+1)]
            out["ti"+str(i)] = rawtemp[(Neval)*i:(Neval)*(i+1)]

        return out

    def _eval_neutral(self, r, phi, z, t):
        """Evaluate neutral input quantities at given coordinates.

        Parameters
        ----------
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].

        Returns
        -------
        out : dict [str, array_like (n,)]
            Evaluated quantities in a dictionary.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield", "neutral")
        Neval = r.size
        out = {}
        m3 = unyt.m**3
        out["n0"] = (np.zeros(r.shape, dtype="f8") + np.nan) / m3

        fun = _LIBASCOT.libascot_neutral_eval_density
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL]

        fun(ctypes.byref(self._sim), self._neutral_offload_array,
            Neval, r, phi, z, t, out["n0"])

        return out

    def _eval_boozer(self, r, phi, z, t, evalfun=False):
        """Evaluate boozer coordinate quantities at given coordinates.

        Parameters
        ----------
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].

        Returns
        -------
        out : dict [str, array_like (n,)]
            Evaluated quantities in a dictionary.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield", "boozer")
        Neval = r.size
        out = {}
        temp = np.zeros(r.shape, dtype="f8")
        nodim = unyt.dimensionless
        T = unyt.T; Wb = unyt.Wb; m = unyt.m; rad = unyt.rad
        if evalfun:
            out["qprof"]   = (np.copy(temp) + np.nan) * nodim
            out["bjac"]    = (np.copy(temp) + np.nan) / T**2
            out["bjacxb2"] = (np.copy(temp) + np.nan) * nodim

            fun = _LIBASCOT.libascot_boozer_eval_fun
            fun.restype  = ctypes.c_int
            fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._bfield_offload_array,
                self._boozer_offload_array,
                Neval, r, phi, z, t, out["qprof"], out["bjac"],
                out["bjacxb2"])
        else:
            out["psi (bzr)"]      = (np.copy(temp) + np.nan) * Wb
            out["theta"]          = (np.copy(temp) + np.nan) * rad
            out["zeta"]           = (np.copy(temp) + np.nan) * rad
            out["dpsidr (bzr)"]   = (np.copy(temp) + np.nan) * Wb/m
            out["dpsidphi (bzr)"] = (np.copy(temp) + np.nan) * Wb
            out["dpsidz (bzr)"]   = (np.copy(temp) + np.nan) * Wb/m
            out["dthetadr"]       = (np.copy(temp) + np.nan) / m
            out["dthetadphi"]     = (np.copy(temp) + np.nan) * nodim
            out["dthetadz"]       = (np.copy(temp) + np.nan) / m
            out["dzetadr"]        = (np.copy(temp) + np.nan) / m
            out["dzetadphi"]      = (np.copy(temp) + np.nan) * nodim
            out["dzetadz"]        = (np.copy(temp) + np.nan) / m
            out["rho (bzr)"]      = (np.copy(temp) + np.nan) * nodim

            fun = _LIBASCOT.libascot_boozer_eval_psithetazeta
            fun.restype  = ctypes.c_int
            fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), self._boozer_offload_array,
                self._bfield_offload_array,
                Neval, r, phi, z, t, out["psi (bzr)"], out["theta"],
                out["zeta"], out["dpsidr (bzr)"], out["dpsidphi (bzr)"],
                out["dpsidz (bzr)"], out["dthetadr"], out["dthetadphi"],
                out["dthetadz"], out["dzetadr"], out["dzetadphi"],
                out["dzetadz"], out["rho (bzr)"])

        return out

    def _eval_mhd(self, r, phi, z, t, evalpot=False):
        """Evaluate MHD quantities at given coordinates.

        Parameters
        ----------
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].
        evalpot : bool, optional
            Evaluate eigenfunctions as well.

        Returns
        -------
        out : dict [str, array_like (n,)]
            Evaluated quantities in a dictionary.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield", "boozer", "mhd")
        Neval = r.size
        out = {}
        T = unyt.T; V = unyt.V; m = unyt.m; s = unyt.s
        temp = np.zeros(r.shape, dtype="f8")
        out["mhd_br"]   = (np.copy(temp) + np.nan) * T
        out["mhd_bphi"] = (np.copy(temp) + np.nan) * T
        out["mhd_bz"]   = (np.copy(temp) + np.nan) * T
        out["mhd_er"]   = (np.copy(temp) + np.nan) * V/m
        out["mhd_ephi"] = (np.copy(temp) + np.nan) * V/m
        out["mhd_ez"]   = (np.copy(temp) + np.nan) * V/m
        out["mhd_phi"]  = (np.copy(temp) + np.nan) * V

        fun = _LIBASCOT.libascot_mhd_eval_perturbation
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim), self._bfield_offload_array,
            self._boozer_offload_array, self._mhd_offload_array,
            Neval, r, phi, z, t, out["mhd_br"], out["mhd_bphi"], out["mhd_bz"],
            out["mhd_er"], out["mhd_ephi"], out["mhd_ez"], out["mhd_phi"])

        if evalpot:
            out["alphaeig"] = (np.copy(temp) + np.nan) * T*m
            out["phieig"]   = (np.copy(temp) + np.nan) * V
            out["dadt"]     = (np.copy(temp) + np.nan) * T*m/s
            out["dadr"]     = (np.copy(temp) + np.nan) * T
            out["dadphi"]   = (np.copy(temp) + np.nan) * T*m
            out["dadz"]     = (np.copy(temp) + np.nan) * T
            out["dphidt"]   = (np.copy(temp) + np.nan) * V/s
            out["dphidr"]   = (np.copy(temp) + np.nan) * V/m
            out["dphidphi"] = (np.copy(temp) + np.nan) * V
            out["dphidz"]   = (np.copy(temp) + np.nan) * V/m

            fun = _LIBASCOT.libascot_mhd_eval
            fun.restype  = ctypes.c_int
            fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR, PTR_ARR,
                            ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL]

            fun(ctypes.byref(self._sim), self._boozer_offload_array,
                self._mhd_offload_array, self._bfield_offload_array,
                Neval, r, phi, z, t, out["alphaeig"],
                out["dadr"], out["dadphi"], out["dadz"], out["dadt"],
                out["phieig"], out["dphidr"], out["dphidphi"], out["dphidz"],
                out["dphidt"])

        return out

    @parseunits(ma="kg", qa="C", r="m", phi="rad", z="m", t="s", va="m/s")
    def input_eval_collcoefs(self, ma, qa, r, phi, z, t, va):
        """Evaluate Coulomb collision coefficients for a given test particle.

        Collision coefficients are evaluated by interpolating the plasma
        parameters on given coordinates. The coefficients are returned as
        a function of the test particle velocity.

        Parameters
        ----------
        ma : float
            Test particle mass.
        qa : float
            Test particle charge.
        r : array_like, (n,)
            R coordinates where data is evaluated.
        phi : array_like, (n,)
            phi coordinates where data is evaluated.
        z : array_like, (n,)
            z coordinates where data is evaluated.
        t : array_like, (n,)
            Time coordinates where data is evaluated.
        va : array_like, (nv,)
            Test particle velocity.

        Returns
        -------
        out : dict [str, array_like (nion+1, n, nv)]
            Evaluated collision coefficients in a dictionary.

            The first dimension is the background plasma species where the first
            index is for the electrons followed by ions.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield", "plasma")
        ma  = float(ma)
        qa  = float(qa)
        Neval = va.size

        n_species = self.input_getplasmaspecies()[0] + 1

        out = {}
        temp = np.zeros((n_species, r.size, va.size), dtype="f8")
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

        fun = _LIBASCOT.libascot_eval_collcoefs
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

        for i in range(r.size):
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
                va.size, va, r[i], phi[i], z[i], t[i],
                ma, qa, F, Dpara, Dperp,
                K, nu, Q, dQ, dDpara, clog, mu0, mu1, dmu0)
            out["F"][:,i,:]      = F[:,:]
            out["Dpara"][:,i,:]  = Dpara[:,:]
            out["Dperp"][:,i,:]  = Dperp[:,:]
            out["K"][:,i,:]      = K[:,:]
            out["nu"][:,i,:]     = nu[:,:]
            out["Q"][:,i,:]      = Q[:,:]
            out["dQ"][:,i,:]     = dQ[:,:]
            out["dDpara"][:,i,:] = dDpara[:,:]
            out["clog"][:,i,:]   = clog[:,:]
            out["mu0"][:,i,:]    = mu0[:,:]
            out["mu1"][:,i,:]    = mu1[:,:]
            out["dmu0"][:,i,:]   = dmu0[:,:]

        return out

    @parseunits(ma="kg", r="m", phi="rad", z="m", t="s", va="m/s")
    def input_eval_atomicsigma(self, ma, anum, znum, r, phi, z, t, va, ion,
                               reaction):
        """Evaluate atomic reaction rate cross-sections for a given test
        particle.

        This function is a work in progress.

        Parameters
        ----------
        ma : float
            Test particle mass
        anum : int
            Test particle atomic mass number.
        znum : int
            Test particle charge number.
        r : array_like, (n,)
            R coordinates where data is evaluated [m].
        phi : array_like (n,)
            phi coordinates where data is evaluated [rad].
        z : array_like (n,)
            z coordinates where data is evaluated [m].
        t : array_like (n,)
            Time coordinates where data is evaluated [s].
        va : array_like (n,)
            Test particle velocities where data is evaluated in each grid point.
        ion : int
            Index number of the background ion species in plasma input.
        reaction : int
            Reaction.

        Returns
        -------
        sigmav : array_like, (n,nv)
            Reaction cross-section.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        """
        self._requireinit("bfield", "plasma", "neutral", "asigma")

        fun = _LIBASCOT.libascot_eval_sigmav
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_ARR, PTR_ARR, PTR_ARR,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        ctypes.c_int, PTR_REAL, ctypes.c_int, ctypes.c_int,
                        ctypes.c_double, ctypes.c_int, ctypes.c_int, PTR_REAL]

        Neval = r.size
        Nv    = va.size

        out = {}
        out["sigmav"] = (np.zeros(r.shape, dtype="f8") + np.nan) / unyt.m**2
        fun(ctypes.byref(self._sim), self._bfield_offload_array,
            self._plasma_offload_array, self._neutral_offload_array,
            self._asigma_offload_array, Neval, r, phi, z, t, Nv, va,
            anum, znum, ma, ion, reaction, out["sigmav"])

        return out["sigmav"]


    def input_getplasmaspecies(self):
        """Get species present in plasma input (electrons excluded).

        Returns
        -------
        nspecies : int
            Number of species.
        mass : array_like
            Species' mass.
        charge : array_like
            Species' charge state.
        anum : array_like
            Species' atomic mass number.
        znum : array_like
            Species' charge number.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        fun = _LIBASCOT.libascot_plasma_get_n_species
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_ARR]

        out = {}
        out["nspecies"] = fun(
            ctypes.byref(self._sim), self._plasma_offload_array) - 1

        out["mass"]   = np.zeros((out["nspecies"],), dtype="f8")
        out["charge"] = np.zeros((out["nspecies"],), dtype="f8")
        out["anum"]   = np.zeros((out["nspecies"],), dtype="i8")
        out["znum"]   = np.zeros((out["nspecies"],), dtype="i8")

        fun = _LIBASCOT.libascot_plasma_get_species_mass_and_charge
        fun.restype  = None
        fun.argtypes = [PTR_SIM, PTR_ARR, PTR_REAL, PTR_REAL, PTR_LONG,
                        PTR_LONG]
        fun(ctypes.byref(self._sim), self._plasma_offload_array,
            out["mass"], out["charge"], out["anum"], out["znum"])

        return out["nspecies"], out["mass"], out["charge"], out["anum"],\
            out["znum"]

    @parseunits(rhovals="1", theta="rad", phi="deg", time="s")
    def input_rhotheta2rz(self, rhovals, theta, phi, time):
        """Convert rho coordinate to (R,z) position.

        Conversion is done at fixed (theta, phi).

        Parameters
        ----------
        rhovals : array_like (n,)
            Normalized poloidal flux coordinates to be converted.
        theta : float
            Poloidal angle.
        phi : float
            Toroidal angle.
        time : float
            Time slice.

        Returns
        -------
        r : array_like (n,)
            R coordinates at (rhovals, theta, phi)
        z : array_like (n,)
            z coordinates at (rhovals, theta, phi)

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield")
        rhovals = np.asarray(rhovals).ravel().astype(dtype="f8")

        ngrid = 1000
        r   = np.zeros((ngrid,), dtype="f8")
        z   = np.zeros((ngrid,), dtype="f8")
        rho = np.zeros((ngrid,), dtype="f8")

        fun = _LIBASCOT.libascot_B_field_eval_rhovals
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

        return (fr(rhovals) * unyt.m, fz(rhovals) * unyt.m)
