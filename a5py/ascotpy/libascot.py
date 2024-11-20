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

    def _ndpointerwithnull(*args, **kwargs):
        """This produces ndpointer which can also be a NULL pointer.

        This wrapper is required since otherwise it is not possible to
        supply None (indicating NULL) as an argument that requires
        a numpy array.
        """
        base = ndpointer(*args, **kwargs)
        def from_param(cls, obj):
            if obj is None:
                return obj
            return base.from_param(obj)
        return type(base.__name__, (base,),
                    {'from_param': classmethod(from_param)})

    PTR_REAL = _ndpointerwithnull(ctypes.c_double, flags="C_CONTIGUOUS")
    PTR_INT  = _ndpointerwithnull(ctypes.c_int,    flags="C_CONTIGUOUS")
    PTR_SIM  = ctypes.POINTER(ascot2py.struct_c__SA_sim_data)
    PTR_ARR  = ctypes.POINTER(ctypes.c_double)
    STRUCT_DIST5D        = ascot2py.struct_c__SA_dist_5D_data
    STRUCT_AFSITHERMAL   = ascot2py.struct_c__SA_afsi_thermal_data
    STRUCT_AFSIDATA      = ascot2py.struct_c__SA_afsi_data
    AFSI_REACTIONS       = ascot2py.Reaction__enumvalues
    _LIBASCOT = ascot2py._libraries['libascot.so']

except ImportError as error:
    _LIBASCOT = None
    PTR_REAL  = None
    PTR_INT   = None
    PTR_SIM   = None
    PTR_ARR   = None
    STRUCT_DIST5D        = None
    STRUCT_AFSITHERMAL   = None
    STRUCT_AFSIDATA      = None
    AFSI_REACTIONS       = None
    msg = \
        "Failed to load libascot.so: " + str(error) + "\n" \
        "Some functionalities of Ascot are not available"
    warnings.warn(msg, stacklevel=4)

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
            fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim),
                Neval, r, phi, z, t, out["br"], out["bphi"], out["bz"],
                out["brdr"], out["brdphi"], out["brdz"], out["bphidr"],
                out["bphidphi"], out["bphidz"], out["bzdr"], out["bzdphi"],
                out["bzdz"])

        if evalrho:
            out["psi"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.Wb
            out["rho"] = (np.zeros(r.shape, dtype="f8") + np.nan) \
                * unyt.dimensionless
            out["rhodpsi"] = (np.zeros(r.shape, dtype="f8") + np.nan) / unyt.Wb
            out["psidr"] = (np.zeros(r.shape, dtype="f8") + np.nan) \
                * unyt.Wb / unyt.m
            out["psidphi"] = (np.zeros(r.shape, dtype="f8") + np.nan) \
                * unyt.Wb
            out["psidz"] = (np.zeros(r.shape, dtype="f8") + np.nan) \
                * unyt.Wb / unyt.m

            fun = _LIBASCOT.libascot_B_field_eval_rho
            fun.restype  = None
            fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim),
                Neval, r, phi, z, t, out["rho"], out["rhodpsi"], out["psi"],
                out["psidr"], out["psidphi"], out["psidz"])

        if evalaxis:
            out["axisr"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.m
            out["axisz"] = (np.zeros(r.shape, dtype="f8") + np.nan) * unyt.m

            fun = _LIBASCOT.libascot_B_field_get_axis
            fun.restype  = None
            fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim), Neval, phi, out["axisr"], out["axisz"])

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
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim),
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
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim), Neval, r, phi, z, t, rawdens, rawtemp)

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
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim), Neval, r, phi, z, t, out["n0"])

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
            fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim),
                Neval, r, phi, z, t, out["qprof"], out["bjac"], out["bjacxb2"])
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
            fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim),
                Neval, r, phi, z, t, out["psi (bzr)"], out["theta"],
                out["zeta"], out["dpsidr (bzr)"], out["dpsidphi (bzr)"],
                out["dpsidz (bzr)"], out["dthetadr"], out["dthetadphi"],
                out["dthetadz"], out["dzetadr"], out["dzetadphi"],
                out["dzetadz"], out["rho (bzr)"])

        return out

    def _eval_mhd(self, r, phi, z, t, mode=None, evalpot=False):
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
        mode : int
            Pick a specific mode that is included by giving its index or
            None to include all modes.
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
        if mode is None: mode = -1 # equal to MHD_INCLUDE_ALL
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
        fun.restype  = None
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]

        fun(ctypes.byref(self._sim),
            Neval, r, phi, z, t, mode, out["mhd_br"], out["mhd_bphi"],
            out["mhd_bz"], out["mhd_er"], out["mhd_ephi"], out["mhd_ez"],
            out["mhd_phi"])

        if evalpot:
            out["alphaeig"] = (np.copy(temp) + np.nan) * m
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
            fun.restype  = None
            fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, ctypes.c_int, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                            PTR_REAL, PTR_REAL, PTR_REAL]

            fun(ctypes.byref(self._sim),
                Neval, r, phi, z, t, mode, out["alphaeig"],
                out["dadr"], out["dadphi"], out["dadz"], out["dadt"],
                out["phieig"], out["dphidr"], out["dphidphi"], out["dphidz"],
                out["dphidt"])

        return out

    def _input_mhd_modes(self):
        """Return mode numbers, amplitudes and frequencies.

        Returns
        -------
        nmodes : int
            Number of modes present.
        nmode : array_like (nmodes,)
            Toroidal mode number of each mode.
        mmode : array_like (nmodes,)
            Poloidal mode number of each mode.
        amplitude : array_like (nmodes,)
            Amplitude of each mode.
        omega : array_like (nmodes,)
            Frequency of each mode.
        phase : array_like (nmodes,)
            Phase of each mode.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("mhd")
        fun = _LIBASCOT.libascot_mhd_get_n_modes
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM]

        out = {}
        out["nmodes"] = fun(ctypes.byref(self._sim))

        out["nmode"]     = np.zeros((out["nmodes"],), dtype="i4")
        out["mmode"]     = np.zeros((out["nmodes"],), dtype="i4")
        out["amplitude"] = np.zeros((out["nmodes"],), dtype="f8")
        out["omega"] = np.zeros((out["nmodes"],), dtype="f8")
        out["phase"]     = np.zeros((out["nmodes"],), dtype="f8")

        fun = _LIBASCOT.libascot_mhd_get_mode_specs
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_INT, PTR_INT, PTR_REAL, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim), out["nmode"], out["mmode"],
            out["amplitude"], out["omega"], out["phase"])

        return out["nmodes"], out["nmode"], out["mmode"], out["amplitude"],\
            out["omega"], out["phase"]

    @parseunits(ma="kg", qa="C", r="m", phi="rad", z="m", t="s", va="m/s")
    def input_eval_collcoefs(self, ma, qa, r, phi, z, t, va, *coefs, grid=True):
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
        va : array_like, (nv,) or (n,)
            Test particle velocity.
        *coefs : str
            Names of the coefficients to be evaluated.

            "clog"   - Coulomb logarithm.
            "f"      - Drag in particle Fokker-Planck equation.
            "k"      - Drag in guiding-center Fokker-Planck equation.
            "nu"     - Pitch collision frequency.
            "dpara"  - Parallel diffusion.
            "dperp"  - Perpendicular diffusion.
            "ddpara" -  d(Dpara) / dv.
            "q"      - From K = Q + dDpara + 2*Dpara / va.
            "dq"     - d(Q) / dv.
            "mu0"    - One of the special functions needed in evaluation.
            "mu1"    - One of the special functions needed in evaluation.
            "dmu0"   - One of the special functions needed in evaluation.
        grid : bool, optional
            If True, all velocity components are evaluated at each (R,phi,z)
            position, i.e., the returned values are of shape (nion+1,n,nv).

            If False, the coefficients are evaluated at positions (R,phi,z,va),
            e.g. along an orbit, and the returned values have shape (nion+1,n).

        Returns
        -------
        *out : array_like, (nion+1, n, nv) or (nion+1, n)
            Evaluated collision coefficients in same order as declared in
            ``*coefs``.

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
        Neval = r.size
        Nv = va.size if grid else 1
        n_species = self.input_getplasmaspecies()[0]

        m = unyt.m; s = unyt.s
        out = {"clog":None, "f":None, "k":None, "nu":None, "dpara":None,
               "dperp":None, "ddpara":None, "q":None, "dq":None, "mu0":None,
               "mu1":None, "dmu0":None}
        temp = np.zeros((n_species, Neval, Nv), dtype="f8")
        if "clog"   in coefs: out["clog"]   = np.copy(temp)
        if "f"      in coefs: out["f"]      = np.copy(temp) * m/s**2
        if "k"      in coefs: out["k"]      = np.copy(temp) * m/s**2
        if "nu"     in coefs: out["nu"]     = np.copy(temp) / s
        if "dpara"  in coefs: out["dpara"]  = np.copy(temp) * m**2/s**3
        if "dperp"  in coefs: out["dperp"]  = np.copy(temp) * m**2/s**3
        if "ddpara" in coefs: out["ddpara"] = np.copy(temp) * m/s**2
        if "q"      in coefs: out["q"]      = np.copy(temp) * m/s**2
        if "dq"     in coefs: out["dq"]     = np.copy(temp) / s
        if "mu0"    in coefs: out["mu0"]    = np.copy(temp)
        if "mu1"    in coefs: out["mu1"]    = np.copy(temp)
        if "dmu0"   in coefs: out["dmu0"]   = np.copy(temp)

        fun = _LIBASCOT.libascot_eval_collcoefs
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, ctypes.c_int, PTR_REAL, ctypes.c_double,
                        ctypes.c_double, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim), Neval, r, phi, z, t, Nv, va, ma, qa,
            out["f"], out["dpara"], out["dperp"], out["k"], out["nu"], out["q"],
            out["dq"], out["ddpara"], out["clog"], out["mu0"], out["mu1"],
            out["dmu0"])

        for d in list(out.keys()):
            if d not in coefs:
                del out[d]
            else:
                if not grid: out[d] = out[d].reshape((n_species, Neval))

        out = [out[k] for k in coefs]
        if len(out) == 1: out = out[0]
        return out

    def input_eval_atomicsigma(self, ma, anum, znum, r, phi, z, t, va, ion,
                               reaction):
        """Deprecated.
        """
        warnings.warn("input_eval_atomicsigma is deprecated. Use "
                      "input_eval_atomiccoefs instead.", DeprecationWarning)
        return

    @parseunits(ma="kg", r="m", phi="rad", z="m", t="s", va="m/s")
    def input_eval_atomiccoefs(self, ma, anum, znum, r, phi, z, t, va,
                                reaction):
        """Evaluate atomic reaction rates for a given test particle.

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
        reaction : {"ionization", "recombination", "charge-exchange",
        "beamstopping"}
            Reaction whose cross-section is computed.

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
        reactions = \
            {v: k for k, v in ascot2py.asigma_reac_type__enumvalues.items()}
        Neval = r.size
        Nv    = va.size
        out   = (np.zeros((Neval,Nv), dtype="f8") + np.nan) / unyt.s

        if reaction == "ionization":
            reaction = reactions["sigmav_recomb"]
        elif reaction == "recombination":
            reaction = reactions["sigmav_ioniz"]
        elif reaction == "charge-exchange":
            reaction = reactions["sigmav_CX"]
        elif reaction == "beamstopping":
            reaction = reactions["sigmav_BMS"]
        else:
            raise ValueError("Unknown reaction")

        fun = _LIBASCOT.libascot_eval_ratecoeff
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, ctypes.c_int, PTR_REAL, ctypes.c_int,
                        ctypes.c_int, ctypes.c_double, ctypes.c_int, PTR_REAL]

        fun(ctypes.byref(self._sim), Neval, r, phi, z, t, Nv, va,
            anum, znum, ma, reaction, out)

        return out

    def input_getplasmaspecies(self):
        """Get species present in plasma input (electrons first).

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
        self._requireinit("plasma")
        fun = _LIBASCOT.libascot_plasma_get_n_species
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM]

        out = {}
        out["nspecies"] = fun(ctypes.byref(self._sim))

        out["mass"]   = np.zeros((out["nspecies"],), dtype="f8")*unyt.kg
        out["charge"] = np.zeros((out["nspecies"],), dtype="f8")*unyt.C
        out["anum"]   = np.zeros((out["nspecies"],), dtype="i4")
        out["znum"]   = np.zeros((out["nspecies"],), dtype="i4")

        fun = _LIBASCOT.libascot_plasma_get_species_mass_and_charge
        fun.restype  = ctypes.c_int
        fun.argtypes = [PTR_SIM, PTR_REAL, PTR_REAL, PTR_INT, PTR_INT]
        fun(ctypes.byref(self._sim),
            out["mass"], out["charge"], out["anum"], out["znum"])

        return out["nspecies"], out["mass"], out["charge"], out["anum"],\
            out["znum"]

    @parseunits(rho="1", theta="rad", phi="rad", time="s")
    def input_rhotheta2rz(self, rho, theta, phi, time, maxiter=100, tol=1e-5):
        """Convert (rho, theta, phi) coordinates to (R,z) positions.

        Parameters
        ----------
        rhovals : array_like (n,)
            Normalized poloidal flux coordinates to be converted.
        theta : array_like (n,)
            Poloidal angle coordinates to be converted.
        phi : array_like (n,)
            Toroidal angle coordinates to be converted.
        time : float
            Time slice (same for all).

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
        rho = np.asarray(rho).ravel().astype(dtype="f8")
        Neval = rho.size
        r   = np.nan * np.zeros((Neval,), dtype="f8") * unyt.m
        z   = np.nan * np.zeros((Neval,), dtype="f8") * unyt.m

        if theta.size == 1:
            theta = theta * np.ones(rho.shape).astype(dtype="f8")
        if phi.size == 1:
            phi = phi * np.ones(rho.shape).astype(dtype="f8")

        fun = _LIBASCOT.libascot_B_field_rhotheta2rz
        fun.restype  = None
        fun.argtypes = [PTR_SIM, ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        ctypes.c_double, ctypes.c_int, ctypes.c_double,
                        PTR_REAL, PTR_REAL]
        fun(ctypes.byref(self._sim),
            Neval, rho, theta, phi, time, maxiter, tol, r, z)

        return (r, z)

    @parseunits(rho="1", theta="rad", phi="rad", time="s")
    def input_findpsi0(self, psi1, nphi=None, phimin=None, phimax=None):
        """Find poloidal flux on axis value numerically.

        Before this function is called, the magnetic field data should contain
        initial guess for the position of the magnetic axis. The algorithm then
        uses the gradient descent method to find psi0. The interpolation is done
        using Ascot's magnetic field interpolation and a little bit of padding
        is added to psi0 so the value can be used as an input parameter for
        the magnetic field without any errors.

        Parameters
        ----------
        psi1 : float
            Poloidal flux at the separatrix.

            This value is used to deduce whether the algorithm searches minimum
            or maximum value when finding psi0.
        nphi : int, optional
            Number of B field grid points in the phi direction between phimin
            and phimax including the end points of this interval.

            Needed for 3D fields.
            For clarity: if the first and the last grid point in the interval
            are the same point, this point is counted twice.
            Example: you have three grid points at 0 deg, 120 deg and 240 deg
            and you use phimin=0, phimax=2pi. The input argument nphi is then
            four instead of three because at 2pi you count again the first grid
            point.
        phimin : float, optional
            Minimum of the phi interval.

            Needed for 3D fields.
        phimax : float, optional
            Maximum of the phi interval.

            Needed for 3D fields.

        Returns
        -------
        r : float
            Axis R-coordinate.
        z : float
            Axis z-coordinate.
        psi0 : float
            Poloidal flux on axis.

        Raises
        ------
        AssertionError
            If required data has not been initialized.
        RuntimeError
            If evaluation in libascot.so failed.
        """
        self._requireinit("bfield")
        tol  = 1e-8
        step = 1e-3
        maxiter = 10**6
        if nphi is None and phimin is None and phimax is None:
            #2D case
            ax = self._eval_bfield(
                1.0*unyt.m, 0.0*unyt.rad, 0.0*unyt.m, 0.0*unyt.s, evalaxis=True)
            psi0 = self._eval_bfield(
                ax["axisr"], 0.0*unyt.rad, ax["axisz"], 0.0*unyt.s, evalrho=True)
            ascent  = int(psi1 < psi0["psi"])

            psi = np.nan * np.zeros((1,), dtype="f8") * unyt.Wb
            rz  = np.zeros((2,), dtype="f8") * unyt.m
            rz[0] = ax["axisr"]
            rz[1] = ax["axisz"]

            fun = _LIBASCOT.libascot_B_field_gradient_descent
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_REAL, PTR_REAL, ctypes.c_double,
                            ctypes.c_double, ctypes.c_int, ctypes.c_int]
            fun(ctypes.byref(self._sim), psi, rz, step, tol, maxiter, ascent)

            if np.isnan(psi[0]):
                raise RuntimeError("Failed to converge.")

            return (rz[0], rz[1], psi)
        elif nphi is not None and phimin is not None and phimax is not None:
            #3D case

            # Divide the 3D field into sectors (phi slices) and find the minimum
            # inside each sector.
            sectoredges = np.linspace(phimin,phimax,nphi)
            nsector = len(sectoredges)-1                      #number of sectors
            psi = np.nan * np.zeros((1,), dtype="f8") * unyt.Wb
            rzphi  = np.zeros((3,), dtype="f8")
            psiconverged = np.nan * np.zeros((1,), dtype="f8") * unyt.Wb

            fun = _LIBASCOT.libascot_B_field_gradient_descent_3d
            fun.restype  = None
            fun.argtypes = [PTR_SIM, PTR_REAL, PTR_REAL, ctypes.c_double,
                            ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.c_int, ctypes.c_int]

            for i in range(nsector):
                phiminsector = sectoredges[i] * unyt.radian
                phimaxsector = sectoredges[i+1] * unyt.radian
                phi = 0.5*(phiminsector + phimaxsector)
                ax = self._eval_bfield(
                    1.0*unyt.m, phi, 0.0*unyt.m, 0.0*unyt.s,
                    evalaxis=True)
                psi0 = self._eval_bfield(
                    ax["axisr"], phi, ax["axisz"], 0.0*unyt.s,evalrho=True)
                ascent  = int(psi1 < psi0["psi"])

                rzphi[0] = ax["axisr"]
                rzphi[1] = ax["axisz"]
                rzphi[2] = phi        #using sector average phi as initial guess

                fun(ctypes.byref(self._sim),
                    psi, rzphi, phiminsector, phimaxsector, step, tol,
                    maxiter, ascent)

                if np.isnan(psi[0]):
                    raise RuntimeError("Failed to converge.")

                if i==0:
                    #Best solution thus far
                    psiconverged[0]=psi[0]
                    rzphiconverged = rzphi
                elif (psi[0] < psiconverged[0]) and ascent == 0:
                    #Hold up, fount something better
                    psiconverged[0] = psi[0]
                    rzphiconverged[:] = rzphi[:]
                elif (psi[0] > psiconverged[0]) and ascent == 1:
                    #Hold up, fount something better
                    psiconverged[0] = psi[0]
                    rzphiconverged[:] = rzphi[:]

                psi *= np.nan    #reset for next iteration

            # Note: rzphiconverged[2] (phi) is not returned!
            return (rzphiconverged[0]*unyt.m, rzphiconverged[1]*unyt.m, psiconverged)
        else:
            #Missing inputs for 3D or unnecessary inputs for 2D
            raise ValueError("All arguments (nphi, phimin, phimax) are needed "
                             "for 3D fields. For 2D fields, none of these "
                             "should be provided.\n\nYou did an oopsie. Search "
                             "your feelings. You know it to be true.")
