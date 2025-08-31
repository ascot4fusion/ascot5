"""Routines that interpolate the input data."""
import ctypes

import unyt
import numpy as np

from ..libascot import LIBASCOT
from .. import physlib


def _init_return_arrays(shape, dtype="f8", units=None):
    arr = np.zeros(shape, dtype=dtype) + np.nan
    if units is None:
        return arr
    return arr * unyt.unyt_quantity.from_string(units)

def _parse_coordinates(r, phi, z, t, grid=False):

    r = r.astype(dtype="f8").ravel()
    z = z.astype(dtype="f8").ravel()
    t = t.astype(dtype="f8").ravel()
    phi = phi.astype(dtype="f8").ravel().to("rad")

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
    r, phi, z, t = _parse_coordinates()

    v = {}

    bnorm = lambda bx, by, bz: np.sqrt(bx**2 + by**2 + bz**2)

    # 2. Extract argument names
    def get_args(func):
        return inspect.signature(func).parameters.keys()

    # 3. Evaluation using a values dictionary
    def evaluate(func, values, key=None):
        args = get_args(func)
        arg_values = [values[arg] for arg in args]
        result = func(*arg_values)
        if key is None:
            key = func.__name__ if hasattr(func, '__name__') else 'result'
        values[key] = result
        return result

    def resolve_dependent_quantities(requested):
        """Find which quantities are needed to compute the requested quantities.
        """

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


        resolved = set()
        to_process = list(requested)
        while to_process:
            q = to_process.pop()
            if q in resolved:
                continue
            required = deps_map.get(q, [])
            to_process.extend(required)
            resolved.add(q)
        return resolved

    def isin(*available):
        return any(q in qnt for q in available)


    if isin("rho", "psi", "rhodpsi", "psidr", "psidphi", "psidz"):
        v.update(**eval_bfield(r, phi, z, t, evalrho=True))
    if isin("br", "bphi", "bz", "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
            "bphidz", "bzdr", "bzdphi", "bzdz", "divb", "bnorm", "jnorm", "jr",
            "jphi", "jz", "gradbr", "gradbphi", "gradbz", "curlbr", "curlbphi",
            "curlbz"):
        v.update(**eval_bfield(r, phi, z, t, evalb=True))
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
    if any(q in qnt for q in ["axisr", "axisz", "rminor"]):
        out.update(self._eval_bfield(r, phi, z, t, evalaxis=True))
        out["rminor"] = np.sqrt(  ( out["axisr"] - r )**2
                                + ( out["axisz"] - z )**2 )
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


def eval_bfield(bfield, r, phi, z, t, evalb=False, evalrho=False, evalaxis=False):
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
    n, shape = r.size, r.shape

    v = {}
    if evalb:
        v.update(
            **_init_return_arrays(
                ("br", "bphi", "bz", "brdr", "brdphi", "brdz", "bphidr",
                 "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz"),
                [shape]*12,
                ["f8"]*12,
                ["T", "T", "T", "T/m", "T", "T/m", "T/m", "T", "T/m", "T/m",
                 "T", "T/m"],
            )
        )
        LIBASCOT.libascot_B_field_eval_B_dB(
            ctypes.byref(bfield), n, r, phi, z, t, *list(v.values)
            )
    return v