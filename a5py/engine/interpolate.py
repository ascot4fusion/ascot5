"""Routines that interpolate or rely on interpolating the input data.

Any algorithms that are not related to the simulation itself, but are still
directly linked to libascot, should be defined here.
"""
import ctypes

import unyt
import numpy as np

from a5py.libascot import LIBASCOT
from a5py import physlib

from a5py.physlib.formulas import resolve_quantities, compute_quantities
from a5py.data.bfield import Bfield
from a5py.data.efield import Efield
from a5py.data.plasma import Plasma
from a5py.data.neutral import Neutral
from a5py.data.boozer import Struct as BoozerMap
from a5py.data.mhd import Mhd
from .functions import init_fun, PTR_DOUBLE, PTR_INT


init_fun(
    "ascot_interpolate",
    ctypes.POINTER(Bfield),
    ctypes.POINTER(Efield),
    ctypes.POINTER(Plasma),
    ctypes.POINTER(Neutral),
    ctypes.POINTER(BoozerMap),
    ctypes.POINTER(Mhd),
    ctypes.c_size_t,
    ctypes.c_int32,
    *([PTR_DOUBLE]*20)
    )

init_fun(
    "ascot_map_rhotheta_to_rz",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.c_double,
    ctypes.c_double,
    *([PTR_DOUBLE]*5)
    )

init_fun(
    "ascot_find_psi_on_axis_2d",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.c_double,
    ctypes.c_double,
    *([PTR_DOUBLE]*2)
    )

init_fun(
    "ascot_find_psi_on_axis_3d",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    *([PTR_DOUBLE]*2)
    )


qnt_mapping = {
    "b": ("br", "bphi", "bz"),
    "bjac": ("brdr", "brdphi", "brdz", "bphidr", "bphidphi", "bphidz",
             "bzdr", "bzdphi", "bzdz"),
    "psi": ("psi", "psidr", "psidphi", "psidz"),
    "rho": ("rho", "rhodpsi"),
    "e": ("er", "ez", "ephi"),
    "n": ("ne", "ni"),
    "t": ("te", "ti"),
    "n0": ("n0",),
    "t0": ("t0",),
    "theta": ("theta", "dthetadr", "dthetadphi", "dthetadz"),
    "zeta": ("zeta", "dzetadr", "dzetadphi", "dzetadz"),
    "alpha": ("alpha", "alphadr", "alphadphi", "alphadz", "alphadt"),
    "Phi": ("Phi", "dPhidr", "dPhidphi", "dPhidz", "dPhidt"),
    "mhd_b": ("br_mhd", "bphi_mhd", "bz_mhd"),
    "mhd_e": ("er_mhd", "ephi_mhd", "ez_mhd"),
    "mhd_phi": ("phi", "dphidr", "dphidphi", "dphidz"),
}


@physlib.parse_units(r="m", phi="deg", z="m", t="s")
def evaluate(
    r,
    phi,
    z,
    t,
    *qnt,
    grid=False,
    bfield=None,
    efield=None,
    plasma=None,
    neutral=None,
    boozer=None,
    mhd=None,
    atomic=None,
    nmode=0,
    ):
    """Evaluate input quantities at given coordinates.

    This method uses ASCOT5 C-routines for evaluating inputs exactly as
    they are evaluated during a simulation.

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
    """
    r, phi, z, t = process_query_points(r, phi, z, t, grid)
    available = []
    for q in qnt_mapping.values():
        available.extend(q)
    available.extend(["r"])
    fundamentals_needed, all_qnt = resolve_quantities(available, qnt)
    if "r" in fundamentals_needed:
        fundamentals_needed.remove("r")
    fundamentals = interpolate(
        r, phi, z, t, nmode, *fundamentals_needed, bfield=bfield, efield=efield,
        plasma=plasma, neutral=neutral, boozer=boozer, mhd=mhd, atomic=atomic
        )
    fundamentals["r"] = r
    evaluated = compute_quantities(fundamentals, all_qnt)
    results = [evaluated[q] for q in qnt]
    for i in range(len(results)):
        try:
            results[i] = results[i].to("A/m**2")
        except:
            pass

    if len(results) == 1:
        return results[0]
    return results


def init_return_array(shape, dtype="f8", units=None):
    """Initialize an array with NaN values and given units.

    This function is a shortcut to initialize suitable arrays to be passed to
    libascot.

    The library itself doesn't care about units, but the units don't
    interfere with it so they can be added here already.

    NaN values are used in initialization since libascot doesn't produce an
    error if data is interpolated outside the grid; it just doesn't return a
    value there. This way the array will have NaNs where interpolation failed
    and actual values otherwise.

    Parameters
    ----------
    shape : tuple of int
        Shape of the array to be initialized.
    dtype : str or np.dtype, optional
        Data type of the array.
    units : str or Unit, optional
        Units of the array.

    Returns
    -------
    arr : ndarray or unyt_array
        Initialized array with NaN values and optional units.
    """
    arr = np.full(shape, np.nan, dtype=dtype)
    if units is None:
        return arr
    return unyt.unyt_array(arr, units)


def process_query_points(r, phi, z, t, grid):
    r"""Process query points so that they can be passed to libascot.so.

    This routine ensures the coordinates are double precision and in correct
    units. If the coordinates form a grid, then return the grid point
    coordinates. If not, then check that the coordinates have same size.

    Parameters
    ----------
    r, phi, z, t : ndarray
        Coordinates.
    grid : bool
        Flag indicating whether the coordinates form a grid.

    Returns
    -------
    r, phi, z, t : ndarray
        Processed coordinates.
    """
    def ensure_units_and_dtype(x, units):
        return x.to(units).astype(dtype="f8")
    r, phi, z, t = (
        ensure_units_and_dtype(r, "m"),
        ensure_units_and_dtype(phi, "rad"),
        ensure_units_and_dtype(z, "m"),
        ensure_units_and_dtype(t, "s")
        )

    if grid:
        if any(arr.squeeze().ndim > 1 for arr in (r, phi, z, t)):
            raise ValueError("Grid inputs must be 1D arrays")
        R, PHI, Z, T = np.meshgrid(
            r.ravel(), phi.ravel(), z.ravel(), t.ravel(), indexing="ij"
            )
        # Meshgrid doesn't ensure that arrays are C_CONTIGUOUS, and in fact Z
        # isn't, so we make a C_CONTIGUOUS copy.
        Z = np.ascontiguousarray(Z)
        return R, PHI, Z, T

    shapes = [arr.shape for arr in (r, phi, z, t)]
    nonsingletons = [s for s in shapes if s != ()]
    if len(nonsingletons) > 1 and len(set(nonsingletons)) != 1:
        raise ValueError(
            f"Input arrays (r, phi, z, t) must have same shape or be singletons. "
            f"Got {r.shape}, {phi.shape}, {z.shape}, {t.shape}."
        )

    arrshape = () if not len(nonsingletons) else nonsingletons[0]
    def broadcast(arr):
        if arr.size > 1:
            return arr
        return np.full(arrshape, arr.ravel()[0]) * arr.units
    return tuple(broadcast(arr) for arr in (r, phi, z, t))


def interpolate(r, phi, z, t, nmode, *qnt, **input_variants):
    """Evaluates fundamental quantities at given coordinates.

    Parameters
    ----------
    r, phi, z, t : ndarray
        Coordinates.
    nmode : int
        Mode number.
    *qnt : str
        Quantities to evaluate (e.g., "br", "psi", "rho", ...).
    **input_variants : dict
        Input variants to be used in interpolation (bfield, efield, plasma,
        neutral, boozer, mhd).

    Returns
    -------
    dict[str, ndarray]
        Mapping of requested quantities to evaluated arrays with units.
    """
    unit_mapping = {
        "T": ("br", "bz", "bphi", "mhd_br", "mhd_bz", "mhd_bphi", "bphidphi",
              "brdphi", "bzdphi",),
        "T/m": ("brdr", "brdz", "bphidr", "bphidz", "bzdr", "bzdz",),
        "V/m": ("er", "ez", "ephi", "mhd_er", "mhd_ez", "mhd_ephi",),
        "Wb/rad": ("psi", "dpsidphi",),
        "(Wb/rad)/m": ("dpsidr", "dpsidz",),
        "1": ("rho",),
        "1/(Wb/rad)": ("drhodpsi",),
        "m**-3": ("ne", "ni", "n0",),
        "eV": ("te", "ti", "t0",),
    }


    def init_return_arrays(qnt, query_shape, input_variants):
        """Initialize return arrays.

        Use NULL pointer for quantities that are not requested to avoid
        unnecessary memory allocation. Some quantities are always evaluated
        simultaneously (e.g. magnetic field components) so store them in the
        same array for simpler handling.

        The array shape is always (n_simultaneous, n_query_points) on the C
        side, but we assume that the query arrays have the correct ordering so
        that we can preserve their shape, i.e. the arrays are
        (n_simultaneous, *query_shape). The special cases are plasma and neutral
        evaluations where the array shape depends on the number of species.

        This function also checks that all required inputs are present and
        raises exception otherwise.
        """
        required_map = {
            ("b", "bjac", "psi", "rho"): ["bfield"],
            ("e",): ["bfield", "efield"],
            ("n", "t"): ["bfield", "plasma"],
            ("n0", "t0"): ["bfield", "neutral"],
            ("theta", "zeta"): ["bfield", "boozer"],
            ("alpha", "Phi", "mhd_b", "mhd_e", "mhd_phi"): ["bfield", "boozer", "mhd"],
        }
        return_arrays, required, nullpointer = {}, [], None
        for name, subqnt in qnt_mapping.items():
            qnt_not_requested = not any(q in subqnt for q in qnt)
            if qnt_not_requested:
                return_arrays[name] = nullpointer
                continue

            for keys, deps in required_map.items():
                if name in keys and not deps in required:
                    required.extend(deps)
            for req in required:
                if req not in input_variants or input_variants[req] is None:
                    raise ValueError(
                        f"Input variant {req} is required for this evaluation."
                    )

            if name == "n":
                ns = len(input_variants["plasma"].species)
                shape = (ns + 1, *query_shape)
            elif name in ["n0", "t0"]:
                ns = len(input_variants["neutral"].species)
                shape = (ns, *query_shape)
            else:
                shape = (len(subqnt), *query_shape)
            return_arrays[name] = init_return_array(shape, dtype="f8")
        return return_arrays, required


    def init_inputs(input_variants, required_inputs):
        """Initialize input data structs (only those that are required)."""
        inputs = {}
        for inp in ["bfield", "efield", "plasma", "neutral", "boozer", "mhd"]:
            val = input_variants.get(inp)
            if val is not None and inp in required_inputs:
                cls = globals()[inp.capitalize()]
                interface = cls()
                interface.use(val)
                inputs[inp] = interface
            else:
                inputs[inp] = None
        return inputs

    def set_units(result):
        for k, v in result.items():
            for unit, quantity in unit_mapping.items():
                if k in quantity:
                    result[k] = unyt.unyt_array(v, unit)
                    break

    query_shape, query_size = r.shape, r.size
    return_arrays, required_inputs = (
        init_return_arrays(qnt, query_shape, input_variants)
        )
    inputs = init_inputs(input_variants, required_inputs)
    ptrs = {
        k: ctypes.byref(v) if v is not None else None for k, v in inputs.items()
        }
    LIBASCOT.ascot_interpolate(
        ptrs["bfield"], ptrs["efield"], ptrs["plasma"], ptrs["neutral"],
        ptrs["boozer"], ptrs["mhd"], query_size, nmode, r, phi, z, t,
        *return_arrays.values()
        )
    result = {}
    for qname in qnt:
        for name, mapping in qnt_mapping.items():
            if qname in mapping and name in return_arrays:
                idx = mapping.index(qname)
                result[qname] = return_arrays[name][idx]

    set_units(result)
    return result
