"""Defines all wrapped functions except those used to initialize data.
"""
import ctypes
from numpy.ctypeslib import ndpointer

from ..data.bfield import Bfield
from ..data.efield import Efield
from ..data.plasma import Plasma
from ..data.neutral import Neutral
from ..data.boozer import BoozerMap
from ..data.mhd import Mhd
from ..libascot import LIBASCOT


def _ndpointerornull(*args, **kwargs):
    """Produce ndpointer which can also be a NULL pointer.

    This wrapper is required since otherwise it is not possible to supply None
    (indicating NULL) as an argument that requires a numpy array.
    """
    base = ndpointer(*args, **kwargs)
    def from_param(_, obj):
        if obj is None:
            return obj
        return base.from_param(obj)
    return type(base.__name__, (base,), {"from_param": classmethod(from_param)})


def _init_fun(name, *argtypes, restype=None):
    """Initialize function in libascot.so."""
    fun = getattr(LIBASCOT, name)
    fun.argtypes = argtypes
    fun.restype = restype


PTR_DOUBLE = _ndpointerornull(ctypes.c_double, flags="C_CONTIGUOUS")
"""Double valued numpy array or None (NULL)."""


PTR_INT = _ndpointerornull(ctypes.c_int32, flags="C_CONTIGUOUS")
"""Int valued numpy array or None (NULL)."""


_init_fun(
    "libascot_interpolate",
    ctypes.POINTER(Bfield),
    ctypes.POINTER(Efield),
    ctypes.POINTER(Plasma),
    ctypes.POINTER(Neutral),
    ctypes.POINTER(BoozerMap.Struct),
    ctypes.POINTER(Mhd),
    ctypes.c_int32,
    ctypes.c_int32,
    *([PTR_DOUBLE]*20)
    )

_init_fun(
    "libascot_find_axis",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    *([PTR_DOUBLE]*2)
    )

_init_fun(
    "libascot_rhotheta2rz",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.c_double,
    ctypes.c_double,
    *([PTR_DOUBLE]*5)
    )

_init_fun(
    "libascot_gradient_descent",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.c_double,
    ctypes.c_double,
    *([PTR_DOUBLE]*2)
    )

_init_fun(
    "libascot_gradient_descent_3d",
    ctypes.POINTER(Bfield),
    ctypes.c_int32,
    ctypes.c_int32,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    *([PTR_DOUBLE]*2)
    )

def get_enum_value(enum):
    return ctypes.c_uint.in_dll(LIBASCOT, enum).value


#END_CONDITIONS = _load_end_condions()
END_CONDITIONS = {
    "reached_time_limit":get_enum_value("endcond_tlim"),
    "below_min_energy":get_enum_value("endcond_emin"),
    "thermalized":get_enum_value("endcond_therm"),
    "hit_wall":get_enum_value("endcond_wall"),
    "below_rho_limit":get_enum_value("endcond_rhomin"),
    "above_rho_limit":get_enum_value("endcond_rhomax"),
    "completed_poloidal_orbits":get_enum_value("endcond_polmax"),
    "completed_toroidal_orbits":get_enum_value("endcond_tormax"),
    "simulation_not_finished":get_enum_value("endcond_cpumax"),
    "finished_gc_in_hybrid_mode":get_enum_value("endcond_hybrid"),
    "neutralized":get_enum_value("endcond_neutr"),
    "ionized":get_enum_value("endcond_ioniz"),
}
"""Mapping from end condition name to end condition code."""
