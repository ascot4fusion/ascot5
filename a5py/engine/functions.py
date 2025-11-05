"""Defines all wrapped functions except those used to initialize data.
"""
import ctypes
from numpy.ctypeslib import ndpointer

from ..libascot import LIBASCOT, init_fun


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


PTR_DOUBLE = _ndpointerornull(ctypes.c_double, flags="C_CONTIGUOUS")
"""Double valued numpy array or None (NULL)."""


PTR_INT = _ndpointerornull(ctypes.c_int32, flags="C_CONTIGUOUS")
"""Int valued numpy array or None (NULL)."""


init_fun("get_endcond", ctypes.c_char_p, restype=ctypes.c_int)

def get_endcond(name: str) -> int:
    return LIBASCOT.get_endcond(name.encode("utf-8"))

def get_enum_value(enum):
    return ctypes.c_uint.in_dll(LIBASCOT, enum).value


END_CONDITIONS = None
