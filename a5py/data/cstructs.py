"""Python wrappers for various helper C structures."""
import ctypes

from ..libascot import LIBASCOT


ENUM = ctypes.c_uint32
"""Python wrapper for an enum defined in C."""


class Linear1D(ctypes.Structure):
    """Python wrapper for the 1D linear-interpolator defined in linint.h."""

    #_pack_ = 1
    _fields_ = [
        ('n_x', ctypes.c_int32),
        ('bc_x', ctypes.c_int32),
        ('x_min', ctypes.c_double),
        ('x_max', ctypes.c_double),
        ('x_grid', ctypes.c_double),
        ('c', ctypes.POINTER(ctypes.c_double)),
        ]


class Linear3D(ctypes.Structure):
    """Python wrapper for the 3D linear-interpolator defined in linint.h."""

    #_pack_ = 1
    _fields_ = [
        ('n_x', ctypes.c_int32),
        ('n_y', ctypes.c_int32),
        ('n_z', ctypes.c_int32),
        ('bc_x', ctypes.c_int32),
        ('bc_y', ctypes.c_int32),
        ('bc_z', ctypes.c_int32),
        ('x_min', ctypes.c_double),
        ('x_max', ctypes.c_double),
        ('x_grid', ctypes.c_double),
        ('y_min', ctypes.c_double),
        ('y_max', ctypes.c_double),
        ('y_grid', ctypes.c_double),
        ('z_min', ctypes.c_double),
        ('z_max', ctypes.c_double),
        ('z_grid', ctypes.c_double),
        ('c', ctypes.POINTER(ctypes.c_double)),
        ]


class Spline1D(ctypes.Structure):
    """Python wrapper for the 1D spline-interpolator defined in interp.h."""

    #_pack_ = 1
    _fields_ = [
        ('n_x', ctypes.c_int32),
        ('bc_x', ctypes.c_int32),
        ('x_min', ctypes.c_double),
        ('x_max', ctypes.c_double),
        ('x_grid', ctypes.c_double),
        ('c', ctypes.POINTER(ctypes.c_double)),
        ]


class Spline2D(ctypes.Structure):
    """Python wrapper for the 2D spline-interpolator defined in interp.h."""

    #_pack_ = 1
    _fields_ = [
        ('n_x', ctypes.c_int32),
        ('n_y', ctypes.c_int32),
        ('bc_x', ctypes.c_int32),
        ('bc_y', ctypes.c_int32),
        ('x_min', ctypes.c_double),
        ('x_max', ctypes.c_double),
        ('x_grid', ctypes.c_double),
        ('y_min', ctypes.c_double),
        ('y_max', ctypes.c_double),
        ('y_grid', ctypes.c_double),
        ('c', ctypes.POINTER(ctypes.c_double)),
        ]


class Spline3D(ctypes.Structure):
    """Python wrapper for the 3D spline-interpolator defined in interp.h."""

    #_pack_ = 1
    _fields_ = [
        ('n_x', ctypes.c_int32),
        ('n_y', ctypes.c_int32),
        ('n_z', ctypes.c_int32),
        ('bc_x', ctypes.c_int32),
        ('bc_y', ctypes.c_int32),
        ('bc_z', ctypes.c_int32),
        ('x_min', ctypes.c_double),
        ('x_max', ctypes.c_double),
        ('x_grid', ctypes.c_double),
        ('y_min', ctypes.c_double),
        ('y_max', ctypes.c_double),
        ('y_grid', ctypes.c_double),
        ('z_min', ctypes.c_double),
        ('z_max', ctypes.c_double),
        ('z_grid', ctypes.c_double),
        ('c', ctypes.POINTER(ctypes.c_double)),
        ]


def free(ptr):
    """Deallocate pointer initialized in libascot.so.

    Parameters
    ----------
    ptr : ctypes.POINTER[Any]
        Pointer to deallocate.
    """
    fun = LIBASCOT.libascot_deallocate
    fun.restype = None
    fun.argtypes = [ctypes.POINTER(None)]
    fun(ptr)
