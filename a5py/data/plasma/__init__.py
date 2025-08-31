"""This module contains all plasma input variants and combines their
factory methods to a single class."""
import ctypes

from .radial import Plasma1D, CreatePlasma1DMixin
from .radialdynamic import Plasma1DDynamic, CreatePlasma1DDynamicMixin


class Plasma(ctypes.Structure):
    """Wrapper for the plasma data in libascot.so."""

    _pack_ = 1
    _fields_ = [
        ('plasma_1D', ctypes.POINTER(Plasma1D.Struct)),
        ('plasma_1Dt', ctypes.POINTER(Plasma1DDynamic.Struct)),
        ('plasma_1DS', ctypes.POINTER(None)),
        ('type', ctypes.c_uint32),
        ]


# pylint: disable=too-many-ancestors
class CreatePlasmaMixin(
    CreatePlasma1DMixin, CreatePlasma1DDynamicMixin,
    ):
    """Mixin class used by `Data` to create plasma input.

    This class just combines all the plasma input mixin classes.
    """

__all__  = [
    "CreatePlasmaMixin", "Plasma1D", "Plasma1DDynamic",
    ]
