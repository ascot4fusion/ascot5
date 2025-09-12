"""This module contains all plasma input variants and combines their
factory methods to a single class."""
import ctypes

from .radial import Plasma1D, CreatePlasma1DMixin
from .radialdynamic import Plasma1DDynamic, CreatePlasma1DDynamicMixin


class Plasma(ctypes.Structure):
    """Wrapper for the plasma data in libascot.so."""

    #_pack_ = 1
    _fields_ = [
        ('plasma_1D', ctypes.POINTER(Plasma1D.Struct)),
        ('plasma_1Dt', ctypes.POINTER(Plasma1DDynamic.Struct)),
        ('plasma_1DS', ctypes.POINTER(None)),
        ('type', ctypes.c_uint32),
        ]

    def use(self, variant):
        """Initialize the pointer and set the type corresponding to the data."""
        variant_vs_name_and_enum = {
            "Plasma1D": ("plasma_1D", 0),
            "Plasma1DDynamic": ("plasma_1Dt", 1),
            }
        variant.stage()
        name, enum = variant_vs_name_and_enum[variant.variant]
        setattr(self, name, ctypes.pointer(variant._struct_))
        self.type = ctypes.c_uint32(enum)


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
