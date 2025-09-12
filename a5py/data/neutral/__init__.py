"""This module contains all neutral density input variants and combines their
factory methods to a single class."""
import ctypes

from .radial import Neutral1D, CreateNeutral1DMixin
from .arbitrary import Neutral3D, CreateNeutral3DMixin


class Neutral(ctypes.Structure):
    """Wrapper for the neutral data in libascot.so."""

    #_pack_ = 1
    _fields_ = [
        ('N01D', ctypes.POINTER(Neutral1D.Struct)),
        ('N03D', ctypes.POINTER(Neutral3D.Struct)),
        ('type', ctypes.c_uint32),
        ]

    def use(self, variant):
        """Initialize the pointer and set the type corresponding to the data."""
        variant_vs_name_and_enum = {
            "Neutral1D": ("N01D", 0),
            "Neutral3D": ("N03D", 1),
            }
        variant.stage()
        name, enum = variant_vs_name_and_enum[variant.variant]
        setattr(self, name, ctypes.pointer(variant._struct_))
        self.type = ctypes.c_uint32(enum)


# pylint: disable=too-many-ancestors
class CreateNeutralMixin(
    CreateNeutral1DMixin, CreateNeutral3DMixin,
    ):
    """Mixin class used by `Data` to create neutral density input.

    This class just combines all the neutral density mixin classes.
    """

__all__  = [
    "CreateNeutralMixin", "Neutral1D", "Neutral3D",
    ]
