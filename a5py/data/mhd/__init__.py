"""This module contains all MHD eigenmode input variants and combines their
factory methods to a single class."""
import ctypes

from .stationary import MhdStationary, CreateMhdStationaryMixin
from .dynamic import MhdDynamic, CreateMhdDynamicMixin


class Mhd(ctypes.Structure):
    """Wrapper for the MHD data in libascot.so."""

    #_pack_ = 1
    _fields_ = [
        ('stat', ctypes.POINTER(MhdStationary.Struct)),
        ('nonstat', ctypes.POINTER(MhdDynamic.Struct)),
        ('type', ctypes.c_uint32),
        ]

    def use(self, variant):
        """Initialize the pointer and set the type corresponding to the data."""
        variant_vs_name_and_enum = {
            "MhdStationary": ("stat", 0),
            "MhdDynamic": ("nonstat", 1),
            }
        variant.stage()
        name, enum = variant_vs_name_and_enum[variant.variant]
        setattr(self, name, ctypes.pointer(variant._struct_))
        self.type = ctypes.c_uint32(enum)


# pylint: disable=too-many-ancestors
class CreateMhdMixin(
    CreateMhdStationaryMixin, CreateMhdDynamicMixin,
    ):
    """Mixin class used by `Data` to create MHD eigenmode input.

    This class just combines all the MHD mixin classes.
    """

__all__  = [
    "CreateMhdMixin", "MhdStationary", "MhdDynamic",
    ]
