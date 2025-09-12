"""This module contains all electric field input variants and combines their
factory methods to a single class."""
import ctypes

from .cartesian import EfieldCartesian, CreateEfieldCartesianMixin
from .radialpotential import (
    EfieldRadialPotential, CreateEfieldRadialPotentialMixin,
    )


class Efield(ctypes.Structure):
    """Wrapper for the electric field data in libascot.so."""

    #_pack_ = 1
    _fields_ = [
        ('ETC', ctypes.POINTER(EfieldCartesian.Struct)),
        ('E1DS', ctypes.POINTER(EfieldRadialPotential.Struct)),
        ('type', ctypes.c_uint32),
    ]

    def use(self, variant):
        """Initialize the pointer and set the type corresponding to the data."""
        variant_vs_name_and_enum = {
            "EfieldCartesian": ("ETC", 0),
            "EfieldRadialPotential": ("E1DS", 1),
            }
        variant.stage()
        name, enum = variant_vs_name_and_enum[variant.variant]
        setattr(self, name, ctypes.pointer(variant._struct_))
        self.type = ctypes.c_uint32(enum)


# pylint: disable=too-many-ancestors
class CreateEfieldMixin(
    CreateEfieldCartesianMixin,
    CreateEfieldRadialPotentialMixin,
    ):
    """Mixin class used by `Data` to create electric field input.

    This class just combines all the electric field mixin classes.
    """

__all__  = [
    "CreateEfieldMixin",
    "EfieldCartesian",
    "EfieldRadialPotential",
    ]
