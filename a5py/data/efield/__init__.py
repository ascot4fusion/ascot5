"""This module contains all electric field input variants and combines their
factory methods to a single class."""
import ctypes

from .cartesian import EfieldCartesian, CreateEfieldCartesianMixin
from .radialpotential import (
    EfieldRadialPotential, CreateEfieldRadialPotentialMixin,
    )


class Efield(ctypes.Structure):
    """Wrapper for the electric field data in libascot.so."""

    _pack_ = 1
    _fields_ = [
        ('ETC', ctypes.POINTER(EfieldCartesian.Struct)),
        ('E1DS', ctypes.POINTER(EfieldRadialPotential.Struct)),
        ('type', ctypes.c_uint32),
    ]


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
