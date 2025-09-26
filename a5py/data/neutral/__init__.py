"""This module contains all neutral density input variants and combines their
factory methods to a single class."""
import ctypes

from a5py.libascot import input_category

from . import radial
from . import arbitrary
from .radial import Neutral1D
from .arbitrary import Neutral3D


@input_category
class Neutral(ctypes.Structure):
    """Wrapper for the neutral data in libascot.so."""

    _fields_ = [
        ('N01D', ctypes.POINTER(radial.Struct)),
        ('N03D', ctypes.POINTER(arbitrary.Struct)),
        ('type', ctypes.c_uint32),
        ]


# pylint: disable=too-many-ancestors
class CreateNeutralMixin(
    radial.CreateMixin,
    arbitrary.CreateMixin,
    ):
    """Mixin class used by `Data` to create neutral density input.

    This class just combines all the neutral density mixin classes.
    """

__all__  = [
    "CreateNeutralMixin", "Neutral1D", "Neutral3D",
    ]
