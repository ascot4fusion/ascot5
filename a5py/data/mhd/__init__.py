"""This module contains all MHD eigenmode input variants and combines their
factory methods to a single class."""
import ctypes

from a5py.libascot import input_category

from . import dynamic
from . import stationary
from .dynamic import MhdDynamic
from .stationary import MhdStationary


@input_category
class Mhd(ctypes.Structure):
    """Wrapper for the MHD data in libascot.so."""

    _fields_ = [
        ("stationary", ctypes.POINTER(dynamic.Struct)),
        ("dynamic", ctypes.POINTER(stationary.Struct)),
        ("type", ctypes.c_int32),
        ]


# pylint: disable=too-many-ancestors
class CreateMhdMixin(
    dynamic.CreateMixin,
    stationary.CreateMixin,
    ):
    """Mixin class used by `Data` to create MHD eigenmode input.

    This class just combines all the MHD mixin classes.
    """

__all__  = [
    "CreateMhdMixin", "MhdStationary", "MhdDynamic",
    ]
