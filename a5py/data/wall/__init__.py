"""This module contains all wall input variants and combines their factory
methods to a single class."""
import ctypes

from .axisymmetric import Wall2D, CreateWall2DMixin
from .triangular import Wall3D, CreateWall3DMixin


class Wall(ctypes.Structure):
    """Wrapper for the wall data in libascot.so."""

    _pack_ = 1
    _fields_ = [
        ('w2d', ctypes.POINTER(Wall2D.Struct)),
        ('w3d', ctypes.POINTER(Wall3D.Struct)),
        ('type', ctypes.c_int32),
        ]


# pylint: disable=too-many-ancestors
class CreateWallMixin(
    CreateWall2DMixin, CreateWall3DMixin,
    ):
    """Mixin class used by `Data` to create wall input.

    This class just combines all the wall input mixin classes.
    """

__all__  = [
    "CreateWallMixin", "Wall2D", "Wall3D",
    ]
