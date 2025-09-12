"""This module contains all wall input variants and combines their factory
methods to a single class."""
import ctypes

from .axisymmetric import Wall2D, CreateWall2DMixin
from .triangular import Wall3D, CreateWall3DMixin


class Wall(ctypes.Structure):
    """Wrapper for the wall data in libascot.so."""

    #_pack_ = 1
    _fields_ = [
        ('w2d', ctypes.POINTER(Wall2D.Struct)),
        ('w3d', ctypes.POINTER(Wall3D.Struct)),
        ('type', ctypes.c_int32),
        ]

    def use(self, variant):
        """Initialize the pointer and set the type corresponding to the data."""
        variant_vs_name_and_enum = {
            "Wall2D": ("w2d", 0),
            "Wall3D": ("w3d", 1),
            }
        variant.stage()
        name, enum = variant_vs_name_and_enum[variant.variant]
        setattr(self, name, ctypes.pointer(variant._struct_))
        self.type = ctypes.c_uint32(enum)


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
