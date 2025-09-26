"""This module contains all wall model input variants.

.. autosummary::
    :nosignatures:

    ~Wall2D
    ~CreateWallMixin.create_wall2d

.. rubric:: Classes

.. autoclass:: Wall2D
    :members:

.. autoclass:: CreateWallMixin
    :members:
    :inherited-members:
"""
import ctypes

from a5py.libascot import input_category

from . import axisymmetric
from . import triangular
from .axisymmetric import Wall2D, CreateWall2DMixin
from .triangular import Wall3D, CreateWall3DMixin


@input_category
class Wall(ctypes.Structure):
    """Wrapper for the wall data in libascot.so."""

    _fields_ = [
        ("w2d", ctypes.POINTER(axisymmetric.Struct)),
        ("w3d", ctypes.POINTER(triangular.Struct)),
        ("type", ctypes.c_int32),
        ]


# pylint: disable=too-many-ancestors
class CreateWallMixin(
    CreateWall2DMixin, #CreateWall3DMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create wall input.

    This class just combines all the wall mixin classes.
    """

__all__  = [
    "CreateWallMixin", "Wall2D", "Wall3D",
    ]
