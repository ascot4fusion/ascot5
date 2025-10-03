"""This module contains all wall model input variants.

.. autosummary::
    :nosignatures:

    ~WallContour2D
    ~WallTriangular3D
    ~CreateWallMixin.create_wallcontour2d
    ~CreateWallMixin.create_walltriangular3d

.. rubric:: Classes

.. autoclass:: WallContour2D
    :members:

.. autoclass:: WallTriangular3D
    :members:

.. autoclass:: CreateWallMixin
    :members:
    :inherited-members:
"""
import ctypes

from a5py.libascot import input_category

from . import axisymmetric
from . import triangular
from .axisymmetric import WallContour2D
from .triangular import WallTriangular3D


# pylint: disable=too-few-public-methods
@input_category
class Wall(ctypes.Structure):
    """Wrapper for the wall data in libascot.so."""

    _fields_ = [
        ("contour2d", ctypes.POINTER(axisymmetric.Struct)),
        ("triangular3d", ctypes.POINTER(triangular.Struct)),
        ("type", ctypes.c_int32),
        ]


# pylint: disable=too-many-ancestors
class CreateWallMixin(
    axisymmetric.CreateMixin,
    triangular.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create wall input.

    This class just combines all the wall mixin classes.
    """

__all__  = [
    "CreateWallMixin",
    "WallContour2D",
    "WallTriangular3D",
    ]
