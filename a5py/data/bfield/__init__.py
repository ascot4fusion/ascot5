"""This module contains all magnetic field input variants.

.. autosummary::
    :nosignatures:

    ~BfieldCartesian
    ~BfieldAnalytical
    ~BfieldSpline2D
    ~BfieldSpline3D
    ~BfieldStellarator
    ~CreateBfieldMixin.create_bfieldcartesian
    ~CreateBfieldMixin.create_bfieldanalytical
    ~CreateBfieldMixin.create_bfieldspline2d
    ~CreateBfieldMixin.create_bfieldspline3d
    ~CreateBfieldMixin.create_bfieldstellarator

.. rubric:: Classes

.. autoclass:: BfieldCartesian
    :members:

.. autoclass:: BfieldAnalytical
    :members:

.. autoclass:: BfieldSpline2D
    :members:

.. autoclass:: BfieldSpline3D
    :members:

.. autoclass:: BfieldStellarator
    :members:

.. autoclass:: CreateBfieldMixin
    :members:
    :inherited-members:

"""
import ctypes

from a5py.libascot import input_category

from . import cartesian
from . import analytical
from . import spline2d
from . import spline3d
from . import stellarator
from .cartesian import BfieldCartesian
from .analytical import BfieldAnalytical
from .spline2d import BfieldSpline2D
from .spline3d import BfieldSpline3D
from .stellarator import BfieldStellarator


# pylint: disable=too-few-public-methods
@input_category
class Bfield(ctypes.Structure):
    """Wrapper for the magnetic field data in B_field.h."""

    _fields_ = [
        ("cartesian", ctypes.POINTER(cartesian.Struct)),
        ("analytical", ctypes.POINTER(analytical.Struct)),
        ("spline2d", ctypes.POINTER(spline2d.Struct)),
        ("spline3d", ctypes.POINTER(spline3d.Struct)),
        ("stellarator", ctypes.POINTER(stellarator.Struct)),
        ("type", ctypes.c_int32),
    ]


# pylint: disable=too-many-ancestors
class CreateBfieldMixin(
    cartesian.CreateMixin,
    analytical.CreateMixin,
    spline2d.CreateMixin,
    spline3d.CreateMixin,
    stellarator.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create magnetic field input.

    This class just combines all the magnetic field mixin classes.
    """

__all__  = [
    "CreateBfieldMixin",
    "BfieldCartesian",
    "BfieldAnalytical",
    "BfieldSpline2D",
    "BfieldSpline3D",
    "BfieldStellarator",
    ]
