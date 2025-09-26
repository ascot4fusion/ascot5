"""This module contains all magnetic field input variants.

.. autosummary::
    :nosignatures:

    ~BfieldCartesian
    ~BfieldAnalytical
    ~Bfield2D
    ~Bfield3D
    ~BfieldStellarator
    ~CreateBfieldMixin.create_bfieldcartesian
    ~CreateBfieldMixin.create_bfieldanalytical
    ~CreateBfieldMixin.create_bfield2d
    ~CreateBfieldMixin.create_bfield3d
    ~CreateBfieldMixin.create_bfieldstellarator

.. rubric:: Classes

.. autoclass:: BfieldCartesian
    :members:

.. autoclass:: BfieldAnalytical
    :members:

.. autoclass:: Bfield2D
    :members:

.. autoclass:: Bfield3D
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
from . import axisymmetric
from . import perturbed
from . import stellarator
from .cartesian import BfieldCartesian
from .analytical import BfieldAnalytical
from .axisymmetric import Bfield2D
from .perturbed import Bfield3D
from .stellarator import BfieldStellarator


# pylint: disable=too-few-public-methods
@input_category
class Bfield(ctypes.Structure):
    """Wrapper for the magnetic field data in B_field.h."""

    _fields_ = [
        ("BTC", ctypes.POINTER(cartesian.Struct)),
        ("BGS", ctypes.POINTER(analytical.Struct)),
        ("B2DS", ctypes.POINTER(axisymmetric.Struct)),
        ("B3DS", ctypes.POINTER(perturbed.Struct)),
        ("BSTS", ctypes.POINTER(stellarator.Struct)),
        ("type", ctypes.c_uint32),
    ]


# pylint: disable=too-many-ancestors
class CreateBfieldMixin(
    cartesian.CreateMixin,
    analytical.CreateMixin,
    axisymmetric.CreateMixin,
    perturbed.CreateMixin,
    stellarator.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create magnetic field input.

    This class just combines all the magnetic field mixin classes.
    """

__all__  = [
    "CreateBfieldMixin",
    "BfieldCartesian",
    "BfieldAnalytical",
    "BfieldStellarator",
    "Bfield2D",
    "Bfield3D",
    ]
