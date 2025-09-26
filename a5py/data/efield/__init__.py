"""This module contains all electric field input variants.

.. autosummary::
    :nosignatures:

    ~EfieldCartesian
    ~CreateEfieldMixin.create_efieldcartesian

.. rubric:: Classes

.. autoclass:: EfieldCartesian
    :members:

.. autoclass:: CreateEfieldMixin
    :members:
    :inherited-members:
"""
import ctypes

from a5py.libascot import input_category

from . import cartesian
from . import potential1d
from .cartesian import EfieldCartesian
from .potential1d import EfieldPotential1D


@input_category
class Efield(ctypes.Structure):
    """Wrapper for the electric field data in libascot.so."""

    _fields_ = [
        ("ETC", ctypes.POINTER(cartesian.Struct)),
        ("E1DS", ctypes.POINTER(potential1d.Struct)),
        ("type", ctypes.c_uint32),
    ]


# pylint: disable=too-many-ancestors
class CreateEfieldMixin(
    cartesian.CreateMixin,
    #potential1d.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create electric field input.

    This class just combines all the electric field mixin classes.
    """

__all__  = [
    "CreateEfieldMixin",
    "EfieldCartesian",
    "EfieldRadialPotential",
    ]
