"""This module contains all plasma input variants.

.. autosummary::
    :nosignatures:

    ~PlasmaLinear1D
    ~CreatePlasmaMixin.create_plasmalinear1d

.. rubric:: Classes

.. autoclass:: PlasmaLinear1D
    :members:

.. autoclass:: CreatePlasmaMixin
    :members:
    :inherited-members:
"""
import ctypes

from a5py.libascot import input_category

from . import linear1d
from . import dynamic1d
from .linear1d import PlasmaLinear1D
from .dynamic1d import PlasmaDynamic1D


@input_category
class Plasma(ctypes.Structure):
    """Wrapper for the plasma data in libascot.so."""

    _fields_ = [
        ("linear1d", ctypes.POINTER(linear1d.Struct)),
        ("dynamic1d", ctypes.POINTER(dynamic1d.Struct)),
        ("type", ctypes.c_uint32),
        ]


# pylint: disable=too-many-ancestors
class CreatePlasmaMixin(
    linear1d.CreateMixin,
    #radialdynamic.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create plasma input.

    This class just combines all the plasma mixin classes.
    """

__all__  = [
    "CreatePlasmaMixin",
    "PlasmaLinear1D",
    "PlasmaDynamic1D",
    ]
