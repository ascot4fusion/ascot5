"""This module contains all plasma input variants.

.. autosummary::
    :nosignatures:

    ~Plasma1D
    ~CreatePlasmaMixin.create_plasma1d

.. rubric:: Classes

.. autoclass:: Plasma1D
    :members:

.. autoclass:: CreatePlasmaMixin
    :members:
    :inherited-members:
"""
import ctypes

from a5py.libascot import input_category

from . import radial
from . import radialdynamic
from .radial import Plasma1D
from .radialdynamic import Plasma1DDynamic


@input_category
class Plasma(ctypes.Structure):
    """Wrapper for the plasma data in libascot.so."""

    _fields_ = [
        ("plasma_1D", ctypes.POINTER(radial.Struct)),
        ("plasma_1Dt", ctypes.POINTER(radialdynamic.Struct)),
        ("plasma_1DS", ctypes.POINTER(None)),
        ("type", ctypes.c_uint32),
        ]


# pylint: disable=too-many-ancestors
class CreatePlasmaMixin(
    radial.CreateMixin,
    #radialdynamic.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create plasma input.

    This class just combines all the plasma mixin classes.
    """

__all__  = [
    "CreatePlasmaMixin",
    "Plasma1D",
    "Plasma1DDynamic",
    ]
