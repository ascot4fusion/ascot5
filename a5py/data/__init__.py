"""Contains all data classes for simulation inputs and outputs, and the main
data structure class.
"""
from .access import InputVariant

from .base import AscotData
from .bfield import Bfield
from .efield import Efield
from .plasma import Plasma
from .neutral import Neutral
from .wall import Wall
from .mhd import Mhd
from .boozer import BoozerMap
from .atomic import Atomic
from .nbi import NbiStruct

__all__ = [
    "InputVariant",
    "AscotData",
    "Bfield",
    "Efield",
    "Plasma",
    "Neutral",
    "Wall",
    "Mhd",
    "BoozerMap",
    "Atomic",
    "NbiStruct",
]
