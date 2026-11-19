"""Package for generating inputs from templates and imported data."""

from .transmute import TransmuteAnalyticalBfieldToSplines
from .analyticalinputs import (
    FlatPlasma,
    PremadeMagneticField,
)
from .import_atomic import ImportAtomicData
from .poincare import PoincareMarker, PoincareOptions
from .optionsxml import make_simple_type, make_element_block, doc, make_schema

__all__ = [
    "TransmuteAnalyticalBfieldToSplines",
    "PremadeMagneticField",
    "FlatPlasma",
    "PoincareMarker",
    "PoincareOptions",
    "ImportAtomicData",
]
