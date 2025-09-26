"""Utility tools that are i) not related to physics, ii) not specific to ASCOT5,
and iii) which don't fit anywhere else.
"""

from .imstorage import ImmutableStorage
from .optionaldep import OptionalDependency
from .stringop import decorate, undecorate, format2universaldate
from .varop import (
    ArrayLike, Numerical, Scalar, validate_variables, check_abscissa, size,
    to_array,
    )

__all__ = [
    "Scalar",
    "ArrayLike",
    "Numerical",
    "ImmutableStorage",
    "OptionalDependency",
    "decorate",
    "undecorate",
    "check_abscissa",
    "validate_variables",
    "format2universaldate",
    ]
