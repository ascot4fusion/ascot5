"""Utility tools that are i) not related to physics, ii) not specific to ASCOT5,
and iii) which don't fit anywhere else.
"""

from .imstorage import ImmutableStorage
from .optionaldep import OptionalDependency
from .stringop import decorate, undecorate, format2universaldate
from .varop import (
    validate_variables, check_abscissa, ArrayLike, Numerical, Scalar,
    )
#dt = datetime.strptime(s, "%Y-%m-%d %H:%M:%S.%f")
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
