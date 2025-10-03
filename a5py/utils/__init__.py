"""Utility tools.

This module contains any tools that are

- not directly related to physics
- not specific to ASCOT5
- and which don't fit anywhere else.

.. rubric:: Array manipulation

.. autofunction:: validate_variables

.. autofunction:: validate_abscissa

.. autofunction:: size

.. autofunction:: scalar2array


.. rubric:: String manipulation

.. autofunction:: decorate

.. autofunction:: undecorate

.. autofunction:: format2universaldate


.. rubric:: Misc

.. autodata:: Scalar

.. autodata:: ArrayLike

.. autoclass:: OptionalDependency
    :members:

.. autoclass:: ImmutableStorage
    :members:

"""
from .imstorage import ImmutableStorage
from .optionaldep import OptionalDependency
from .stringop import decorate, undecorate, format2universaldate
from .arrayop import (
    ArrayLike, Scalar, validate_variables, validate_abscissa, size, scalar2array,
    )

__all__ = [
    "Scalar",
    "ArrayLike",
    "ImmutableStorage",
    "OptionalDependency",
    "size",
    "decorate",
    "undecorate",
    "scalar2array",
    "validate_abscissa",
    "validate_variables",
    "format2universaldate",
    ]
