"""Utility tools that are i) not related to physics, ii) not specific to ASCOT5,
and iii) which don't fit anywhere else.
"""
import re
import importlib
from datetime import datetime
from typing import Tuple, List, Union

import numpy as np
import unyt


Scalar = Union[int, float]
"""Type of a scalar."""

ArrayLike = Union[List[Scalar], np.ndarray, unyt.unyt_array]
"""Type of a numerical array."""

Numerical = Union[Scalar, ArrayLike]
"""Type of numerical data."""


def format2universaldate(date : datetime) -> str:
    """Convert a datetime object to a string in the format ASCOT5 uses.

    Parameters
    ----------
    date : datetime.datetime
        The datetime object to be converted.

    Returns
    -------
    formatted_date : str
        The datetime object converted to a string in the format
        "YYYY-MM-DD HH:MM:SS".
    """
    return date.strftime("%Y-%m-%d %H:%M:%S")


def validate_variables(
        variables : List[Numerical] | Numerical,
        names : List[str] | str,
        shape : List[Tuple[int]] | Tuple[int],
        dtype : List[str] | str,
        ) -> List[Numerical] | Numerical:
    """Validate the shape and dtype of multiple variables.

    This function checks that a variable can be inferred to have the expected
    format (e.g. it can be cast into the expected shape) and returns the
    variables cast exactly in the expected format.

    Parameters
    ----------
    variables : array_like or list[array_like]
        The variable(s) to validate.

        If dictionary is provided, the variables are cast in place.
    names : string or list of strings
        Names of the variables (required to generate error message).
    shape : tuple of ints or a list of tuples of ints
        The expected shape for all variables or for each variable separately.
    dtype : type
        The expected type for all variables or for each variable separately.

    Returns
    -------
    *cast_variables : array_like
        The variables cast to the expected format.

    Raises
    ------
    ValueError
        If not able to infer a variable to have the expected format.
    """
    if not isinstance(variables, list):
        variables = [variables]
        if not isinstance(names, list):
            names = [names]
    if not isinstance(shape, list):
        shape = [shape] * len(variables)
    if not isinstance(dtype, list):
        dtype = [dtype] * len(variables)

    cast_variables = []
    for var, name, expshape, exptype in zip(variables, names, shape, dtype):
        varshape = var.shape
        try:
            var = np.reshape(var, expshape)
        except ValueError:
            raise ValueError(
                f"Argument {name} has incompatible shape {varshape}, "
                f"expected: {expshape}."
                ) from None
        vartype = var.dtype
        var_original, var = var, var.astype(exptype)
        equivalent = ( var_original == var if varshape == ()
                      else (var_original == var).all() )
        if not equivalent:
            raise ValueError(
                f"Argument {name} has incompatible type {vartype}, "
                f"expected: {exptype}."
                )
        cast_variables.append(var)
    if len(cast_variables) == 1:
        return [cast_variables[0]]
    return tuple(cast_variables)


def decorate(
        string: str,
        color: str | None = None,
        bold: bool = False,
        underline: bool = False,
        ) -> str:
    """Make text underlined, bold, and/or colourful.

    Parameters
    ----------
    string : str
        String to be decorated.
    color : {"green", "purple"}, optional
        Color the text.
    bold : bool, optional
        Bold the text.
    underline : bool, optional
        Underline the text.

    Returns
    -------
    decorated : str
        The decorated string.
    """
    colors = {
        "green":"\033[32m",
        "purple":"\033[35m",
    }
    effects = {
        "reset":"\033[0m",
        "bold":"\033[01m",
        "underline":"\033[04m",
    }

    if color:
        if color not in colors:
            raise ValueError(
                f"Unknown color: {color}. "
                f"Available colors are: {colors.keys()}."
                )
        color = colors[color]
    else:
        color = ""

    bolded = effects["bold"] if bold else ""
    underlined = effects["underline"] if underline else ""
    return f"{bolded}{underlined}{color}{string}{effects['reset']}"


def undecorate(string: str) -> str:
    """Remove decorations (ANSI escape sequences) from a string.

    Parameters
    ----------
    string : str
        A string to undecorate.

    Returns
    -------
    undecorated : str
        The undecorated string.
    """
    return re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])").sub("", string)


def check_abscissa(abscissa, name, uniform=True, periodic=False) -> None:
    """Check if the abscissa is increasing and otherwise valid.

    Parameters
    ----------
    abscissa : array_like
        The abscissa to check.
    name : str
        Name of the abscissa (for the possible error message).
    uniform : bool, optional
        If True, check if the abscissa is also uniform.
    periodic : bool, optional
        If True, check if the abscissa is also periodic.

    Raises
    ------
    ValueError
        If the abscissa is invalid.
    """
    nx, dx = abscissa.size, np.diff(abscissa)
    if nx < 2:
        raise ValueError(f"{name} must have at least two points.")
    if np.any(dx < 0):
        raise ValueError(f"{name} must be increasing.")
    dx_expected = np.diff( abscissa[[0,-1]] ) / ( nx - 1 )
    if uniform and not np.allclose(dx.v, dx_expected.v, atol=0):
        raise ValueError(f"{name} must be uniform.")
    xmax, xmin = abscissa[-1] + dx[0], abscissa[0]
    nperiod = ( 360 * unyt.deg / (xmax - xmin) ).v
    if periodic and not np.isclose(nperiod, np.round(nperiod)):
        raise ValueError(
            f"{name} must be periodic (the current interval was [{xmin:.2f}, "
            f"`{xmax:.2f}`] i.e. `nperiod = {nperiod:.2f}`)."
            )


class OptionalDependency:
    """A lazy, optional import wrapper for third-party packages.

    This class allows you to declare an optional dependency without causing
    an ImportError at module import time. Instead, the dependency is loaded
    only when first accessed. If the dependency is not available in the
    environment, an informative ImportError is raised *at the point of use*.

    Example
    -------
    Instead of::

        import optionalpackage as optional

    use::

        optional = OptionalDependency("optionalpackage")

    Then, when calling::

        optional.do_stuff()

    an ImportError will be raised if ``optionalpackage`` is missing.
    If the package is present, it will be imported on demand and cached for
    subsequent use.

    Benefits
    --------
    - **Graceful degradation**: code can be imported and executed even if
      optional dependencies are not installed, as long as they are not used.
    - **Lazy loading**: dependencies are imported only when needed, which
      can improve startup time.
    - **Caching**: once imported, the module object is cached, so repeated
      accesses incur no extra overhead.
    """

    def __init__(self, name: str):
        self._name = name
        self._module = None
        self._loaded = False

    def _load(self):
        if not self._loaded:
            try:
                self._module = importlib.import_module(self._name)
            except ImportError:
                raise ImportError(
                    f"Optional dependency '{self._name}' is not installed. "
                    f"Please install it to use this feature."
                )
            self._loaded = True
        return self._module

    def __getattr__(self, item):
        return getattr(self._load(), item)
