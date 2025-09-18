from typing import Union

import numpy as np
import unyt

Scalar = Union[int, float]
"""Type of a scalar."""

ArrayLike = Union[list[Scalar], np.ndarray, unyt.unyt_array]
"""Type of a numerical array."""

Numerical = Union[Scalar, ArrayLike]
"""Type of numerical data."""

def validate_variables(
        variables : list[Numerical] | Numerical,
        names : list[str] | str,
        shape : list[tuple[int]] | tuple[int],
        dtype : list[str] | str,
        ) -> list[Numerical] | Numerical:
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