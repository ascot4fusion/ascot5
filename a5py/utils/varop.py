import warnings
from typing import Union, Sequence, Any, Sized
from contextlib import contextmanager

import unyt
import numpy as np
from numpy.typing import NDArray

Scalar = Union[int, float]
"""Type of a scalar."""

ArrayLike = Union[list[Scalar], np.ndarray, unyt.unyt_array]
"""Type of a numerical array."""

Numerical = Union[
    int, float,
    Sequence[int], Sequence[float],
    NDArray[np.float64],
    unyt.unyt_array]
"""Type of numerical data."""


class VariableValidator:
    def __init__(self):
        self.unitless = []

    def validate(
            self,
            name: str,
            value: int | float | list | np.ndarray | unyt.unyt_array | None,
            expected_shape: tuple[int],
            units: str=None,
            dtype: str="f8",
            default: np.ndarray=None
            ) -> np.ndarray | unyt.unyt_array:
        """Validate a variable.

        This function checks that the variable has or can be cast into the
        expected format (shape, units, dtype) and returns the variable cast.

        Parameters
        ----------
        name : str
            Name of the variable (for generating error messages).
        value : int | float | list | np.ndarray | unyt.unyt_array
            Value of the variable.
        expected_shape : tuple[int], *optional*
            Expected shape of the variable if applicable.

            Special case is ``(1,)`` for when the argument is expected to be 1D
            vector with undefined length.
        units : str, *optional*
            Units of the variable.
        dtype : str, *optional*
            Data type of the variable.
        default : int | float | list | np.ndarray | unyt.unyt_array, *optional*
            If the value is None, this value cast in the expected format is
            returned.

        Returns
        -------
        *cast_variable : np.ndarray | unyt.unyt_array
            The variable cast to the expected format.
        """
        if value is None:
            if default is None:
                raise ValueError(
                    "Default value must be provided for optional arguments: "
                    f"{name}"
                    )
            if units is not None:
                default = unyt.unyt_array(default, units)
            value = default

        if isinstance(value, (int, float, list, tuple)):
            arr = np.array(value, dtype=dtype)
        elif isinstance(value, (np.ndarray, unyt.unyt_array)):
            arr = value.astype(dtype)
        else:
            raise TypeError(
                f"Argument '{name}' must be int, float, list, numpy array, or "
                f"unyt array, not '{type(value)}'"
                )

        if expected_shape == ():
            if arr.size != 1:
                raise ValueError(
                    f"Argument '{name}' has incompatible shape '{arr.shape}', "
                    f"expected scalar."
                    )
        elif expected_shape == (-1,):
            if arr.ndim != 1 and not ( arr.ndim == 2 and 1 in arr.shape ):
                raise ValueError(
                    f"Argument '{name}' has incompatible shape '{arr.shape}', "
                    f"expected 1D vector."
                    )
        elif( arr.size > 1 and
            arr.shape != expected_shape and
            arr.T.shape != expected_shape ):
            raise ValueError(
                f"Argument '{name}' has incompatible shape '{arr.shape}', "
                f"expected: {expected_shape}."
                )
        arr = arr.reshape(expected_shape)

        if units is not None:
            if not hasattr(value, "units"):
                if units != "1":
                    self.unitless.append((name, units))
                arr = unyt.unyt_array(arr, units)
            else:
                try:
                    arr.convert_to_units(units)
                except unyt.exceptions.UnitConversionError:
                    raise ValueError(
                        f"Argument '{name}' has incompatible units "
                        f"'{arr.units}', expected '{units}'."
                        ) from None
        return arr

    def flush(self):
        if self.unitless:
            names = ", ".join([f"'{name}'" for name, _ in self.unitless])
            units = ", ".join([f"'{unit}'" for _, unit in self.unitless])
            warnings.warn(
                f"Argument(s) {names} given without dimensions. Assumed {units}.",
                stacklevel=5,
                )
            self.unitless.clear()


@contextmanager
def validate_variables():
    validator = VariableValidator()
    try:
        yield validator
    finally:
        validator.flush()


def size(x: Any) -> int:
    """Return the size of x.

    Separate function is needed because unyt.array.quantity has size but
    not len.
    """
    try:
        return x.size
    except AttributeError:
        pass
    return len(x)


def to_array(x: Any, n: int) -> unyt.unyt_array | np.ndarray:
    """Convert x to an array of length n."""
    if x is None:
        return None
    if not isinstance(x, Sized):
        return np.full(n, x)

    if size(x) == 1:
        if hasattr(x, "units"):
            return np.full(n, x) * x.units
        return np.full(n, x)
    return x


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