from functools import wraps
import warnings

import unyt
import numpy as np

import warnings
from typing import Union, Tuple, List, Dict, Any
from functools import wraps

from ..utils import Numerical
from ..exceptions import AscotUnitWarning


def parse_units(strip: bool = False, **expected: str) -> Any:
    """Check the units of the arguments before they are passed to the wrapped
    function.

    Parameters
    ----------
    strip : bool, optional
        Strip units after they have been checked so that the wrapped function
        receives unitless arguments (but implicitly in correct units).
    **expected : dict[str, str]
        Name of the argument and the expected unit, e.g. `r="m"`.
    """
    def actualdecorator(func):
        """Check that the arguments that are expected to have units have the
        correct units.

        If an argument is dimensionless, but it is expected to have units, the
        units are assigned and a warning is displayed.

        Parameters
        ----------
        fun : callable
            The wrapped function.
        """
        argument_names = list(func.__code__.co_varnames)
        @wraps(func)
        def wrapper(*args, **kwargs):
            """Wrapper that parses the units of the arguments before they are
            passed.
            """
            args, kwargs, assumed_units = match_units_arguments(
                args, kwargs, argument_names, expected, strip
                )
            if assumed_units:
                names, units = "", ""
                for name, unit in assumed_units.items():
                    if unit == "1":
                        # No need to warn if the quantity is dimensionless
                        continue
                    names += f"'{name}', "
                    units += f"'{unit}', "
                if names:
                    names, units = names[:-2], units[:-2]
                    warnings.warn(
                        f"Argument(s) {names} given without dimensions "
                        f"(assumed {units})", AscotUnitWarning, stacklevel=2,
                        )
            return func(*args, **kwargs)
        return wrapper
    return actualdecorator


def match_units_arguments(
        args: List[np.ndarray],
        kwargs: Dict[str, np.ndarray],
        argument_names: List[str],
        expected: Dict[str, str],
        strip: bool,
        ) -> Tuple[List[unyt.unyt_array], Dict[str, unyt.unyt_array],
                   Dict[str, str]]:
    """Check that the function arguments have the expected units.

    This function is mostly meant to be used internally. Consider using the
    `parseunits` decorator instead.

    Parameters
    ----------
    args : List[np.ndarray]
        Positional arguments passed to the function.
    kwargs : Dict[str, np.ndarray]
        Keyword arguments passed to the function.
    argument_names : List[str]
        Names of args or all arguments (as they appear in the function
        signature).
    expected : Dict[str, str]
        Name of the argument and the expected unit.
    strip : bool
        After checking and converting the value to expected units, strip units
        and return an unitless nd.array instead.

    Returns
    -------
    parsed_args : Union[List[unyt.unyt_array], List[ndarray]]
        `*args` with the expected units.
    parsed_kwargs : Union[Dict[str, unyt.unyt_array], Dict[str, ndarray]]
        `**kwargs` with the expected units.
    assumedunits : Dict[str, str]
        Names of the arguments, that had no units, and which units they were
        assumed to be in (i.e. the expected units).
    """
    assumedunits = {}
    parsed_args = []
    for arg, name in zip(args, argument_names):
        if name in expected:
            unit = expected[name]
            parsed_args.append(match_units(arg, unit, strip=strip))
            if not hasattr(arg, "units"):
                assumedunits[name] = unit
        else:
            parsed_args.append(arg)
    parsed_kwargs = {}
    for name, arg in kwargs.items():
        if name in expected:
            unit = expected[name]
            parsed_kwargs[name] = match_units(arg, unit, strip=strip)
            if not hasattr(arg, "units"):
                assumedunits[name] = unit
        else:
            parsed_kwargs[name] = arg
    return parsed_args, parsed_kwargs, assumedunits


def match_units(
        value: Numerical,
        unit: Union[str, unyt.Unit],
        strip: bool = False,
        auto_assign: bool = True,
        ) -> np.ndarray:
    """Check that the units of the value are what is expected.

    In case the value is given without units, it is assumed to be in the
    expected units already.

    Parameters
    ----------
    value : ndarray
        Numerical data with or without units.
    unit : str or unyt.Unit
        Expected units.
    strip : bool, optional
        After checking and converting the value to expected units, strip units
        and return an unitless nd.array instead.
    auto_assign : bool, optional
        Automatically assign the expected units to the value if it does not
        have any units.

    Returns
    -------
    matched : ndarray
        The value in expected units.

    Raises
    ------
    ValueError
        If the units of the value does not match the expected.
    """
    if isinstance(unit, str):
        unit = unyt.unyt_quantity.from_string(unit).units
    if not hasattr(value, "units"):
        if auto_assign:
            value *= unit
        else:
            raise ValueError(
                f"Value has incorrect units: expected {unit} but got "
                f"dimensionless"
            )
    if value.units.dimensions != unit.dimensions:
        raise ValueError(
            f"Value has incorrect dimensions: expected {unit} but got "
            f"{value.units.dimensions}"
        )
    # Converting integer quantities might fail silently (bug?), so we use floats
    value = value.astype("f8", copy=False)
    if strip:
        return value.to(unit).v
    return value.to(unit)
