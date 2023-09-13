from functools import wraps
import warnings

import unyt
import numpy as np

def parseunits(strip=False, **units):
    """Prepare arguments that are expected to have physical units.

    This decorator:

    - Makes sure every argument has expected dimensions.
    - Assigns units if units were not provided but they were expected.
    - Strips units if asked (after checking/assignment).

    Examples:

    :: code-block python

       @parseunits(x="m", strip=True)
       def fun(x):
           pass

    Parameters
    ----------
    strip : bool, optional
        Strip units so that dimensionless quantities (but in expected units)
        are passed to the function.
    **units
        Argument in wrapped function and the expected unit as a string.
    """
    for k in units.keys():
        # Convert unit strings to unyts
        units[k] = unyt.unyt_quantity.from_string(units[k]).units

    def actualdecorator(fun):
        """Parse fun arguments that are expected to have units.
        """
        argnames = list(fun.__code__.co_varnames)
        if "self" in argnames:
            argnames.remove("self")

        def checkandstrip(val, unit, name, assignedunits):
            """Check units of a given argument and strip if needed.
            """
            dim = unit.dimensions
            try:
                # Try to get argument units
                valdim = val.units.dimensions
            except AttributeError:
                # Argument doesn't have units, assign and add warning
                val = val*unit
                valdim = dim
                assignedunits[name] = unit

            if valdim != dim:
                raise ValueError(
                    "\"%s\" has incorrect dimensions: expected %s but got %s" %
                    (name, dim, valdim))
            if strip:
                return val.to(unit).v
            return val

        @wraps(fun)
        def wrapper(*args, **kwargs):
            """Replace args and kwargs with parsed arguments when necessary.
            """
            parsedargs = [None]*len(args)
            assignedunits = {} # Dimensionless units with assigned units
            for i, name in enumerate(argnames):
                if i == len(args): break
                if name in units:
                    parsedargs[i] = checkandstrip(args[i], units[name],
                                                  name, assignedunits)
                else:
                    parsedargs[i] = args[i]
            for name, val in kwargs.items():
                if name in units:
                    kwargs[name] = checkandstrip(kwargs[name], units[name],
                                                 name, assignedunits)
            if len(assignedunits) > 0:
                msg1 = "Argument(s) "
                msg2 = " given without dimensions (assumed "
                for name, unit in assignedunits.items():
                    msg1 += name + ", "
                    msg2 += str(unit) + ", "
                msg = msg1[:-2] + msg2[:-2] + ")"
                warnings.warn(msg, UserWarning)

            return fun(*parsedargs, **kwargs)

        return wrapper


    return actualdecorator
