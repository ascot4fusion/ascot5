from functools import wraps
from itertools import chain
import wrapt
import warnings

import unyt
import numpy as np
import inspect


def prepareunitarrays(ravel=False, **arg_units):
    """
    Prepare argument arrays. Use as a decorator.

    This function ravels any numerical input into a 1D numpy array. If units
    are assigned to this function, those units are applied automatically
    (and warning is issued) if the inputs arrays don't have units already.

    Args:
        ravel : bool, optional <br>
            Flag indicating whether the argument arrays/scalars are converted to
            1D arrays.
        arg_units : unyt.unit <br>
            Name of each argument (identical to what are used in the function
            this decorator is applied to) to which this decorator applies. As a
            value, use the default unit (unyt.unit) for that argument. Use None
            if default units are not applied for that argument.

    Examples:

    Set default unit "meters" for R, ravel R and ids, and leave endcond
    unmodified.
    >>> @checkunits(R=unyt.m, ids=None)
    ... fun(R, ids, endcond)
    """
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):

        names_of_args = list(wrapped.__code__.co_varnames)
        if "self" in names_of_args:
            names_of_args.remove("self")

        def _execute(*_args, **_kwargs):
            preparedargs = [None]*len(_args)

            for i in range(len(_args)):
                preparedargs[i] = _args[i]

                if names_of_args[i] in arg_units.keys():
                    if ravel:
                        preparedargs[i] = \
                        np.asarray(preparedargs[i]).ravel().astype(dtype="f8")

                    if arg_units[names_of_args[i]] is not None \
                       and type(preparedargs[i]) != unyt.array.unyt_array:
                        preparedargs[i] *= arg_units[names_of_args[i]]
                        warnings.warn(
                            "Units not given for " + names_of_args[i] +
                            ".\nUsing default unit " +
                            str(arg_units[names_of_args[i]]) + " instead.")

            return wrapped(*preparedargs, **_kwargs)

        return _execute(*args, **kwargs)

    return wrapper


def stripunits(base, *arr):
    """
    Convert to desired unit system and strip units.

    This decorator converts all arguments that have units to a given unit system
    (e.g. SI/mks) and strips units so that the arguments can be treated as an
    ordinary numpy arrays.

    Args:
        base : str <br>
            Unit base to which arguments are converted before stripping.
        arr : str, optional <br>
            Explicitly name the arguments to be stripped. Others are left as is.
    """
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        names_of_args = list(wrapped.__code__.co_varnames)

        def _execute(*_args, **_kwargs):
            strippedargs = [None]*len(_args)
            for i in range(len(_args)):
                strippedargs[i] = _args[i]
                if not (len(arr) == 0 or names_of_args[i] in arr):
                    continue
                if type(_args[i]) == unyt.array.unyt_array \
                   or type(_args[i]) == unyt.unyt_quantity:
                    strippedargs[i].convert_to_base(base)
                    strippedargs[i] = np.array(strippedargs[i])

            return wrapped(*strippedargs, **_kwargs)

        return _execute(*args, **kwargs)

    return wrapper



def accepts(**arg_units):
    """
    Raise error if incorrect units are used. Wrapt package compatible unyt.accepts().

    This function is identical to unyt.accepts() except this one is wrapped
    using the wrapt package making it compatible with other decorators that
    also wrap functions using the wrapt package.
    """

    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        names_of_args = list(wrapped.__code__.co_varnames)
        if "self" in names_of_args:
            names_of_args.remove("self")

        def _execute(*_args, **_kwargs):
            for arg_name, arg_value in chain(zip(names_of_args, _args), _kwargs.items()):
                if arg_name in arg_units:  # function argument needs to be checked
                    dimension = arg_units[arg_name]
                    if not _has_dimensions(arg_value, dimension):
                        raise TypeError(
                            "arg '%s=%s' does not match %s"
                            % (arg_name, arg_value, dimension)
                            )
            return wrapped(*_args, **_kwargs)

        return _execute(*args, **kwargs)

    return wrapper


def prepareandcheckunits(ravel=False, base=None, **arg_units):
    """
    Combines accepts(), preparedunitarrays(), and stripunits().

    ravel : bool, optional <br>
            Flag indicating whether the argument arrays/scalars are converted to
            1D arrays.
    base : str <br>
            Unit base to which arguments are converted before stripping. If
            None, the stripping is not performed.
    arg_units : unyt.unit <br>
            Name of each argument (identical to what are used in the function
            this decorator is applied to) to which this decorator applies. As a
            value, use the default unit (unyt.unit) for that argument. Use None
            if default units are not applied for that argument.
    """
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):

        names_of_args = list(wrapped.__code__.co_varnames)
        if "self" in names_of_args:
            names_of_args.remove("self")

        def _execute(*_args, **_kwargs):
            preparedargs = [None]*len(_args)

            # Go through _args
            for i in range(len(_args)):
                preparedargs[i] = _args[i]

                if names_of_args[i] in arg_units.keys():
                    if ravel:
                        unit = 1
                        try:
                            unit = preparedargs[i].units
                        except:
                            pass
                        preparedargs[i] = \
                        np.asarray(preparedargs[i]).ravel().astype(dtype="f8")*unit

                    if arg_units[names_of_args[i]] is not None \
                       and (type(preparedargs[i]) != unyt.array.unyt_array
                            and type(preparedargs[i]) != unyt.unyt_quantity):
                        preparedargs[i] *= arg_units[names_of_args[i]]
                        warnings.warn(
                            "Units not given for " + names_of_args[i] +
                            ".\nUsing default unit " +
                            str(arg_units[names_of_args[i]]) + " instead.")

            # Same process for _kwargs, there is probably a cleaner way to do this
            for k in _kwargs.keys():

                if k in arg_units.keys():
                    if ravel:
                        _kwargs[k] = \
                        np.asarray(_kwargs[k]).ravel().astype(dtype="f8")

                    if arg_units[k] is not None and _kwargs[k] is not None \
                       and (type(_kwargs[k]) != unyt.array.unyt_array
                            and type(_kwargs[k]) != unyt.unyt_quantity):
                        _kwargs[k] *= arg_units[k]
                        warnings.warn(
                            "Units not given for " + k +
                            ".\nUsing default unit " +
                            str(arg_units[k]) + " instead.")

            for arg_name, arg_value in chain(zip(names_of_args, preparedargs),
                                             _kwargs.items()):
                if arg_name in arg_units and arg_value is not None:
                    dimension = arg_units[arg_name].dimensions
                    if not _has_dimensions(arg_value, dimension):
                        raise TypeError(
                            "arg '%s=%s' does not match %s"
                            % (arg_name, arg_value, dimension)
                            )

            if base is not None:
                for i in range(len(_args)):
                    if type(preparedargs[i]) == unyt.array.unyt_array \
                       or type(preparedargs[i]) == unyt.unyt_quantity:
                        preparedargs[i].convert_to_base(base)
                        preparedargs[i] = np.array(preparedargs[i])

                for k in _kwargs.keys():
                    if type(_kwargs[k]) == unyt.array.unyt_array \
                       or type(_kwargs[k]) == unyt.unyt_quantity:
                        _kwargs[k].convert_to_base(base)
                        _kwargs[k] = np.array(_kwargs[k])

            return wrapped(*preparedargs, **_kwargs)

        return _execute(*args, **kwargs)

    return wrapper


def _has_dimensions(quant, dim):
    """
    Helper function for accepts().
    """
    try:
        arg_dim = quant.units.dimensions
    except AttributeError:
        arg_dim = None
    return arg_dim == dim
