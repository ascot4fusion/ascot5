"""Base classes for input and run variants."""
import inspect
from typing import Tuple, Any

import unyt
import numpy as np

from ... import utils
from ... import physlib

from . import metadata
from .treeparts import Leaf, MetaData
from .tree import RunVariant
from .dataholder import DataHolder


class InputVariant(Leaf, DataHolder):
    """Base class for input variants."""

    def _read_hdf5(self, *names):
        """Read dataset from HDF5 file.

        This is just a helper method which can be called when implementing the
        properties in subclasses.

        Parameters
        ----------
        *names : str
            Names of the datasets to read.

        Returns
        -------
        *data : np.ndarray or unyt.unyt_array
            The data.
        """
        hdf5manager = self._treemanager.hdf5manager
        if len(names) > 1:
            data = []
            for name in names:
                data.append(
                    hdf5manager.read_datasets(self.qid, self.variant, name)
                    )
            return tuple(data)
        return hdf5manager.read_datasets(self.qid, self.variant, names[0])

    def _from_struct_(self, name, shape=None, units=None):
        """Return a copy of an array in struct.

        This is just a helper method which can be called when implementing the
        properties in subclasses.

        Parameters
        ----------
        name : str
            Name of the array in struct.
        shape : tuple, optional
            Shape of the data.
        units : str, optional
            Units of the data.

        Returns
        -------
        data : np.ndarray or unyt.unyt_array
            A copy of the data.
        """
        unit = unyt.unyt_quantity.from_string(units) if units else 1
        arr = np.ctypeslib.as_array(getattr(self._struct_, name)).copy()
        if shape is not None:
            arr = arr.reshape(shape)
        return arr * unit


def new_metadata(variant, note):
    qid, date, default_note = metadata.generate_metadata()
    if note is None:
        note = default_note
    return MetaData(qid, date, note, variant)


def parse_parameters(*args):
    """Convert parameters to a dictionary with properly cast values.

    Helper function to parse input parameters from all parameters that are
    given to `create_*` functions. The name of the parameters are deduced
    by inspecting the caller function. The values are cast to `np.ndarray` or
    `unyt.unyt_array` if they have units.

    Parameters
    ----------
    *args
        Parameters in the same order as they are in the `create_*` function.

    Returns
    -------
    parameters : dict[str, np.ndarray or unyt.unyt_array]
        The parsed input parameters.
    """
    names = inspect.currentframe().f_back.f_code.co_varnames
    notparameters = ["self", "note", "activate", "dryrun", "store_hdf5"]
    names = [name for name in names if name not in notparameters]
    def cast(x):
        """Ensure parameter is a numpy or unyt array."""
        return x if isinstance(x, unyt.unyt_array) else np.asarray(x)
    return {
        arg:(val if val is None else cast(val)) for arg, val in zip(names, args)
        }


def validate_parameters(
        parameters, names, units, shape, dtype, default, required,
        ):
    """Validate parameters.

    Parameters
    ----------
    *args :
        The parameters to validate.

        The order of `units`, `shape`, `dtype`, and `default` should match the
        order the parameters are given here.
    units : list[strings]
        Expected units.
    shape : list[int]
        Expected shapes.
    dtype : list[str]
        Expected dtype.
    default : list[unyt.unyt_quantity]
        The default values to be used if the parameter is None.
    required : bool
        Flag indicating whether the parameters are required.

        For required parameters default values are used only if every parameter
        is None. If only some are None, exception is raised instead.

    Returns
    -------
    parameters : dict[str, unyt.unyt_quantity]
        The validated parameters.
    """
    dummyvals = {}
    for name, dummy, unit in zip(names, default, units):
        if parameters[name] is None:
            dummyvals[name] = dummy * unyt.unyt_array.from_string(unit)
    if required and len(dummyvals) > 0 and len(dummyvals) < len(names):
        missing = ", ".join([f"{name}" for name in dummyvals])
        raise ValueError(f"Missing required parameter(s) {missing}")
    parameters.update(dummyvals)
    variables = []
    for name, unit in zip(names, units):
        variables.append(
            physlib.match_units(parameters[name], unit, strip=True)
        )
    utils.validate_variables(variables, names, shape, dtype)
    for name, val in zip(names, variables):
        parameters[name] = val
    return parameters


validate_required_parameters = (
    lambda *args, **kwargs: validate_parameters(*args, **kwargs, required=True)
)
"""Call `validate_parameters` assuming that the parameters are required."""

validate_optional_parameters = (
    lambda *args, **kwargs: validate_parameters(*args, **kwargs, required=False)
)
"""Call `validate_parameters` assuming that the parameters are optional."""
