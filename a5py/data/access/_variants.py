"""Base classes for input and run variants."""
import ctypes
import inspect
from typing import Tuple, Any

import unyt
import numpy as np

from .. import cstructs
from ... import utils
from ... import physlib

from . import leaves
from .leaves import Leaf, OutputVariant
from .dataholder import DataHolder

NOTPARAMETERS = ["self", "note", "activate", "dryrun", "store_hdf5"]
"""List of arguments that are passed to the `create_*` functions but which
do not specify the actual data."""


class RunVariant(OutputVariant):
    """Base class for run variants."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class Diagnostic():
    """Base class for diagnostics."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._staged: bool = True

    def read_hdf5(self, hdf5manager):
        """"""
        self._staged = False

    def write_hdf5(self, hdf5manager):
        """"""
        self._staged = False


class InputVariant(Leaf, DataHolder):
    """Base class for input variants."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._staged = False

    def _read_hdf5(self, *names):
        """Read dataset from HDF5 file.

        This is a helper method which can be called when implementing the
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

    def _from_struct_(self, name, shape=None, idx=None, units=None):
        """Return a copy of an array in struct.

        This is just a helper method which can be called when implementing the
        properties in subclasses.

        Parameters
        ----------
        name : str
            Name of the array in struct.
        shape : tuple, optional
            Shape of the data.
        idx : int, optional
            In case the attribute is an array of structs, this is the index on
            the array which is read.
        units : str, optional
            Units of the data.

        Returns
        -------
        data : np.ndarray or unyt.unyt_array
            A copy of the data.
        """
        unit = unyt.unyt_quantity.from_string(units) if units else 1
        attribute = getattr(self._struct_, name)
        if idx is not None:
            attribute = attribute[idx]
        if isinstance(attribute, cstructs.linint1D_data):
            s1 = attribute.n_x
            arr = np.ctypeslib.as_array(attribute.c[:s1]).copy()
            arr = arr.reshape((s1,))
        elif isinstance(attribute, cstructs.linint3D_data):
            s1, s2, s3 = attribute.n_x, attribute.n_y, attribute.n_z
            arr = np.ctypeslib.as_array(attribute.c[:s1*s2*s3]).copy()
            arr = arr.reshape((s1,s2,s3))
        elif isinstance(attribute, cstructs.interp1D_data):
            s1 = attribute.n_x
            arr = np.ctypeslib.as_array(attribute.c[:s1*2:2]).copy()
            arr = arr.reshape((s1,))
        elif isinstance(attribute, cstructs.interp2D_data):
            s1, s2 = attribute.n_x, attribute.n_y
            arr = np.ctypeslib.as_array(attribute.c[:s1*s2*4:4]).copy()
            arr = arr.reshape((s1,s2))
        elif isinstance(attribute, cstructs.interp3D_data):
            s1, s2, s3 = attribute.n_x, attribute.n_y, attribute.n_z
            arr = np.ctypeslib.as_array(attribute.c[:s1*s2*s3*8:8]).copy()
            arr = arr.reshape((s1,s2,s3))
        else:
            arr = np.ctypeslib.as_array(attribute, shape).copy()
        if shape is not None:
            arr = arr.reshape(shape)
        return arr * unit

    def destroy(self):
        if self._staged:
            self.unstage()
        self.destroy()


def new_metadata(variant, note):
    qid, date, default_note = leaves.generate_metadata()
    if note is None:
        note = default_note
    return MetaData(qid, date, note, variant)


def get_parameternames():
    """Get parameter names from the caller function.

    Helper function to separate the parameters that specify input data from
    the common parameters that are passed to the `create_*` functions.

    Parameters
    ----------
    allparameters : list[str]
        List of all parameters that are passed to the `create_*` function.

    Returns
    -------
    names : list[str]
        The names of the parameters.
    """
    names = inspect.currentframe().f_back.f_code.co_varnames
    return [name for name in names if name not in NOTPARAMETERS]


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
    names = [name for name in names if name not in NOTPARAMETERS]
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

        Use empty strings for parameters where units don't apply (e.g. strings).
        The unitless parameters must be listed last.
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
    flag_unitless_encountered = False
    for name, dummy, unit in zip(names, default, units):
        if parameters[name] is None:
            if unit == "":
                dummyvals[name] = dummy
                flag_unitless_encountered = True
            else:
                dummyvals[name] = dummy * unyt.unyt_array.from_string(unit)
                if flag_unitless_encountered:
                    raise ValueError(
                        "Unitless parameters must be listed last."
                    )
    if required and len(dummyvals) > 0 and len(dummyvals) < len(names):
        missing = ", ".join([f"{name}" for name in dummyvals])
        raise ValueError(f"Missing required parameter(s) {missing}")
    parameters.update(dummyvals)
    variables = []
    for name, unit in zip(names, units):
        if unit != "":
            variables.append(
                physlib.match_units(parameters[name], unit, strip=False)
                )
    variables = utils.validate_variables(variables, names, shape, dtype)
    for name, val in zip(names, variables):
        parameters[name] = val
    return parameters


validate_required_parameters = (
    lambda parameters, names, units, shape, dtype, default:
    validate_parameters(
        parameters, names, units, shape, dtype, default, required=True
    )
)
"""Call `validate_parameters` assuming that the parameters are required."""

validate_optional_parameters = (
    lambda parameters, names, units, shape, dtype, default:
    validate_parameters(
        parameters, names, units, shape, dtype, default, required=False
    )
)
"""Call `validate_parameters` assuming that the parameters are optional."""
