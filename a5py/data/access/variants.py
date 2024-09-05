"""Base classes for input and run variants."""
import unyt
import numpy as np

from . import metadata
from .treeparts import Leaf, MetaData
from .tree import RunVariant
from .dataholder import DataHolder
from ... import utils
from ... import physlib


class InputVariant(Leaf, DataHolder):
    """Base class for input variants."""

    def _read_hdf5(self, name):
        """Read dataset from HDF5 file.

        This is just a helper method which can be called when implementing the
        properties in subclasses.

        Parameters
        ----------
        name : str
            Name of the dataset to read.

        Returns
        -------
        data : np.ndarray or unyt.unyt_array
            The data.
        """
        hdf5manager = self._treemanager.hdf5manager
        return hdf5manager.read_datasets(self.qid, self.variant, name)

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


def validate_parameters(
        parameters, names, units, shape, dtype, required, default,
        ):
    """Validate parameters.

    Parameters
    ----------
    parameters : dict
        The parameters to validate.
    names : list[str]
        Names of the parameters.

        The names should match the keys in the `parameters` dictionary. The
        order of `units`, `required`, `shape`, `dtype`, and `default` should
        match the order of `names`.
    units : list[strings]
        Expected units.
    shape : list[int]
        Expected shapes.
    dtype : list[str]
        Expected dtype.
    required : bool
        Flag indicating whether the parameters are required.

        For required parameters dummy values are used only if every parameter is
        None. If only some are None, exception is raised instead.
    default : list[unyt.unyt_quantity]
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
