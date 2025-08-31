"""Common interface to access data independently from how it is stored."""
import ctypes
from enum import Enum
from typing import Dict, TypeVar
from abc import ABC, abstractmethod

import numpy as np

T = TypeVar('T', bound='DataHolder')


class Format(Enum):
    """Data storage formats."""
    HDF5: int = 1
    """Data is stored in HDF5 file."""
    MEMORY: int = 2
    """Data is stored in memory."""


# pylint: disable=too-few-public-methods
class DataHolder(ABC):
    """An abstract base class to abstract the data storage.

    The data is accessed through this class's properties and the actual data can
    be stored either in memory or in the HDF5 file.

    This is an abstract class. The child of this class should specify the
    `ctypes.Structure` class representing the stored data. The properties of
    the child specify how the data is read from the struct (and the offload
    array) and from the file. The class methods specify how the object is
    initialized from the file or from the given `nd.array`'s.

    Attributes
    ----------
    _staged : bool
        Whether the data is initialized in C.
    _format : Format
        The format in which the data is currently stored.
    _struct_ : ctypes.Structure
        The struct defining the data stored in C.
    """

    def __init__(self, struct : ctypes.Structure, *args, **kwargs) -> None:
        """Initialize an object that does not yet contain any data.

        Parameters
        ----------
        struct : ctypes.Structure
            Python wrapper for the C header corresponding to this data.
        """
        super().__init__(*args, **kwargs)
        self._staged: bool = False
        self._format: Format = Format.MEMORY
        self._struct_: ctypes.Structure = struct

    def _slice_array(self, name, start: int, end: int) -> np.ndarray:
        """Return a copied slice of data array.

        Parameters
        ----------
        name : str
            Name of the array in the struct.
        start : int
            Start index of the slice (inclusive).
        end : int
            End index of the slice (non-inclusive).
        """
        array = getattr(self._struct_, name)
        return np.ctypeslib.as_array(array[start:end]).copy()

    @abstractmethod
    def _export_hdf5(self) -> None:
        """Export the data to the HDF5 file."""

    @abstractmethod
    def export(self) -> Dict[str, np.ndarray]:
        """Return a dictionary with sufficient data to duplicate this instance.

        Returns
        -------
        data : dict[str, np.ndarray or unyt.unyt_array]
            Data that can be passed to the create method to duplicate this
            instance.
        """
        return {}

    @abstractmethod
    def stage(self):
        """Make this data ready for simulation or evaluation.

        The memory consumption may increase significantly, so remember to
        unstage afterwards with :meth:`unstage`.
        """

    @abstractmethod
    def unstage(self):
        """Undo the effect of :meth:`stage` and free consumed memory."""