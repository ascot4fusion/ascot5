"""Common interface to access data independently from how it is stored."""
import ctypes
from enum import IntFlag, auto, Enum
from typing import Dict, TypeVar
from abc import ABC, abstractmethod

import numpy as np
import unyt

T = TypeVar('T', bound='DataHolder')

class Status(IntFlag):
    SAVED = auto()
    STAGED = auto()

    def __contains__(self, item):
        """Allow `item in Status.BOTH` checks."""
        # If this status is a group (tuple of values)
        if isinstance(self.value, tuple):
            return item.value in self.value
        return item == self


class Format(Enum):
    """Data storage formats."""
    HDF5: int = 1
    """Data is stored in HDF5 file."""
    MEMORY: int = 2
    """Data is stored in memory."""

class DataStruct(ctypes.Structure):

    def readonly_carray(self, name, shape, units=None):
        arr = np.ctypeslib.as_array(getattr(self, name), shape=shape)
        arr.flags.writeable = False
        if not units is None:
            return unyt.unyt_array(arr, units)
        return arr


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
