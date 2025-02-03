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
    CSTRUCT: int = 2
    """Data is stored in memory using C struct and offload arrays."""


# pylint: disable=too-few-public-methods
class DataHolder(ABC):
    """An abstract base class to abstract the data storage.

    The data is accessed through this class's properties and the actual data can
    be stored either in memory (as C struct and offload array) or in the HDF5
    file.

    This is an abstract class. The child of this class should specify the
    `ctypes.Structure` class representing the stored data. The properties of
    the child specify how the data is read from the struct (and the offload
    array) and from the file. The class methods specify how the object is
    initialized from the file or from the given `nd.array`'s.

    Attributes
    ----------
    _format : Format
        The format in which the data is currently stored.
    _struct_ : ctypes.Structure
        The struct defining the data stored in C.
    _offload_array_ : ctypes.
        Allocatable array of double type data.
    _int_offload_array_ : ctypes.
        Allocatable array of int type data.
    _offload_array_count : int
        The current index on the _offload_array_ where the next array of data
        is amended.
    _in_offload_array_count : int
        The current index on the _int_offload_array_ where the next array of
        data is amended.
    """

    def __init__(self, struct : ctypes.Structure, *args, **kwargs) -> None:
        """Initialize an object that does not yet contain any data.

        Parameters
        ----------
        struct : ctypes.Structure
            Python wrapper for the C header corresponding to this data.
        """
        super().__init__(*args, **kwargs)
        self._format: Format = Format.CSTRUCT
        self._struct_: ctypes.Structure = struct
        self._offload_array_: ctypes.Array[ctypes.c_double]
        self._int_offload_array_: ctypes.Array[ctypes.c_int64]
        self._offload_array_count: int = 0
        self._int_offload_array_count: int = 0

    def _allocate_offload_arrays(self, doublesize: int, intsize: int) -> None:
        """Allocate the offload array(s).

        Parameters
        ----------
        doublesize : int
            Number of elements to allocate in the offload array.
        intsize : int
            Number of elements to allocate in the integer offload array.
        """
        if doublesize:
            if self._struct_.offload_array_size is None:
                raise ValueError("The offload array is already allocated.")
            self._offload_array_ = (ctypes.c_double * doublesize)()
        self._struct_.offload_array_size = doublesize

        if intsize:
            if self._struct_.offload_array_size is None:
                raise ValueError("The int offload array is already allocated.")
            self._int_offload_array_ = (ctypes.c_int64 * intsize)()
        self._struct_.int_offload_array_size = intsize

    def _append_offload_array(self, *data: np.ndarray) -> None:
        """Append data to the offload array.

        Parameters
        ----------
        *data : np.ndarray
            Array(s) to append.
        """
        if len(data) == 0:
            self._offload_array_count = 0
            return
        for datum in data:
            double = ctypes.sizeof(ctypes.c_double)
            if(   self._offload_array_count + datum.nbytes
                > double * self._struct_.offload_array_size ):
                raise ValueError(
                    f"The data does not fit into the offload array "
                    f"(size={self._struct_.offload_array_size}, "
                    f"counter={self._offload_array_count // double}, "
                    F"data={datum.nbytes // double})"
                    )

            pointer = ctypes.byref(self._offload_array_, self._offload_array_count)
            ctypes.memmove(pointer, datum.ctypes.data, datum.nbytes)
            self._offload_array_count += datum.nbytes

    def _slice_offload_array(self, start: int, end: int) -> np.ndarray:
        """Return a copied slice of the offload array.

        Parameters
        ----------
        start : int
            Start index of the slice (inclusive).
        end : int
            End index of the slice (non-inclusive).
        """
        return np.ctypeslib.as_array(self._offload_array_[start:end]).copy()

    def _append_int_offload_array(self, *data: np.ndarray) -> None:
        """Append data to the integer offload array.

        Parameters
        ----------
        *data : np.ndarray
            Array(s) to append.
        """
        if len(data) == 0:
            self._int_offload_array_count = 0
            return
        for datum in data:
            int32 = ctypes.sizeof(ctypes.c_int32)
            if(   self._int_offload_array_count + datum.nbytes
                > int32 * self._struct_.int_offload_array_size ):
                raise ValueError(
                    f"The data does not fit into the offload array "
                    f"(size={self._struct_.int_offload_array_size}, "
                    f"counter={self._int_offload_array_count // int32}, "
                    F"data={datum.nbytes // int32})"
                    )

            pointer = ctypes.byref(
                self._int_offload_array_, self._int_offload_array_count
                )
            ctypes.memmove(pointer, datum.ctypes.data, datum.nbytes)
            self._int_offload_array_count += datum.nbytes

    def _slice_int_offload_array(self, start: int, end: int) -> np.ndarray:
        """Return a copied slice of the integer offload array.

        Parameters
        ----------
        start : int
            Start index of the slice (inclusive).
        end : int
            End index of the slice (non-inclusive).
        """
        return np.ctypeslib.as_array(self._int_offload_array_[start:end]).copy()

    @abstractmethod
    def _export_hdf5(self) -> None:
        """Export the data to the HDF5 file."""

    @abstractmethod
    def export(self) -> Dict[str, np.ndarray]:
        """Return a dictionary with sufficient data to duplicate this instance.
        """
        return {}
