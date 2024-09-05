"""Common interface to access data independently from how it is stored."""
from __future__ import annotations

from enum import Enum
from typing import Type, TypeVar
from abc import ABC, abstractmethod

import ctypes
import numpy as np

T = TypeVar('T', bound='DataHolder')


class Format(Enum):
    """Data storage formats."""
    HDF5 = 1
    """Data is stored in HDF5 file."""
    CSTRUCT = 2
    """Data is stored in memory using C struct and offload arrays."""


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

    def __init__(self, struct, *args, **kwargs) -> None:
        """Initialize an object that does not yet contain any data.

        Parameters
        ----------
        struct : ctypes.Structure
            The struct defining the data stored in C.
        """
        super().__init__(*args, **kwargs)
        self._format = Format.CSTRUCT
        self._struct_ = struct
        self._offload_array_ = None
        self._int_offload_array_ = None
        self._offload_array_count = 0
        self._int_offload_array_count = 0

    def _append_offload_array(self, *data):
        """Append data to the offload array."""
        if not len(data):
            self._offload_array_count = 0
            return
        for datum in data:
            DOUBLE = ctypes.sizeof(ctypes.c_double)
            if(   self._offload_array_count + datum.nbytes
                > DOUBLE * self._struct_.offload_array_size ):
                raise ValueError(
                    f"The data does not fit into the offload array "
                    f"(size={self._struct_.offload_array_size}, "
                    f"counter={self._offload_array_count // DOUBLE}, "
                    F"data={datum.nbytes // DOUBLE})"
                    )

            pointer = ctypes.byref(self._offload_array_, self._offload_array_count)
            ctypes.memmove(pointer, datum.ctypes.data, datum.nbytes)
            self._offload_array_count += datum.nbytes

    def _slice_offload_array(self, start, end):
        """Return a copied slice of the offload array."""
        return np.ctypeslib.as_array(self._offload_array_[start:end]).copy()

    def _append_int_offload_array(self, *data):
        """Append data to the integer offload array."""
        if not len(data):
            self._int_offload_array_count = 0
            return
        for datum in data:
            INT32 = ctypes.sizeof(ctypes.c_int32)
            if(   self._int_offload_array_count + datum.nbytes
                > INT32 * self._struct_.int_offload_array_size ):
                raise ValueError(
                    f"The data does not fit into the offload array "
                    f"(size={self._struct_.int_offload_array_size}, "
                    f"counter={self._int_offload_array_count // INT32}, "
                    F"data={datum.nbytes // INT32})"
                    )

            pointer = ctypes.byref(
                self._int_offload_array_, self._int_offload_array_count
                )
            ctypes.memmove(pointer, datum.ctypes.data, datum.nbytes)
            self._int_offload_array_count += datum.nbytes

    def _slice_int_offload_array(self, start, end):
        """Return a copied slice of the integer offload array."""
        return np.ctypeslib.as_array(self._int_offload_array_[start:end]).copy()

    @abstractmethod
    def _export_hdf5(self, *args, **kwargs):
        """Export the data to the HDF5 file."""
        pass

    @abstractmethod
    def export(cls, *args, **kwargs):
        """Export the data to `nd.array`s."""
        return {}

    @classmethod
    def init_hdf5(cls: Type[T]) -> T:
        """Initialize an instance from HDF5 file."""
        obj = cls()
        obj._format = Format.HDF5
        return obj
