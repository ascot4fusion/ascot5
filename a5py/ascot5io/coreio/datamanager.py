"""Abstract classes for objects containing data.
"""
from __future__ import annotations

from typing import Any, List

import h5py
import numpy as np
from enum import Enum

from a5py.exceptions import AscotIOException

class DataManager():
    """A class to abstract the data storage.

    This class provides a unified interface to store and access data via same
    references, regardless of whether the actual data is stored in HDF5 files,
    IMAS IDS, or in memory.

    For example, `get("a")` will retrieve the value of "a" using a method
    based on this object's `_format`.
    """

    class Format(Enum):
        """Data storage formats."""
        HDF5 = 1
        MEMORY = 2
        IMAS_IDS = 3

    HDF5 = Format.HDF5
    MEMORY = Format.MEMORY
    IMAS_IDS = Format.IMAS_IDS

    def __init__(self,
                 hdf5_filename: str=None,
                 path_within_hdf5: str=None,
                 imas_ids: str=None,
                 **kwargs: Any,
                 ) -> None:
        """Initializes an object which by default stores the data in memory.

        Parameters
        ----------
        hdf5_filename : str, optional
            Name of the HDF5 file to store the data.
        path_within_hdf5 : str, optional
            The path of the group within the HDF5 that has the datasets.
        **kwargs
            Key value pairs for data if the data is being stored in the memory.
        """
        super().__init__()
        self._imas_ids: str = imas_ids
        self._hdf5_filename: str = hdf5_filename
        self._path_within_hdf5: str = path_within_hdf5

        self._format: DataManager.Format = DataManager.MEMORY

        if imas_ids:
            self._format = DataManager.IMAS_IDS
        elif hdf5_filename and path_within_hdf5:
            self._format = DataManager.HDF5

        if self._format == DataManager.MEMORY:
            self._storage = kwargs

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return (
            f"{self.__class__.__name__}(hdf5_filename={self._hdf5_filename}, "
            f"path_within_hdf5={self._path_within_hdf5}, "
            f"imas_ids={self._imas_ids})"
        )

    def get(self, key: str) -> np.array:
        """Get the value of the given attribute.

        The value is read from wherever the data is stored.
        """
        if self._format == DataManager.MEMORY:
            return self._storage[key]
        elif self._format == DataManager.IMAS_IDS:
            return self._imas_ids.get(key)
        elif self._format == DataManager.HDF5:
            with h5py.File(self._hdf5_filename, "r") as f:
                return f[self._path_within_hdf5][key][:]

    def migrate(self,
                newformat: DataManager.Format,
                hdf5_filename: str=None,
                path_within_hdf5: str=None,
                imas_ids: str=None,
                ):
        """Change the storage locationand write the existing data to the new
        location.
        """
        if newformat == self._format:
            return
        if newformat == DataManager.MEMORY:
            pass
        elif newformat == DataManager.HDF5:
            self._hdf5_filename = hdf5_filename
            self._path_within_hdf5 = path_within_hdf5
            with h5py.File(self._hdf5_filename, "r") as f:
                for key, data in self._storage:
                    f[self._path_within_hdf5][key] = data.value
        else:
            raise AscotIOException(f"Unsupported format: {newformat}")

class DataContainer():
    pass


class DataGroup(DataContainer):
    pass
