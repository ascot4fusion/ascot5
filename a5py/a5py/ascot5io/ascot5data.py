"""
Parent class for input and output classes and common methods for those.

File: ascot5data.py
"""

import h5py

from . import ascot5file

class AscotData():
    """
    Abstract class for ASCOT5 input and output data.
    """

    def __init__(self, hdf5):
        self._file   = hdf5.file.filename
        self._group  = hdf5.name.split("/")[2]
        self._path   = hdf5.name
        self._opened = None

    def __enter__(self):
        return self._open()

    def __exit__(self, exception_type, exception_value, traceback):
        self._close()

    def get_desc(self):
        f = self._open()
        val = ascot5file.get_desc(f.file, self._group)
        self._close()
        return val

    def get_date(self):
        f = self._open()
        val = ascot5file.get_date(f.file, self._group)
        self._close()
        return val

    def get_qid(self):
        return ascot5file.get_qid(self._group)

    def get_type(self):
        return self._path.split("/")[-1].split("_")[0]

    def get_name(self):
        return self._group

    def _open(self):
        if self._opened is not None:
            return None

        self._opened = h5py.File(self._file, "a")[self._path]
        return self._opened

    def _close(self):
        self._opened.file.close()
        self._opened = None
