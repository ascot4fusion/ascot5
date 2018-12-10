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
        self._isopen = False

    def get_desc(self):
        f = self._open()
        val = ascot5file.get_desc(f.file, self._group)
        self._close(f)
        return val

    def get_date(self):
        f = self._open()
        val = ascot5file.get_date(f.file, self._group)
        self._close(f)
        return val

    def get_qid(self):
        return ascot5file.get_qid(self._group)

    def get_type(self):
        return self._path.split("/")[-1].split("-")[0]

    def _open(self):
        if self._isopen == True:
            return None

        self._isopen = True
        return h5py.File(self._file, "r")[self._path]

    def _close(self, f):
        self._isopen = False
        f.file.close()
