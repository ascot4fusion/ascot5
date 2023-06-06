"""
Trivial cartesian electric field IO.

File: E_TC.py
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from . E import E

def write_hdf5(fn, exyz, desc=None):
    """
    Write trivial cartesian electric field input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        exyz : array_like (3,1) <br>
            Electric field value in cartesian coordinates
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "efield"
    group  = "E_TC"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("exyz", (3,1), data=exyz, dtype="f8")

    return gname


def read_hdf5(fn, qid):
    """
    Read Cartesian electric field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "efield/E_TC_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


def write_hdf5_dummy(fn, desc="Dummy"):
    """
    Write dummy data.
    """
    return write_hdf5(fn, np.array([0,0,0]), desc=desc)


class E_TC(E):
    """
    Object representing E_TC data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
