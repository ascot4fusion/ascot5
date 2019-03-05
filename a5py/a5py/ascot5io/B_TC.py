"""
Trivial Cartesian magnetic field HDF5 IO.

File: B_TC.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, bxyz, jacobian, rhoval, psival=rhoval, axisr=1, axisz=0,
               desc=None):
    """
    Write trivial cartesian magnetic field input to HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        bxyz : array_like (3,1) <br>
            Magnetic field in cartesian coordinates at origo.
        jacobian : array_like (3,3) <br>
            Magnetic field Jacobian, jacobian[i,j] = dB_i/dx_j
        rhoval: float <br>
            Constant rho value.
        psival: float, optional <br>
            Constant psi value.
        axisr: float, optional <br>
            Magnetic axis R coordinate.
        axisz: real, optional <br>
            Magnetic axis z coordinate.
        desc : str, optional <br>
            Input description.

    Returns:
        QID of the new input that was written.
    """

    parent = "bfield"
    group  = "B_TC"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("bxyz",     (3,1), data=bxyz,     dtype="f8")
        g.create_dataset("jacobian", (3,3), data=jacobian, dtype="f8")
        g.create_dataset("rhoval",   (1,),  data=rhoval,   dtype="f8")
        g.create_dataset("psival",   (1,),  data=psival,   dtype="f8")
        g.create_dataset("axisr",    (1,),  data=axisR,    dtype="f8")
        g.create_dataset("axisz",    (1,),  data=axisz,    dtype="f8")

    return g.name


def read_hdf5(fn, qid):
    """
    Read Cartesian magnetic field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "bfield/B_TC_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


class B_TC(AscotData):
    """
    Object representing B_TC data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
