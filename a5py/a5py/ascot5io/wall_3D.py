"""
3D wall IO.

File: wall_3D.py
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5file import add_group
from a5py.ascot5io.base import AscotInput

def write_hdf5(fn, n, x1x2x3, y1y2y3, z1z2z3, flag, desc=None):
    """
    Write 3D wall input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    n : int
        Number of wall triangles
    x1x2x3, y1y2y3, z1z2z3 : real n x 3 numpy array
        Wall triangle vertices' xyz coordinates.
    flag : int n x 1 numpy array
        Indicates which part of the wall (e.g. divertor) triangle belongs to.
    """

    parent = "wall"
    group  = "wall_3D"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group, desc=desc)

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + '/x1x2x3', (n,3), dtype='f8', data=x1x2x3)
        f.create_dataset(path + '/y1y2y3', (n,3), dtype='f8', data=y1y2y3)
        f.create_dataset(path + '/z1z2z3', (n,3), dtype='f8', data=z1z2z3)
        f.create_dataset(path + '/flag',  (n,1), dtype='i4', data=flag)
        f.create_dataset(path + '/n',     (1,), dtype='i4', data=n)
        f.create_dataset(path + '/min_x', (1,), dtype='f8', data=np.amin(x1x2x3))
        f.create_dataset(path + '/max_x', (1,), dtype='f8', data=np.amax(x1x2x3))
        f.create_dataset(path + '/min_y', (1,), dtype='f8', data=np.amin(y1y2y3))
        f.create_dataset(path + '/max_y', (1,), dtype='f8', data=np.amax(y1y2y3))
        f.create_dataset(path + '/min_z', (1,), dtype='f8', data=np.amin(z1z2z3))
        f.create_dataset(path + '/max_z', (1,), dtype='f8', data=np.amax(z1z2z3))


def read_hdf5(fn, qid):
    """
    Read 3D wall input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the wall to be read.

    Returns
    -------

    Dictionary containing wall data.
    """

    path = "wall" + "/wall_3D-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["n"]      = f[path]["n"][:]
        out["x1x2x3"] = f[path]["x1x2x3"][:]
        out["y1y2y3"] = f[path]["y1y2y3"][:]
        out["z1z2z3"] = f[path]["z1z2z3"][:]
        out["flag"]   = f[path]["flag"][:]

    return out

class wall_3D(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
