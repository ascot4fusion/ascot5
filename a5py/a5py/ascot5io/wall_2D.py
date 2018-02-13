"""
2D wall IO.
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, n, R, z):
    """
    Write 2D wall input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    n : int
        Number of wall segments
    R, z : real n x 1 numpy array
        Wall segment vertices' R and z coordinates

    Notes
    -----
    
    First (R,z) point doesn't have to correspond to last point
    as the wall is closed automatically. TODO check this
    """

    mastergroup = "wall"
    subgroup    = "wall_2D"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # TODO Check that inputs are consistent.

    f.create_dataset(path + "/n", (1,), data=n, dtype='i4')
    f.create_dataset(path + "/r", data=R, dtype='f8')
    f.create_dataset(path + "/z", data=z, dtype='f8')

    f.close()


def read_hdf5(fn):
    """
    Read 2D wall input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing wall data.
    """

    path = "wall/2D"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

    # Actual data.
    out["n"] = f[path]["n"][:]
    out["R"] = f[path]["r"][:]
    out["z"] = f[path]["z"][:]
    
    f.close()

    return out
