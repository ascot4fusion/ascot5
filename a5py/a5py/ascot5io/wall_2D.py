"""
2D wall IO.
"""
import h5py
import numpy as np
import random
import datetime

from . ascot5group import replacegroup, setgrouptype, setmetadata

def write_hdf5(fn, n, R, z):
    """
    Write 2D wall input in HDF5 file.

    TODO Not compatible with new HDF5 format.

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

    group = "wall"
    type_ = "2D"
    path = "wall/2D"

    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    setgrouptype(f, group, type_)
    replacegroup(f, path)
    setmetadata(f[path])

    # TODO Check that inputs are consistent.

    f[path].attrs['n'] = n
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
