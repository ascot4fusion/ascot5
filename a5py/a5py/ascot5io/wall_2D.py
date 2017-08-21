"""
2D wall IO.
"""
import h5py
import numpy as np
import random
import datetime

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
    if not group in f:
        o = f.create_group(group)
        o.attrs["type"] = np.string_(type_)
    else:
        o = f[group]
        del o.attrs["type"]
        o.attrs["type"] = np.string_(type_)

    # Remove group if one is already present.
    if path in f:
        del f[path]
    f.create_group(path)

    # TODO Check that inputs are consistent.

    # Metadata.
    qid = random.getrandbits(64)
    f[path].attrs["qid"]  = np.int64_(qid)
    f[path].attrs["date"] = np.string_(datetime.datetime.now())


    f.create_dataset(path + "/n", data=n)
    f.create_dataset(path + "/r", data=R)
    f.create_dataset(path + "/z", data=z)

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
