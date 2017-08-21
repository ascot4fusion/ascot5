"""
2D wall IO.
"""
import h5py
import numpy as np
import random
import datetime

def write_hdf5(fn, n, x1x2x3, y1y2y3, z1z2z3, flag):
    """
    Write 3D wall input in HDF5 file.

    TODO Not compatible with new HDF5 format.

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

    group = "wall"
    type_ = "3D"
    path = "wall/3D"

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

    # Actual data.
    f.create_dataset('wall/3D/x1x2x3', (n,3), dtype='f8', data=x1x2x3)
    f.create_dataset('wall/3D/y1y2y3', (n,3), dtype='f8', data=y1y2y3)
    f.create_dataset('wall/3D/z1z2z3', (n,3), dtype='f8', data=z1z2z3)
    f.create_dataset('wall/3D/flag', (n,1), dtype='i4', data=flag)
    f['wall/3D'].attrs['n_elements'] = n
    f['wall/3D'].attrs['min_x'] = np.amin(x1x2x3)
    f['wall/3D'].attrs['max_x'] = np.amax(x1x2x3)
    f['wall/3D'].attrs['min_y'] = np.amin(y1y2y3)
    f['wall/3D'].attrs['max_y'] = np.amax(y1y2y3)
    f['wall/3D'].attrs['min_z'] = np.amin(z1z2z3)
    f['wall/3D'].attrs['max_z'] = np.amax(z1z2z3)

    f.close()


def read_hdf5(fn):
    """
    Read 3D wall input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing wall data.
    """

    path = "wall/3D"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

    # Actual data.
    out["n"]      = f[path]["n"][:]
    out["x1x2x3"] = f[path]["x1x2x3"][:]
    out["y1y2y3"] = f[path]["y1y2y3"][:]
    out["z1z2z3"] = f[path]["z1z2z3"][:]
    out["flag"]   = f[path]["flag"][:]
    
    f.close()
