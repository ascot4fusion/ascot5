"""
Options IO.
"""
import h5py
import numpy as np
import random
import datetime

def write_hdf5(fn,options):
    """
    Write options.

    Unlike most other "write" functions, this one takes dictionary
    as an argument. The dictionary should have exactly the same format
    as given by the "read" function in this module.

    TODO not compatible with new format

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    options : dictionary
        Options to be written in dictionary format.
    """

    group = "options"
    path = "options/"
        
    f = h5py.File(fn, "a")

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
    for opt in options:
        f.create_dataset(path + opt, data=f[path][opt])
    
    f.close()


def read_hdf5(fn):
    """
    Read options from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing options.
    """
    
    path = "options"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

    # Actual data.
    for opt in f[path]:
        out[opt] = f[path][opt][:]

    f.close()

    return out
