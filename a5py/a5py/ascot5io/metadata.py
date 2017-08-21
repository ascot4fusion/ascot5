"""
Metadata IO.
"""
import numpy as np
import h5py

def write_hdf5(fn, meta):
    """
    Write metadata in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    meta : str
        Unspecified meta information as a single string.
    """
    path = "metadata/"

    f = h5py.File(fn,"a")
        
    # Remove group if one is already present.
    if path in f:
        del f[path]
    f.create_group(path)

    f.create_dataset(path + "meta", data=meta, dtype="S10")

    f.close()

def read_hdf5(fn):
    """
    Read metadata from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing metadata.
    """
    
    path = "metadata/"

    f = h5py.File(fn,"r")

    out = {}
    out["meta"] = f[path + "meta"][:]

    f.close()


