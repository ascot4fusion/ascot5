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

    with h5py.File(fn,"a") as f:
        # Remove group if one is already present.
        if path in f:
            del f[path]
        f.create_group(path)

        f.create_dataset(path + "meta", data=meta, dtype="S10")


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

    with h5py.File(fn,"r") as f:
        out = {}
        out["meta"] = f[path + "meta"][:]

    return out
