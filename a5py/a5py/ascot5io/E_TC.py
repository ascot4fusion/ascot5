"""
Trivial cartesian electric field IO.
The trivial cartesian electric field is a constant
electric field defined in cartesian coordinates.
"""
import h5py
import numpy as np
import random
import datetime
	
from . ascot5group import replacegroup, setgrouptype, setmetadata

def write_hdf5(fn, Exyz):
    """
    Write trivial cartesian electric field input in HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Exyz : real 3 x 1 numpy array
        Electric field value in cartesian coordinates
    """

    group = "efield"
    type_ = "E_TC"
    path = "efield/E_TC"

    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    setgrouptype(f, group, type_)
    replacegroup(f, path)
    setmetadata(f[path])

    # Check that input is a valid array with three elements.
    if Exyz.shape != (3,1) and Exyz.shape != (1,3) and Exyz.shape != (3,):
        raise Exception('Exyz has invalid format.')

    # Actual data.
    f.create_dataset(path + "/Exyz", data = Exyz, dtype="f8")
    f.close()
    

def read_hdf5(fn):
    """
    Read trivial cartesian electric field input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing electric field data.
    """

    path = "efield/E_TC"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

    # Actual data.
    out["Exyz"] = f[path]["Exyz"][:]

    f.close()

    return out
