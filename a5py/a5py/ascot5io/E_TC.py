"""
Trivial cartesian electric field IO.
The trivial cartesian electric field is a constant
electric field defined in cartesian coordinates.
"""
import h5py
import numpy as np
import random
import datetime
	
from . ascot5group import creategroup

def write_hdf5(fn, Exyz):
    """
    Write trivial cartesian electric field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Exyz : real 3 x 1 numpy array
        Electric field value in cartesian coordinates
    """

    mastergroup = "efield"
    subgroup    = "E_TC"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # Check that input is a valid array with three elements.
    if Exyz.shape != (3,1) and Exyz.shape != (1,3) and Exyz.shape != (3,):
        raise Exception('Exyz has invalid format.')

    # Actual data.
    f.create_dataset(path + "/Exyz", (3,1), data = Exyz, dtype="f8")
    f.close()
    

def read_hdf5(fn, qid):
    """
    Read trivial cartesian electric field input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the efield to be read.

    Returns
    -------

    Dictionary containing electric field data.
    """

    path = "efield" + "/E_TC-" + qid

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = qid
    out["date"] = f[path].attrs["date"]
    out["description"] = f[path].attrs["description"]
    
    # Actual data.
    out["Exyz"] = f[path]["Exyz"][:]

    f.close()

    return out
