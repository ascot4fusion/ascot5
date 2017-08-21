"""
Trivial Cartesian magnetic field HDF5 IO.
TC is magnetic field with constant Jacobian matrix. The field is
defined in terms of cartesian coordinates so it does not support
all ASCOT5 features, and is only intended for debugging.
"""
import h5py
import numpy as np
import random
import datetime

def write_hdf5(fn, Bxyz, J, rhoval):
    """
    Write trivial cartesian magnetic field input in HDF5 file.

    TODO Not compatible with new HDF5 format. (Or old...)

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Bxyz : real 3 x 1 numpy array
        Magnetic field in cartesian coordinates at origo.
    J : real 9 x 9 numpy array
        Magnetic field Jacobian, J(i,j) = dB_i/dx_j
    rhoval: real
        Constant rho (simultaneously psi) value.
    """
    
    group = "bfield"
    type_ = "B_TC"
    path = "bfield/B_TC"
    
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
    f.create_dataset(path + "/Bxyz",   data=Bxyz, dtype="f8")
    f.create_dataset(path + "/J",      data=J, dtype="f8")
    f.create_dataset(path + "/rhoval", data=rhoval, dtype="f8")
    f.close()


def read_hdf5(fn)
    """
    Read trivial cartesian magnetic field input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing magnetic field data.
    """

    path = "bfield/B_TC"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

    # Actual data.
    out["Bxyz"]   = f[path]["Bxyz"][:]
    out["J"]      = f[path]["J"][:]
    out["rhoval"] = f[path]["rhoval"][:]

    f.close()
