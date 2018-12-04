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

from . ascot5group import creategroup, setdescription

def write_hdf5(fn, Bxyz, J, rhoval, psival=0, axisR=1, axisz=0, desc=None):
    """
    Write trivial cartesian magnetic field input in HDF5 file.

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
    psival: real, optional
        Constant psi value. Same as rho by default.
    axisR: real, optional
        Magnetic axis R coordinate. Default value 1.
    axisz: real, optional
        Magnetic axis z coordinate. Default value 0.
    """

    mastergroup = "bfield"
    subgroup    = "B_TC"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = creategroup(f, mastergroup, subgroup, desc=desc)

        # TODO Check that inputs are consistent.
        if psival == 0:
            psival = rhoval

        # Actual data.
        f.create_dataset(path + "/Bxyz",         data=Bxyz,   dtype="f8")
        f.create_dataset(path + "/J",            data=J,      dtype="f8")
        f.create_dataset(path + "/rhoval", (1,), data=rhoval, dtype="f8")
        f.create_dataset(path + "/psival", (1,), data=psival, dtype="f8")
        f.create_dataset(path + "/axisr", (1,),  data=axisR,  dtype="f8")
        f.create_dataset(path + "/axisz", (1,),  data=axisz,  dtype="f8")


def read_hdf5(fn, qid):
    """
    Read trivial cartesian magnetic field input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the bfield to be read.

    Returns
    -------

    Dictionary containing magnetic field data.
    """

    path = "bfield" + "/B_TC-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["Bxyz"]   = f[path]["Bxyz"][:]
        out["J"]      = f[path]["J"][:]
        out["rhoval"] = f[path]["rhoval"][:]
        out["psival"] = f[path]["psival"][:]
        out["axisR"]  = f[path]["axisr"][:]
        out["axisz"]  = f[path]["axisz"][:]

    return out
