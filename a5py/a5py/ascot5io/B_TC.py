"""
Trivial Cartesian magnetic field HDF5 IO.

TC is magnetic field with constant Jacobian matrix. The field is
defined in terms of cartesian coordinates so it does not support
all ASCOT5 features, and is only intended for debugging.

File: B_TC.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, Bxyz, J, rhoval, psival=None, axisR=1, axisz=0, desc=None):
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

    parent = "bfield"
    group  = "B_TC"

    # TODO Check that inputs are consistent.
    if psival == None:
        psival = rhoval

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        # Actual data.
        g.create_dataset("bxyz",     (3,1), data=Bxyz,   dtype="f8")
        g.create_dataset("jacobian", (3,3), data=J,      dtype="f8")
        g.create_dataset("rhoval",   (1,),  data=rhoval, dtype="f8")
        g.create_dataset("psival",   (1,),  data=psival, dtype="f8")
        g.create_dataset("axisr",    (1,),  data=axisR,  dtype="f8")
        g.create_dataset("axisz",    (1,),  data=axisz,  dtype="f8")


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

class B_TC(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
