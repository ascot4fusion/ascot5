"""
Stellarator magnetic field IO.
"""
import numpy as np
import h5py
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               B_R, B_phi, B_z, s, n_periods,
               axisR, axisphi, axisz):
    """
    Write stellarator magnetic field input in HDF5 file.

    TODO fill documentation

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rlim, Rmax, phimin, phimax, zmin, zmax : real
        Edges of the uniform Rphiz-grid.
    nR, nphi, nz : int
        Number of Rphiz-grid points.
    B_R, B_phi, B_z : real R x phi x z numpy array
        Magnetic field components in Rphiz-grid
    s : real
        TODO
    n_periods : int
        TODO
    axisR, axisphi, axisz : real
        Magnetic axis Rphiz-location.
    """

    mastergroup = "bfield"
    subgroup    = "B_STS"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # TODO Check that inputs are consistent.

    # Actual data.
    f.create_dataset(path + "/r_min", (1,), data=Rmin, dtype="f8")
    f.create_dataset(path + "/r_max", (1,), data=Rmax, dtype="f8")
    f.create_dataset(path + "/n_r", (1,),   data=nR, dtype="i8")

    f.create_dataset(path + "/phi_min", (1,), data=phimin, dtype="f8")
    f.create_dataset(path + "/phi_max", (1,), data=phimax, dtype="f8")
    f.create_dataset(path + "/n_phi", (1,),   data=nphi, dtype="i8")

    f.create_dataset(path + "/z_min", (1,), data=zmin, dtype="f8")
    f.create_dataset(path + "/z_max", (1,), data=zmax, dtype="f8")
    f.create_dataset(path + "/n_z", (1,),   data=nz, dtype="i8")

    # Magnetic field data
    f.create_dataset(path + "/B_r",   data=B_R, dtype="f8")
    f.create_dataset(path + "/B_phi", data=B_phi, dtype="f8")
    f.create_dataset(path + "/B_z",   data=B_z, dtype="f8")
    f.create_dataset(path + "/s",     data=s, dtype="f8")

    # Magnetic axis (TODO not implemented)
    f.create_dataset(path + "/axis_r",   data=axisR, dtype="f8")
    f.create_dataset(path + "/axis_phi", data=axisphi, dtype="f8")
    f.create_dataset(path + "/axis_z",   data=axisz, dtype="f8")

    # Toroidal periods
    f.create_dataset(path + "/toroidalPeriods", (1,), data=n_periods, dtype="i4")
    
    f.close()


def read_hdf5(fn, qid):
    """
    Read stellarator magnetic field input from HDF5 file.

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

    path = "bfield" + "/B_STS-" + qid

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = qid
    out["date"] = f[path].attrs["date"]
    out["description"] = f[path].attrs["description"]

    # Actual data.
    out["Rmin"] = f[path]["r_min"][:]
    out["Rmax"] = f[path]["r_max"][:]
    out["nR"]   = f[path]["n_r"][:]

    out["phimin"] = f[path]["phi_min"][:]
    out["phimax"] = f[path]["phi_max"][:]
    out["nphi"]   = f[path]["n_phi"][:]

    out["zmin"] = f[path]["z_min"][:]
    out["zmax"] = f[path]["z_max"][:]
    out["nz"]   = f[path]["n_z"][:]

    out["s"]     = f[path]["s"][:]
    out["B_R"]   = f[path]["B_r"][:]
    out["B_phi"] = f[path]["B_phi"][:]
    out["B_z"]   = f[path]["B_z"][:]

    out["axisR"]   = f[path]["axis_r"][:]
    out["axisphi"] = f[path]["axis_phi"][:]
    out["axisz"]   = f[path]["axis_z"][:]

    out["n_periods"] = f[path]["toroidalPeriods"][:]

    f.close()

    return out
