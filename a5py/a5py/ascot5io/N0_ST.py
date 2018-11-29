"""
Stellarator neutral density HDF5 IO
"""
import numpy as np
import h5py
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz,
               phimin, phimax, nphi, n_periods, n0):
    """
    Write stellarator neutral density input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rmin, Rmax, phimin, phimax, zmin, zmax : real
        Edges of the uniform Rphiz-grid.
    nR, nphi, nz : int
        Number of Rphiz-grid points.
    n_periods : int
        Number of toroidal periods in the device.
    n0 : real R x phi x z numpy array
        Neutral density in Rphiz-grid for half a period.

    """

    mastergroup = "neutral"
    subgroup    = "N0_ST"

    # Transpose n0 from (r, phi, z) to (z, phi, r)
    n0 = np.transpose(n0,(2,1,0))

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = creategroup(f, mastergroup, subgroup)

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

        # Toroidal periods
        f.create_dataset(path + "/toroidalPeriods", data=n_periods, dtype="i4")

        # Neutral data
        f.create_dataset(path + "/n0",   data=n0, dtype="f8")


def read_hdf5(fn):
    """
    Read stellarator neutral density input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing neutral density data.
    """

    path = "neutral/N0_ST"

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = f[path].attrs["qid"]
        out["date"] = f[path].attrs["date"]

        # Actual data.
        out["Rmin"] = f[path]["R_min"][:]
        out["Rmax"] = f[path]["R_max"][:]
        out["nR"]   = f[path]["n_R"][:]

        out["phimin"] = f[path]["phi_min"][:]
        out["phimax"] = f[path]["phi_max"][:]
        out["nphi"]   = f[path]["n_phi"][:]

        out["zmin"] = f[path]["z_min"][:]
        out["zmax"] = f[path]["z_max"][:]
        out["nz"]   = f[path]["n_z"][:]

        out["n0"]   = f[path]["n0"]

        out["n_periods"] = f[path]["toroidalPeriods"][:]

    return out
