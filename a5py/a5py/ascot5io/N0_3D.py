"""
Non-axisymmetric tokamak neutral density HDF5 IO
"""
import numpy as np
import h5py
import random
import datetime

from . ascot5group import replacegroup, setgrouptype, setmetadata

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi, n0):
    """
    Write 3D neutral density input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rmin, Rmax, phimin, phimax, zmin, zmax : real
        Edges of the uniform Rphiz-grid.
    nR, nphi, nz : int
        Number of Rphiz-grid points.
    n0 : real R x phi x z numpy array
        Neutral density in Rphiz-grid.

    Notes
    -------
    
    Rphiz coordinates form a right-handed cylindrical coordinates.
    phi is in degrees and is considered as a periodic coordinate.
    phimin is where the period begins and phimax is the last data point,
    i.e. not the end of period. E.g if you have a n = 2 symmetric system
    with nphi = 180 deg and phimin = 0 deg, then phimax should be 179 deg.

    """

    group = "neutral"
    type_ = "N0_3D"
    path = "neutral/N0_3D"
    
    # Create group and set the type to this one.
    f = h5py.File(fn, "a")
    setgrouptype(f, group, type_)
    replacegroup(f, path)
    setmetadata(f[path])

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

    f.create_dataset(path + "/n0",   data=n0, dtype="f8")
    f.close()


def read_hdf5(fn):
    """
    Read 3D neutral density input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing neutral density data.
    """

    path = "neutral/N0_3D"

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = f[path].attrs["qid"]
    out["date"] = f[path].attrs["date"]

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

    out["n0"]   = f[path]["n0"]

    f.close()

    return out
