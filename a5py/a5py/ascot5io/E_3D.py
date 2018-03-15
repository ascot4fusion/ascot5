"""
Non-axisymmetric tokamak magnetic field HDF5 IO
"""
import numpy as np
import h5py
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               E_R, E_phi, E_z):
    """
    Write 3D electric field input in HDF5 file for trilinear interpolation.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rlim, Rmax, phimin, phimax, zmin, zmax : real
        Edges of the uniform Rphiz-grid.
    nR, nphi, nz : int
        Number of Rphiz-grid points.
    E_R, E_phi, E_z : real R x phi x z numpy array
        Electric field components in Rphiz-grid.

    Notes
    -------
    
    Rphiz coordinates form a right-handed cylindrical coordinates.
    phi is in degrees and is considered as a periodic coordinate.
    phimin is where the period begins and phimax is the last data point,
    i.e. not the end of period. E.g if you have a n = 2 symmetric system
    with nphi = 180 deg and phimin = 0 deg, then phimax should be 179 deg.

    """

    mastergroup = "efield"
    subgroup    = "E_3D"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # Transpose grids
    E_R = np.transpose(E_R,(1,0,2))
    E_phi = np.transpose(E_phi,(1,0,2))
    E_z = np.transpose(E_z,(1,0,2))
    
    # TODO Check that inputs are consistent.
    
    # Actual data.
    f.create_dataset(path + "/R_min", (1,), data=Rmin, dtype="f8")
    f.create_dataset(path + "/R_max", (1,), data=Rmax, dtype="f8")
    f.create_dataset(path + "/n_R", (1,),   data=nR, dtype="i8")

    f.create_dataset(path + "/phi_min", (1,), data=phimin, dtype="f8")
    f.create_dataset(path + "/phi_max", (1,), data=phimax, dtype="f8")
    f.create_dataset(path + "/n_phi", (1,),   data=nphi, dtype="i8")

    f.create_dataset(path + "/z_min", (1,), data=zmin, dtype="f8")
    f.create_dataset(path + "/z_max", (1,), data=zmax, dtype="f8")
    f.create_dataset(path + "/n_z", (1,),   data=nz, dtype="i8")

    f.create_dataset(path + "/E_R",   data=E_R,   dtype="f8")
    f.create_dataset(path + "/E_phi", data=E_phi, dtype="f8")
    f.create_dataset(path + "/E_z",   data=E_z,   dtype="f8")

    f.close()


def read_hdf5(fn):
    """
    Read 3D electric field input from HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing magnetic field data.
    """

    path = "efield/E_3D"

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

    out["E_R"]   = f[path]["E_r"][:]
    out["E_phi"] = np.reshape(f[path]["E_phi"][:], (out["nz"][0], out["nphi"][0], out["nR"][0]))#f[path]["B_phi"][:]
    out["E_z"]   = f[path]["E_z"][:]
    
    f.close()

    return out
