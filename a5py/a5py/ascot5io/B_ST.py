"""
Stellarator magnetic field IO.
"""
import numpy as np
import h5py
import random
import datetime

from . ascot5group import creategroup

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               B_R, B_phi, B_z, psi, n_periods,
               axismin, axismax, naxis, axisR, axisz,
               pRmin=None, pRmax=None, pnR=None, pzmin=None, pzmax=None, pnz=None,
               symmetry_mode=0, psiaxis=0, psisepx=1):
    """
    Write stellarator magnetic field input in HDF5 file.

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
    psi : real
        Normalized toroidal flux in Rphiz-grid
    n_periods : int
        Number of toroidal periods.
    naxis : int
        Number of axis grid points.
    axisR, axisz : real
        Magnetic axis R- and z-location as a function of phi.
    pRmin, pRmax, pnR, pphimin, pphimax, pnphi, pzmin, pzmax, pnz : opt
        Optional parameters that define a separate grid for psi.
    symmetry_mode : opt
        Mode of symmetry used. 0 = stellarator symmetry, 1 = toroidal periodic
    psiaxis, psisepx : real    
        Psi values at magnetic axis and separatrix
    """

    mastergroup = "bfield"
    subgroup    = "B_STS"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)

    # TODO Check that inputs are consistent.

    # Define psigrid to be same as Bgrid if not stated otherwise.
    if(pRmin is None or pRmax is None or pnR is None
       or pphimin is None or pphimax is None or pnphi is None 
       or pzmin is None or pzmax is None or pnz is None):
        pRmin   = Rmin
        pRmax   = Rmax
        pnR     = nR
        pphimin = phimin
        pphimax = phimax
        pnphi   = nphi
        pzmin   = zmin
        pzmax   = zmax
        pnz     = nz

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

    f.create_dataset(path + "/psigrid_R_min", (1,), data=pRmin, dtype="f8")
    f.create_dataset(path + "/psigrid_R_max", (1,), data=pRmax, dtype="f8")
    f.create_dataset(path + "/psigrid_n_R", (1,),   data=pnR, dtype="i8")
    
    f.create_dataset(path + "/psigrid_phi_min", (1,), data=pphimin, dtype="f8")
    f.create_dataset(path + "/psigrid_phi_max", (1,), data=pphimax, dtype="f8")
    f.create_dataset(path + "/psigrid_n_phi", (1,),   data=pnphi, dtype="i8")

    f.create_dataset(path + "/psigrid_z_min", (1,), data=pzmin, dtype="f8")
    f.create_dataset(path + "/psigrid_z_max", (1,), data=pzmax, dtype="f8")
    f.create_dataset(path + "/psigrid_n_z", (1,),   data=pnz, dtype="i8")

    # Magnetic field data
    f.create_dataset(path + "/B_r",   data=B_R, dtype="f8")
    f.create_dataset(path + "/B_phi", data=B_phi, dtype="f8")
    f.create_dataset(path + "/B_z",   data=B_z, dtype="f8")
    f.create_dataset(path + "/psi",   data=psi, dtype="f8")

    # Magnetic axis
    f.create_dataset(path + "/axis_min", (1,), data=axismin, dtype="f8")
    f.create_dataset(path + "/axis_max", (1,), data=axismax, dtype="f8")
    f.create_dataset(path + "/n_axis", (1,),   data=naxis, dtype="i8")

    f.create_dataset(path + "/axis_r",   data=axisR, dtype="f8")
    f.create_dataset(path + "/axis_z",   data=axisz, dtype="f8")

    f.create_dataset(path + "/psi0", (1,), data=psiaxis, dtype="f8")
    f.create_dataset(path + "/psi1", (1,), data=psisepx, dtype="f8")

    # Toroidal periods
    f.create_dataset(path + "/toroidalPeriods", (1,), data=n_periods, dtype="i4")

    # Symmetry mode
    f.create_dataset(path + "/symmetry_mode", (1,), data=symmetry_mode, dtype="i4")

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

    out["psi"]   = f[path]["psi"][:]
    out["B_R"]   = f[path]["B_r"][:]
    out["B_phi"] = f[path]["B_phi"][:]
    out["B_z"]   = f[path]["B_z"][:]

    out["psi0"] = f[path]["psi0"][:]
    out["psi1"] = f[path]["psi1"][:]

    out["axisr"] = f[path]["axis_r"][:]
    out["axisz"] = f[path]["axis_z"][:]

    out["axismin"] = f[path]["axis_min"][:]
    out["axismax"] = f[path]["axis_max"][:]
    out["naxis"]   = f[path]["n_axis"][:]

    out["n_periods"] = f[path]["toroidalPeriods"][:]

    f.close()

    return out
