"""
Non-axisymmetric tokamak magnetic field HDF5 IO

File: B_3DS.py
"""
import numpy as np
import h5py
import random
import datetime

from . ascot5file import add_group
from a5py.ascot5io.base import AscotInput

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               axisR, axisz, psiRz, psiaxis, psisepx,
               B_R, B_phi, B_z,
               pRmin=None, pRmax=None, pnR=None, pzmin=None, pzmax=None,
               pnz=None, desc=None):
    """
    Write 3DS magnetic field input in HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rlim, Rmax, phimin, phimax, zmin, zmax : real
        Edges of the uniform R x phi x z-grid [m x deg x m].
    nR, nphi, nz : int
        Number of Rphiz-grid points.
    axisR, axisz : real
        Magnetic axis Rz-location.
    psiRz : real R x z numpy array
        Psi values in the Rz-grid.
    psiaxis, psisepx : real
        Psi values at magnetic axis and separatrix
    B_R, B_phi, B_z : real R x phi x z numpy array
        Magnetic field components in Rphiz-grid.
    pRmin, pRmax, pnR, pzmin, pzmax, pnz : opt
        Optional parameters that define a separate grid for psi.

    Notes
    -------

    Rphiz coordinates form a right-handed cylindrical coordinates.
    phi is in degrees and is considered as a periodic coordinate.
    phimin is where the period begins and phimax is the last data point,
    i.e. not the end of period. E.g if you have a n = 2 symmetric system
    with nphi = 180 deg and phimin = 0 deg, then phimax should be 179 deg.

    Within ASCOT5, the magnetic field is evaluated as:

    B_R = B_R' + dPsi/dz,
    B_phi = B_phi',
    B_z = B_z' + dPsi/dR,

    where ' notates input fields.
    """

    parent = "bfield"
    group  = "B_3DS"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        path = add_group(f, parent, group, desc=desc)

        # Transpose grids
        B_R = np.transpose(B_R,(1,0,2))
        B_phi = np.transpose(B_phi,(1,0,2))
        B_z = np.transpose(B_z,(1,0,2))

        # Define psigrid to be same as Bgrid if not stated otherwise.
        if(pRmin is None or pRmax is None or pnR is None or pzmin is None or
           pzmax is None or pnz is None):
            pRmin = Rmin
            pRmax = Rmax
            pnR   = nR
            pzmin = zmin
            pzmax = zmax
            pnz   = nz

        # TODO Check that inputs are consistent.

        # Actual data.
        f.create_dataset(path + "/R_min",         (1,), data=Rmin,   dtype="f8")
        f.create_dataset(path + "/R_max",         (1,), data=Rmax,   dtype="f8")
        f.create_dataset(path + "/n_R",           (1,), data=nR,     dtype="i8")

        f.create_dataset(path + "/phi_min",       (1,), data=phimin, dtype="f8")
        f.create_dataset(path + "/phi_max",       (1,), data=phimax, dtype="f8")
        f.create_dataset(path + "/n_phi",         (1,), data=nphi,   dtype="i8")

        f.create_dataset(path + "/z_min",         (1,), data=zmin,   dtype="f8")
        f.create_dataset(path + "/z_max",         (1,), data=zmax,   dtype="f8")
        f.create_dataset(path + "/n_z",           (1,), data=nz,     dtype="i8")

        f.create_dataset(path + "/psigrid_R_min", (1,), data=pRmin,  dtype="f8")
        f.create_dataset(path + "/psigrid_R_max", (1,), data=pRmax,  dtype="f8")
        f.create_dataset(path + "/psigrid_n_R",   (1,), data=pnR,    dtype="i8")

        f.create_dataset(path + "/psigrid_z_min", (1,), data=pzmin,  dtype="f8")
        f.create_dataset(path + "/psigrid_z_max", (1,), data=pzmax,  dtype="f8")
        f.create_dataset(path + "/psigrid_n_z",   (1,), data=pnz,    dtype="i8")

        f.create_dataset(path + "/psi",                 data=psiRz,  dtype="f8")
        f.create_dataset(path + "/B_R",                 data=B_R,    dtype="f8")
        f.create_dataset(path + "/B_phi",               data=B_phi,  dtype="f8")
        f.create_dataset(path + "/B_z",                 data=B_z,    dtype="f8")

        f.create_dataset(path + "/axis_R",        (1,), data=axisR,  dtype="f8")
        f.create_dataset(path + "/axis_z",        (1,), data=axisz,  dtype="f8")

        f.create_dataset(path + "/psi0",         (1,), data=psiaxis, dtype="f8")
        f.create_dataset(path + "/psi1",         (1,), data=psisepx, dtype="f8")


def read_hdf5(fn, qid):
    """
    Read 3D magnetic field input from HDF5 file.

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

    path = "bfield" + "/B_3DS-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["R_min"] = f[path]["R_min"][:]
        out["R_max"] = f[path]["R_max"][:]
        out["n_R"]   = f[path]["n_R"][:]

        out["phi_min"] = f[path]["phi_min"][:]
        out["phi_max"] = f[path]["phi_max"][:]
        out["n_phi"]   = f[path]["n_phi"][:]

        out["z_min"] = f[path]["z_min"][:]
        out["z_max"] = f[path]["z_max"][:]
        out["n_z"]   = f[path]["n_z"][:]

        out["psi"]   = f[path]["psi"][:]
        out["B_R"]   = f[path]["B_R"][:]
        out["B_phi"] = f[path]["B_phi"][:]
        out["B_z"]   = f[path]["B_z"][:]

        out["axis_R"] = f[path]["axis_R"][:]
        out["axis_z"] = f[path]["axis_z"][:]

        out["psiaxis"] = f[path]["psi0"][:]
        out["psisepx"] = f[path]["psi1"][:]

    return out

class B_3DS(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
