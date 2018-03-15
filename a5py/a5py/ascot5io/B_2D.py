"""
Axisymmetric magnetic field HDF5 IO
"""
import numpy as np
import h5py

from . ascot5group import creategroup

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz,
               axisR, axisz, psiRz, psiaxis, psisepx,
               B_R, B_phi, B_z):
    """
    Write 2DS magnetic field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rlim, Rmax, zmin, zmax : real
        Edges of the uniform Rz-grid.
    nR, nz : int
        Number of Rz-grid points.
    axisR, axisz : real
        Magnetic axis Rz-location.
    psiRz : real R x z numpy array
        Psi values in the Rz-grid.
    psiaxis, psisepx : real    
        Psi values at magnetic axis and separatrix
    B_R, B_phi, B_z : real R x z numpy array
        Magnetic field components in Rz-grid.

    Notes
    -------

    Within ASCOT5, the magnetic field is evaluated as:
 
    B_R = B_R' + dPsi/dz,
    B_phi = B_phi',
    B_z = B_z' + dPsi/dR,

    where ' notates input fields.
    """

    mastergroup = "bfield"
    subgroup    = "B_2DS"
    
    # Create a group for this input.
    f = h5py.File(fn, "a")
    path = creategroup(f, mastergroup, subgroup)
    
    # TODO Check that inputs are consistent.

    # Actual data.
    f.create_dataset(path + "/R_min", (1,), data=Rmin, dtype="f8")
    f.create_dataset(path + "/R_max", (1,), data=Rmax, dtype="f8")
    f.create_dataset(path + "/n_R", (1,),   data=nR, dtype="i8")

    f.create_dataset(path + "/z_min", (1,), data=zmin, dtype="f8")
    f.create_dataset(path + "/z_max", (1,), data=zmax, dtype="f8")
    f.create_dataset(path + "/n_z", (1,),   data=nz, dtype="i8")

    f.create_dataset(path + "/psi",   data=psiRz, dtype="f8")
    f.create_dataset(path + "/B_R",   data=B_R, dtype="f8")
    f.create_dataset(path + "/B_phi", data=B_phi, dtype="f8")
    f.create_dataset(path + "/B_z",   data=B_z, dtype="f8")

    f.create_dataset(path + "/axis_R", (1,), data=axisR, dtype="f8")
    f.create_dataset(path + "/axis_z", (1,), data=axisz, dtype="f8")

    f.create_dataset(path + "/psi0", (1,), data=psiaxis, dtype="f8")
    f.create_dataset(path + "/psi1", (1,), data=psisepx, dtype="f8")
    f.close()


def read_hdf5(fn, qid):
    """
    Read 2D magnetic field input from HDF5 file.

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

    path = "bfield" + "/B_2DS-" + qid

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"]  = qid
    out["date"] = f[path].attrs["date"]
    out["description"] = f[path].attrs["description"]

    # Actual data.
    out["R_min"] = f[path]["R_min"][:]
    out["R_max"] = f[path]["R_max"][:]
    out["n_R"]   = f[path]["n_R"][:]

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

    f.close()

    return out
