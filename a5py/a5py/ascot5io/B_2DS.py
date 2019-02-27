"""
Axisymmetric magnetic field HDF5 IO

File: B_2DS.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz,
               axisR, axisz, psiRz, psiaxis, psisepx,
               B_R, B_phi, B_z, desc=None):
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

    parent = "bfield"
    group  = "B_2DS"

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        # TODO Check that inputs are consistent.

        # Actual data.
        g.create_dataset("rmin",  (1,), data=Rmin,    dtype="f8")
        g.create_dataset("rmax",  (1,), data=Rmax,    dtype="f8")
        g.create_dataset("nr",    (1,), data=nR,      dtype="i8")
        g.create_dataset("zmin",  (1,), data=zmin,    dtype="f8")
        g.create_dataset("zmax",  (1,), data=zmax,    dtype="f8")
        g.create_dataset("nz",    (1,), data=nz,      dtype="i8")
        g.create_dataset("axisr", (1,), data=axisR,   dtype="f8")
        g.create_dataset("axisz", (1,), data=axisz,   dtype="f8")
        g.create_dataset("psi0",  (1,), data=psiaxis, dtype="f8")
        g.create_dataset("psi1",  (1,), data=psisepx, dtype="f8")

        g.create_dataset("psi",  (nz, nR), data=psiRz, dtype="f8")
        g.create_dataset("bR",   (nz, nR), data=B_R,   dtype="f8")
        g.create_dataset("bphi", (nz, nR), data=B_phi, dtype="f8")
        g.create_dataset("bz",   (nz, nR), data=B_z,   dtype="f8")


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

class B_2DS(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
