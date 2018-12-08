"""
Non-axisymmetric tokamak neutral density HDF5 IO

File: N0_3D.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotInput

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi, n0,
               desc=None):
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

    parent = "neutral"
    group  = "N0_3D"

    # Transpose n0 from (r, phi, z) to (phi, z, r)
    n0 = np.transpose(n0,(1,2,0))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("r_min",   (1,), data=Rmin,   dtype="f8")
        g.create_dataset("r_max",   (1,), data=Rmax,   dtype="f8")
        g.create_dataset("n_r",     (1,), data=nR,     dtype="i8")
        g.create_dataset("phi_min", (1,), data=phimin, dtype="f8")
        g.create_dataset("phi_max", (1,), data=phimax, dtype="f8")
        g.create_dataset("n_phi",   (1,), data=nphi,   dtype="i8")
        g.create_dataset("z_min",   (1,), data=zmin,   dtype="f8")
        g.create_dataset("z_max",   (1,), data=zmax,   dtype="f8")
        g.create_dataset("n_z",     (1,), data=nz,     dtype="i8")
        g.create_dataset("n0",            data=n0,     dtype="f8")


def read_hdf5(fn, qid):
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

    path = "neutral/N0_3D-" + qid

    with h5py.File(fn,"r") as f:
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

        out["n0"]   = f[path]["n0"]

    return out

class N0_3D(AscotInput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
