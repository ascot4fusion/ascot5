"""
Stellarator neutral density HDF5 IO

File: N0_ST.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz,
               phimin, phimax, nphi, n_periods, n0, desc=None):
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

    parent = "neutral"
    group  = "N0_ST"

    # Transpose n0 from (r, phi, z) to (z, phi, r)
    n0 = np.transpose(n0,(2,1,0))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group)

        g.create_dataset("r_min",           (1,), data=Rmin,      dtype="f8")
        g.create_dataset("r_max",           (1,), data=Rmax,      dtype="f8")
        g.create_dataset("n_r",             (1,), data=nR,        dtype="i8")
        g.create_dataset("phi_min",         (1,), data=phimin,    dtype="f8")
        g.create_dataset("phi_max",         (1,), data=phimax,    dtype="f8")
        g.create_dataset("n_phi",           (1,), data=nphi,      dtype="i8")
        g.create_dataset("z_min",           (1,), data=zmin,      dtype="f8")
        g.create_dataset("z_max",           (1,), data=zmax,      dtype="f8")
        g.create_dataset("n_z",             (1,), data=nz,        dtype="i8")
        g.create_dataset("toroidalPeriods",       data=n_periods, dtype="i4")
        g.create_dataset("n0",                    data=n0,        dtype="f8")


def read_hdf5(fn, qid):
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

    path = "neutral/N0_ST-" + qid

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

        out["n_periods"] = f[path]["toroidalPeriods"][:]

    return out

class N0_STS(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
