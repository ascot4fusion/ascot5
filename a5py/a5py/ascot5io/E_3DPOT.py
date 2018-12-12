"""
Non-axisymmetric electrostatic potential HDF5 IO

File: E_3DPOT.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               tmin, tmax, nt, cyclic_time, E_pot, desc=None):
    """
    Write 3D electrostatic potential input in HDF5 file for interpolation.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rmin, Rmax, phimin, phimax, zmin, zmax, tmin, tmax : real
        Edges of the uniform Rphiz-grid and time grid.
    nR, nphi, nz, nt, cyclic_time : int
        Number of Rphiz-grid points and time, and a flag for using cyclic time data. 
    E_pot : real R x phi x z numpy array
        Electrostatic potential in Rphiz-grid.

    Notes
    -------

    Rphiz coordinates form a right-handed cylindrical coordinates.
    phi is in degrees and is considered as a periodic coordinate.
    phimin is where the period begins and phimax is the last data point,
    i.e. not the end of period. E.g if you have a n = 2 symmetric system
    with nphi = 180 deg and phimin = 0 deg, then phimax should be 179 deg.

    """

    parent = "efield"
    group  = "E_3DPOT"

    # TODO transpose grids if necessary
    E_pot = np.transpose(E_pot,(1,0,2,1))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("R_min",       (1,), data=Rmin,        dtype="f8")
        g.create_dataset("R_max",       (1,), data=Rmax,        dtype="f8")
        g.create_dataset("n_R",         (1,), data=nR,          dtype="i8")
        g.create_dataset("phi_min",     (1,), data=phimin,      dtype="f8")
        g.create_dataset("phi_max",     (1,), data=phimax,      dtype="f8")
        g.create_dataset("n_phi",       (1,), data=nphi,        dtype="i8")
        g.create_dataset("z_min",       (1,), data=zmin,        dtype="f8")
        g.create_dataset("z_max",       (1,), data=zmax,        dtype="f8")
        g.create_dataset("n_z",         (1,), data=nz,          dtype="i8")
        g.create_dataset("t_min",       (1,), data=tmin,        dtype="f8")
        g.create_dataset("t_max",       (1,), data=tmax,        dtype="f8")
        g.create_dataset("n_t",         (1,), data=nt,          dtype="i8")
        g.create_dataset("cyclic_time", (1,), data=cyclic_time, dtype="i8")
        g.create_dataset("E_pot",             data=E_POT,       dtype="f8")


def read_hdf5(fn,qid):
    """
    Read 3D electrostatic potential input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the efield potential to be read.


    Returns
    -------

    Dictionary containing electrostatic potential data.
    """

    path = "efield" + "/E_3DPOT-" + qid

    f = h5py.File(fn,"r")

    out = {}

    # Metadata.
    out["qid"] = qid
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

    out["t_min"] = f[path]["t_min"][:]
    out["t_max"] = f[path]["t_max"][:]
    out["n_t"]   = f[path]["n_t"][:]
    out["cyclic_time"]   = f[path]["cyclic_time"][:]

    out["E_pot"]   = f[path]["E_pot"][:]

    f.close()

    return out

class E_3DPOT(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
