"""
Non-axisymmetric tokamak neutral HDF5 IO

File: N0_3D.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               n_species, anum, znum, n0, t0, desc=None, maxwellian=1):
    """
    Write 3D neutral input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    Rmin, Rmax, phimin, phimax, zmin, zmax : real
        Edges of the uniform Rphiz-grid.
    nR, nphi, nz : int
        Number of Rphiz-grid points.
    n_species : int
        Number of neutral species.
    anum, znum : int n_species numpy array
        Mass number and charge number for the different species
    n0 : real R x phi x z x n_species numpy array
        Neutral density in Rphiz-grid.
    t0 : real R x phi x z x n_species numpy array
        Neutral temperature in Rphiz-grid.
    maxwellian : int n_species numpy array
        Whether species distribution is Maxwellian of monoenergetic

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

    # Transpose n0 and t0 from (spec, r, phi, z) to (spec, phi, z, r)
    n0 = np.transpose(n0,(0,2,3,1))
    t0 = np.transpose(t0,(0,2,3,1))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("r_min",      (1,),      data=Rmin,       dtype="f8")
        g.create_dataset("r_max",      (1,),      data=Rmax,       dtype="f8")
        g.create_dataset("n_r",        (1,),      data=nR,         dtype="i8")
        g.create_dataset("phi_min",    (1,),      data=phimin,     dtype="f8")
        g.create_dataset("phi_max",    (1,),      data=phimax,     dtype="f8")
        g.create_dataset("n_phi",      (1,),      data=nphi,       dtype="i8")
        g.create_dataset("z_min",      (1,),      data=zmin,       dtype="f8")
        g.create_dataset("z_max",      (1,),      data=zmax,       dtype="f8")
        g.create_dataset("n_z",        (1,),      data=nz,         dtype="i8")
        g.create_dataset("n_species",  (1,),      data=n_species,  dtype="i8")
        g.create_dataset("anum",       (n_species,), data=anum,    dtype="i8")
        g.create_dataset("znum",       (n_species,), data=znum,    dtype="i8")
        g.create_dataset("maxwellian", (n_species,), data=maxwellian,dtype="i8")
        g.create_dataset("n0",                    data=n0,         dtype="f8")
        g.create_dataset("t0",                    data=t0,         dtype="f8")

def read_hdf5(fn, qid):
    """
    Read 3D neutral input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.

    Returns
    -------

    Dictionary containing neutral data.
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

        out["n_species"] = f[path]["n_species"][:]

        out["anum"] = f[path]["anum"][:]
        out["znum"] = f[path]["znum"][:]
        out["maxwellian"] = f[path]["znum"][:]

        out["n0"]   = f[path]["n0"]

    return out

class N0_3D(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
