"""
Non-axisymmetric tokamak electric field HDF5 IO

File: E_3D.py
"""
import numpy as np
import h5py

from ._iohelpers.fileapi import add_group
from . E import E

def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
               er, ephi, ez, desc=None):
    """
    Write 3D electric field input in HDF5 file for trilinear interpolation.

    The toroidal angle phi is treated as a periodic coordinate meaning that
    E(phi) = E(phi + n*(b_phimax - b_phimin)). Do note that to avoid duplicate
    data, the last points in phi axis in E data are not at phimax, i.e.
    er[:,-1,:] != ER(phi=phimax).

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        rmin : float <br>
            Minimum value in R grid [m].
        rmax : float <br>
            Maximum value in R grid [m].
        nr : int <br>
            Number of R grid points.
        zmin : float <br>
            Minimum value in z grid [m].
        zmax : float <br>
            Maximum value in z grid [m].
        nz : int <br>
            Number of z grid points.
        phimin : float <br>
            Minimum value in phi grid [deg].
        phimax : float <br>
            Maximum value in phi grid [deg].
        nphi : int <br>
            Number of phi grid points.
        er : array_like (nr,nphi,nz) <br>
            Electric field R component [V/m].
        ephi : array_like (nr,nphi,nz) <br>
            Electric field phi component [V/m].
        ez : array_like (nr,nphi,nz) <br>
            Electric field z component [V/m].    
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "efield"
    group  = "E_3D"
    gname  = ""
    
    assert er.shape   == (nr,nphi,nz)
    assert ephi.shape == (nr,nphi,nz)
    assert ez.shape   == (nr,nphi,nz)
    
    er = np.transpose(er,(2,1,0))
    ephi = np.transpose(ephi,(2,1,0))
    ez = np.transpose(ez,(2,1,0))

    # Create a group for this input.
    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("rmin",          (1,),  data=rmin,   dtype="f8")
        g.create_dataset("rmax",          (1,),  data=rmax,   dtype="f8")
        g.create_dataset("nr",            (1,),  data=nr,     dtype="i4")
        g.create_dataset("phimin",        (1,),  data=phimin, dtype="f8")
        g.create_dataset("phimax",        (1,),  data=phimax, dtype="f8")
        g.create_dataset("nphi",          (1,),  data=nphi,   dtype="i4")
        g.create_dataset("zmin",          (1,),  data=zmin,   dtype="f8")
        g.create_dataset("zmax",          (1,),  data=zmax,   dtype="f8")
        g.create_dataset("nz",            (1,),  data=nz,     dtype="i4")
        g.create_dataset("er",  (nz, nphi, nr),  data=er,     dtype="f8")
        g.create_dataset("ephi",(nz, nphi, nr),  data=ephi,   dtype="f8")
        g.create_dataset("ez",  (nz, nphi, nr),  data=ez,     dtype="f8")

    return gname

def read_hdf5(fn,qid):
    """
    Read 3D electric field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "efield" + "/E_3D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    out["er"]   = np.transpose(out["er"],   (2,1,0))
    out["ephi"] = np.transpose(out["ephi"], (2,1,0))
    out["ez"]   = np.transpose(out["ez"],   (2,1,0))
    return out

class E_3D(E):
    """
    Object representing E_3D data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
