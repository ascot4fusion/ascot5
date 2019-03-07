"""
Non-axisymmetric tokamak neutral HDF5 IO

File: N0_3D.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
               nspecies, anum, znum, density, temperature, maxwellian=1,
               desc=None):
    """
    Write 3D neutral input in HDF5 file.

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
        nspecies : int <br>
            Number of neutral species.
        anum : array_like (nspecies,1) <br>
            Neutral species' atomic mass number.
        znum array_like (nspecies,1) <br>
            Neutral species' charge number.
        density array_like (nspecies,nr,nphi,nz) <br>
            Neutral species-wise density [m^-3].
        temperature array_like (nr,nphi,nz,nspecies) <br>
            Neutral species-wise temperature [eV].
        maxwellian array_like (nspecies,1) <br> :
            Whether species distribution is Maxwellian (1) of monoenergetic (0)
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "neutral"
    group  = "N0_3D"

    # Transpose n0 and t0 from (spec, r, phi, z) to (spec, phi, z, r)
    density     = np.transpose(density,(0,2,3,1))
    temperature = np.transpose(temperature,(0,2,3,1))

    if maxwellian == 1:
        maxwellian = np.ones( (nspecies,1) )

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("rmin",     (1,), data=rmin,     dtype="f8")
        g.create_dataset("rmax",     (1,), data=rmax,     dtype="f8")
        g.create_dataset("nr",       (1,), data=nr,       dtype="i4")
        g.create_dataset("phimin",   (1,), data=phimin,   dtype="f8")
        g.create_dataset("phimax",   (1,), data=phimax,   dtype="f8")
        g.create_dataset("nphi",     (1,), data=nphi,     dtype="i4")
        g.create_dataset("zmin",     (1,), data=zmin,     dtype="f8")
        g.create_dataset("zmax",     (1,), data=zmax,     dtype="f8")
        g.create_dataset("nz",       (1,), data=nz,       dtype="i4")
        g.create_dataset("nspecies", (1,), data=nspecies, dtype="i4")

        g.create_dataset("anum",        (nspecies,),  data=anum,
                         dtype="i4")
        g.create_dataset("znum",        (nspecies,),  data=znum,
                         dtype="i4")
        g.create_dataset("maxwellian",  (nspecies,),  data=maxwellian,
                         dtype="i4")

        g.create_dataset("density",     (nspecies,nphi,nz,nr), data=density,
                         dtype="f8")
        g.create_dataset("temperature", (nspecies,nphi,nz,nr), data=temperature,
                         dtype="f8")

    return g.name


def read_hdf5(fn, qid):
    """
    Read 3D neutral input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "neutral/N0_3D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


class N0_3D(AscotData):
    """
    Object representing N0_3D data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
