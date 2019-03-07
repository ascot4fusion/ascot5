"""
Axisymmetric magnetic field HDF5 IO

File: B_2DS.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz,
               axisr, axisz, psi0, psi1,
               psi, br, bphi, bz, desc=None):
    """
    Write 2DS magnetic field input in HDF5 file.

    Note that br and bz should not include the equilibrium component of the
    magnetic field as that is calculated from psi by ASCOT5 during the
    simulation.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        rmin : float <br>
            R grid min edge [m].
        rmax : float <br>
            R grid max edge [m].
        nr : int <br>
            Number of R grid points.
        zmin : float <br>
            z grid min edge [m].
        zmax : float <br>
            z grid max edge [m].
        nz : int <br>
            Number of z grid points.
        axisr : float <br>
            Magnetic axis R coordinate [m].
        axisz : float <br>
            Magnetic axis z coordinate [m].
        psi0 : float <br>
            On-axis poloidal flux value [Vs/m].
        psi1 : float <br>
            Separatrix poloidal flux value [Vs/m].
        psi : array_like (nz, nr) <br>
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (nz,nr) <br>
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
        bphi : array_like (nz,nr) <br>
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nz,nr) <br>
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "bfield"
    group  = "B_2DS"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("rmin",  (1,), data=rmin,  dtype="f8")
        g.create_dataset("rmax",  (1,), data=rmax,  dtype="f8")
        g.create_dataset("nr",    (1,), data=nr,    dtype="i4")
        g.create_dataset("zmin",  (1,), data=zmin,  dtype="f8")
        g.create_dataset("zmax",  (1,), data=zmax,  dtype="f8")
        g.create_dataset("nz",    (1,), data=nz,    dtype="i4")
        g.create_dataset("axisr", (1,), data=axisr, dtype="f8")
        g.create_dataset("axisz", (1,), data=axisz, dtype="f8")
        g.create_dataset("psi0",  (1,), data=psi0,  dtype="f8")
        g.create_dataset("psi1",  (1,), data=psi1,  dtype="f8")

        g.create_dataset("psi",  (nz, nR), data=psi,  dtype="f8")
        g.create_dataset("br",   (nz, nR), data=br,   dtype="f8")
        g.create_dataset("bphi", (nz, nR), data=bphi, dtype="f8")
        g.create_dataset("bz",   (nz, nR), data=bz,   dtype="f8")

    return g.name


def read_hdf5(fn, qid):
    """
    Read 2D magnetic field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "bfield/B_2DS_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


class B_2DS(AscotData):
    """
    Object representing B_2DS data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
