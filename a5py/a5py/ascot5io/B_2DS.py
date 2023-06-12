"""Axisymmetric magnetic field HDF5 IO
"""
import numpy as np
import h5py

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz,
               axisr, axisz, psi, psi0, psi1,
               br, bphi, bz, desc=None):
    """Write 2DS magnetic field input in HDF5 file.

    Note that br and bz should not include the equilibrium component of the
    magnetic field as that is calculated from psi by ASCOT5 during the
    simulation.

    Parameters
    ----------
    fn : str
        Full path to the HDF5 file.
    rmin : float
        R grid min edge [m].
    rmax : float
        R grid max edge [m].
    nr : int
        Number of R grid points.
    zmin : float
        z grid min edge [m].
    zmax : float
        z grid max edge [m].
    nz : int
        Number of z grid points.
    axisr : float
        Magnetic axis R coordinate [m].
    axisz : float
        Magnetic axis z coordinate [m].
    psi0 : float
        On-axis poloidal flux value [Vs/m].
    psi1 : float
        Separatrix poloidal flux value [Vs/m].
    psi : array_like (nr, nz)
        Poloidal flux values on the Rz grid [Vs/m].
    br : array_like (nr,nz)
        Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
    bphi : array_like (nr,nz)
        Magnetic field phi component on Rz grid [T].
    bz : array_like (nr,nz)
        Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
    desc : str, optional
        Input description.

    Returns
    -------
    name : str
        Name, i.e. "<type>_<qid>", of the new input that was written.

    Raises
    ------
    ValueError
        If inputs were not consistent.
    """

    if psi.shape  != (nr,nz): raise ValueError("Inconsistent shape for psi.")
    if br.shape   != (nr,nz): raise ValueError("Inconsistent shape for br.")
    if bphi.shape != (nr,nz): raise ValueError("Inconsistent shape for bphi.")
    if bz.shape   != (nr,nz): raise ValueError("Inconsistent shape for bz.")

    psi  = np.transpose(psi)
    br   = np.transpose(br)
    bphi = np.transpose(bphi)
    bz   = np.transpose(bz)

    parent = "bfield"
    group  = "B_2DS"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

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

        g.create_dataset("psi",  (nz, nr), data=psi,  dtype="f8")
        g.create_dataset("br",   (nz, nr), data=br,   dtype="f8")
        g.create_dataset("bphi", (nz, nr), data=bphi, dtype="f8")
        g.create_dataset("bz",   (nz, nr), data=bz,   dtype="f8")

    return gname


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

    out["psi"]  = np.transpose(out["psi"])
    out["br"]   = np.transpose(out["br"])
    out["bphi"] = np.transpose(out["bphi"])
    out["bz"]   = np.transpose(out["bz"])
    return out


def write_hdf5_dummy(fn, desc="Dummy"):
    import a5py.ascot5io.B_GS as B_GS
    return B_GS.write_hdf5_dummy(fn, kind="2DS", desc=desc)


class B_2DS(DataGroup):
    """
    Object representing B_2DS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
