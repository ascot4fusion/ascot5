"""Stellarator magnetic field IO.
"""
import numpy as np
import h5py

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
               b_phimin, b_phimax, b_nphi, psi0, psi1,
               br, bphi, bz, psi,
               axis_phimin, axis_phimax, axis_nphi, axisr, axisz,
               psi_rmin=None, psi_rmax=None, psi_nr=None,
               psi_zmin=None, psi_zmax=None, psi_nz=None,
               psi_phimin=None, psi_phimax=None, psi_nphi=None,
               desc=None):
    """Write stellarator magnetic field input in HDF5 file.

    Parameters
    ----------
    fn : str
        Full path to the HDF5 file.
    b_rmin : float
        Magnetic field data R grid min edge [m].
    b_rmax : float
        Magnetic field data R grid max edge [m].
    b_nr : int
        Number of R grid points in magnetic field data.
    b_zmin : float
        Magnetic field data z grid min edge [m].
    b_zmax : float
        Magnetic field data z grid max edge [m].
    b_nz : int
        Number of z grid points in magnetic field data.
    b_phimin : float
        Magnetic field data phi grid min edge [deg].
    b_phimax : float
        Magnetic field data phi grid max edge [deg].
    b_nphi : int
        Number of phi grid points in magnetic field data.
    axis_phimin : float
        Magnetic axis phi grid min value [deg].
    axis_phimax : float
        Magnetic axis phi grid max value [deg].
    axis_nphi : float
        Number of points in magnetic axis phi grid.
    axisr : float
        Magnetic axis R coordinates on axis phi grid [m].
    axisz : float
        Magnetic axis z coordinates on axis phi grid [m].
    psi0 : float
        On-axis poloidal flux value [Vs/m].
    psi1 : float
        Separatrix poloidal flux value [Vs/m].
    psi : array_like (psi_nr,psi_nphi,psi_nz)
        Poloidal flux values on the Rz grid [Vs/m].
    br : array_like (b_nr,b_nphi,b_nz)
        Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
    bphi : array_like (b_nr,b_nphi,b_nz)
        Magnetic field phi component on Rz grid [T].
    bz : array_like (b_nr,b_nphi,b_nz)
        Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
    psi_rmin : float, optional
        Psi data R grid min edge [m].
    psi_rmax : float, optional
        Psi data R grid max edge [m].
    psi_nr : int, optional
        Number of R grid points in psi data.
    psi_zmin : float, optional
        Psi data z grid min edge [m].
    psi_zmax : float, optional
        Psi data z grid max edge [m].
    psi_nz : int, optional
        Number of z grid points in psi data.
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

    parent = "bfield"
    group  = "B_STS"
    gname  = ""

    # Define psigrid to be same as Bgrid if not stated otherwise.
    if(psi_rmin is None or psi_rmax is None or psi_nr is None
       or psi_phimin is None or psi_phimax is None or psi_nphi is None
       or psi_zmin is None or psi_zmax is None or psi_nz is None):
        psi_rmin   = b_rmin
        psi_rmax   = b_rmax
        psi_nr     = b_nr
        psi_phimin = b_phimin
        psi_phimax = b_phimax
        psi_nphi   = b_nphi
        psi_zmin   = b_zmin
        psi_zmax   = b_zmax
        psi_nz     = b_nz

    if psi.shape  != (psi_nr,psi_nphi,psi_nz):
        raise ValueError("Inconsistent shape foor psi.")
    if br.shape   != (b_nr,b_nphi,b_nz):
        raise ValueError("Inconsistent shape foor br.")
    if bphi.shape != (b_nr,b_nphi,b_nz):
        raise ValueError("Inconsistent shape foor bphi.")
    if bz.shape   != (b_nr,b_nphi,b_nz):
        raise ValueError("Inconsistent shape foor bz.")

    psi  = np.transpose(psi,  (2,1,0))
    br   = np.transpose(br,   (2,1,0))
    bphi = np.transpose(bphi, (2,1,0))
    bz   = np.transpose(bz,   (2,1,0))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("b_rmin",      (1,),     data=b_rmin,      dtype="f8")
        g.create_dataset("b_rmax",      (1,),     data=b_rmax,      dtype="f8")
        g.create_dataset("b_nr",        (1,),     data=b_nr,        dtype="i4")
        g.create_dataset("b_phimin",    (1,),     data=b_phimin,    dtype="f8")
        g.create_dataset("b_phimax",    (1,),     data=b_phimax,    dtype="f8")
        g.create_dataset("b_nphi",      (1,),     data=b_nphi,      dtype="i4")
        g.create_dataset("b_zmin",      (1,),     data=b_zmin,      dtype="f8")
        g.create_dataset("b_zmax",      (1,),     data=b_zmax,      dtype="f8")
        g.create_dataset("b_nz",        (1,),     data=b_nz,        dtype="i4")
        g.create_dataset("psi_rmin",    (1,),     data=psi_rmin,    dtype="f8")
        g.create_dataset("psi_rmax",    (1,),     data=psi_rmax,    dtype="f8")
        g.create_dataset("psi_nr",      (1,),     data=psi_nr,      dtype="i4")
        g.create_dataset("psi_phimin",  (1,),     data=psi_phimin,  dtype="f8")
        g.create_dataset("psi_phimax",  (1,),     data=psi_phimax,  dtype="f8")
        g.create_dataset("psi_nphi",    (1,),     data=psi_nphi,    dtype="i4")
        g.create_dataset("psi_zmin",    (1,),     data=psi_zmin,    dtype="f8")
        g.create_dataset("psi_zmax",    (1,),     data=psi_zmax,    dtype="f8")
        g.create_dataset("psi_nz",      (1,),     data=psi_nz,      dtype="i4")
        g.create_dataset("axis_phimin", (1,),     data=axis_phimin, dtype="f8")
        g.create_dataset("axis_phimax", (1,),     data=axis_phimax, dtype="f8")
        g.create_dataset("axis_nphi",   (1,),     data=axis_nphi,   dtype="i4")
        g.create_dataset("psi0",        (1,),     data=psi0,        dtype="f8")
        g.create_dataset("psi1",        (1,),     data=psi1,        dtype="f8")

        g.create_dataset("axisr", (axis_nphi,), data=axisr, dtype="f8")
        g.create_dataset("axisz", (axis_nphi,), data=axisz, dtype="f8")

        g.create_dataset("br",         (b_nz,b_nphi,b_nr),       data=br,
                         dtype="f8")
        g.create_dataset("bphi",       (b_nz,b_nphi,b_nr),       data=bphi,
                         dtype="f8")
        g.create_dataset("bz",         (b_nz,b_nphi,b_nr),       data=bz,
                         dtype="f8")
        g.create_dataset("psi",        (psi_nz,psi_nphi,psi_nr), data=psi,
                         dtype="f8")

    return gname


def read_hdf5(fn, qid):
    """
    Read stellarator magnetic field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "bfield/B_STS_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    out["psi"]  = np.transpose(out["psi"],  (2,1,0))
    out["br"]   = np.transpose(out["br"],   (2,1,0))
    out["bphi"] = np.transpose(out["bphi"], (2,1,0))
    out["bz"]   = np.transpose(out["bz"],   (2,1,0))
    return out


class B_STS(DataGroup):
    """
    Object representing B_STS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
