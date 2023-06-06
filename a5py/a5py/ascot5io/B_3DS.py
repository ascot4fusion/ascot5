"""
Non-axisymmetric tokamak magnetic field HDF5 IO

File: B_3DS.py
"""
import numpy as np
import h5py

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

from . import B_2DS

def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
               b_phimin, b_phimax, b_nphi,
               axisr, axisz, psi, psi0, psi1, br, bphi, bz,
               psi_rmin=None, psi_rmax=None, psi_nr=None,
               psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
    """
    Write 3DS magnetic field input in HDF5 file.

    Note that br and bz should not include the equilibrium component of the
    magnetic field as that is calculated from psi by ASCOT5 during the
    simulation.

    It is possible to use different Rz grids for psi and magnetic field
    components by giving Rz grid for psi separately.

    The toroidal angle phi is treated as a periodic coordinate meaning that
    B(phi) = B(phi + n*(b_phimax - b_phimin)). Do note that to avoid duplicate
    data, the last points in phi axis in B data are not at b_phimax, i.e.
    br[:,-1,:] != BR(phi=b_phimax).

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        b_rmin : float <br>
            Magnetic field data R grid min edge [m].
        b_rmax : float <br>
            Magnetic field data R grid max edge [m].
        b_nr : int <br>
            Number of R grid points in magnetic field data.
        b_zmin : float <br>
            Magnetic field data z grid min edge [m].
        b_zmax : float <br>
            Magnetic field data z grid max edge [m].
        b_nz : int <br>
            Number of z grid points in magnetic field data.
        b_phimin : float <br>
            Magnetic field data phi grid min edge [deg].
        b_phimax : float <br>
            Magnetic field data phi grid max edge [deg].
        b_nphi : int <br>
            Number of phi grid points in magnetic field data.
        axisr : float <br>
            Magnetic axis R coordinate [m].
        axisz : float <br>
            Magnetic axis z coordinate [m].
        psi0 : float <br>
            On-axis poloidal flux value [Vs/m].
        psi1 : float <br>
            Separatrix poloidal flux value [Vs/m].
        psi : array_like (nr, nz) <br>
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (nr,nphi,nz) <br>
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
        bphi : array_like (nr,nphi,nz) <br>
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nr,nphi,nz) <br>
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
        psi_rmin : float, optional <br>
            Psi data R grid min edge [m].
        psi_rmax : float, optional <br>
            Psi data R grid max edge [m].
        psi_nr : int, optional <br>
            Number of R grid points in psi data.
        psi_zmin : float, optional <br>
            Psi data z grid min edge [m].
        psi_zmax : float, optional <br>
            Psi data z grid max edge [m].
        psi_nz : int, optional <br>
            Number of z grid points in psi data.
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "bfield"
    group  = "B_3DS"
    gname  = ""

    # Define psigrid to be same as Bgrid if not stated otherwise.
    if(psi_rmin is None or psi_rmax is None or psi_nr is None or
       psi_zmin is None or psi_zmax is None or psi_nz is None):
        psi_rmin = b_rmin
        psi_rmax = b_rmax
        psi_nr   = b_nr
        psi_zmin = b_zmin
        psi_zmax = b_zmax
        psi_nz   = b_nz

    assert psi.shape  == (psi_nr,psi_nz)
    assert br.shape   == (b_nr,b_nphi,b_nz)
    assert bphi.shape == (b_nr,b_nphi,b_nz)
    assert bz.shape   == (b_nr,b_nphi,b_nz)

    psi  = np.transpose(psi)
    br   = np.transpose(br,   (2,1,0))
    bphi = np.transpose(bphi, (2,1,0))
    bz   = np.transpose(bz,   (2,1,0))

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("b_rmin",   (1,), data=b_rmin,   dtype="f8")
        g.create_dataset("b_rmax",   (1,), data=b_rmax,   dtype="f8")
        g.create_dataset("b_nr",     (1,), data=b_nr,     dtype="i4")
        g.create_dataset("b_phimin", (1,), data=b_phimin, dtype="f8")
        g.create_dataset("b_phimax", (1,), data=b_phimax, dtype="f8")
        g.create_dataset("b_nphi",   (1,), data=b_nphi,   dtype="i4")
        g.create_dataset("b_zmin",   (1,), data=b_zmin,   dtype="f8")
        g.create_dataset("b_zmax",   (1,), data=b_zmax,   dtype="f8")
        g.create_dataset("b_nz",     (1,), data=b_nz,     dtype="i4")
        g.create_dataset("psi_rmin", (1,), data=psi_rmin, dtype="f8")
        g.create_dataset("psi_rmax", (1,), data=psi_rmax, dtype="f8")
        g.create_dataset("psi_nr",   (1,), data=psi_nr,   dtype="i4")
        g.create_dataset("psi_zmin", (1,), data=psi_zmin, dtype="f8")
        g.create_dataset("psi_zmax", (1,), data=psi_zmax, dtype="f8")
        g.create_dataset("psi_nz",   (1,), data=psi_nz,   dtype="i4")
        g.create_dataset("axisr",    (1,), data=axisr,    dtype="f8")
        g.create_dataset("axisz",    (1,), data=axisz,    dtype="f8")
        g.create_dataset("psi0",     (1,), data=psi0,     dtype="f8")
        g.create_dataset("psi1",     (1,), data=psi1,     dtype="f8")

        g.create_dataset("psi",  (psi_nz, psi_nr),     data=psi,  dtype="f8")
        g.create_dataset("br",   (b_nz, b_nphi, b_nr), data=br,   dtype="f8")
        g.create_dataset("bphi", (b_nz, b_nphi, b_nr), data=bphi, dtype="f8")
        g.create_dataset("bz",   (b_nz, b_nphi, b_nr), data=bz,   dtype="f8")

    return gname


def read_hdf5(fn, qid):
    """
    Read 3D magnetic field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "bfield/B_3DS_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    out["psi"]  = np.transpose(out["psi"])
    out["br"]   = np.transpose(out["br"],   (2,1,0))
    out["bphi"] = np.transpose(out["bphi"], (2,1,0))
    out["bz"]   = np.transpose(out["bz"],   (2,1,0))
    return out


def write_hdf5_dummy(fn, desc="Dummy"):
    import a5py.ascot5io.B_GS as B_GS
    return B_GS.write_hdf5_dummy(fn, kind="3DS", desc=desc)


def write_hdf5_B_2DS(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                     b_phimin, b_phimax, b_nphi,
                     axisr, axisz, psi, psi0, psi1, br, bphi, bz,
                     desc=None):
    """
    Convert the 3D data to 2D.
    """
    bphi = bphi=np.mean(bphi,axis=1)
    B_2DS.write_hdf5(fn, rmin=b_rmin, rmax=b_rmax,
                     zmin=b_zmin, zmax=b_zmax, nr=b_nr,
                     nz=b_nz, axisr=axisr, axisz=axisz,
                     psi0=psi0, psi1=psi1, psi=psi,
                     br=bphi*0, bphi=bphi, bz=bphi*0)



class B_3DS(DataGroup):
    """
    Object representing B_3DS data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
