"""
Non-axisymmetric tokamak magnetic field HDF5 IO

File: B_3DS.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
               b_phimin, b_phimax, b_nphi,
               axisr, axisz, psi0, psi1, psi, br, bphi, bz,
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
    B(phi) = B(phi + n*(b_phimax - b_phimin)). Do note that to avoid dublicate
    data, the last points in phi axis in B data are not at b_phimax, i.e.
    br(-1,:,:) != BR(phi=b_phimax).

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
        psi : array_like (nz, nr) <br>
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (nz,nr) <br>
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
        bphi : array_like (nz,nr) <br>
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nz,nr) <br>
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

    # Define psigrid to be same as Bgrid if not stated otherwise.
    if(pRmin is None or pRmax is None or pnR is None or pzmin is None or
       pzmax is None or pnz is None):
        pRmin = Rmin
        pRmax = Rmax
        pnR   = nR
        pzmin = zmin
        pzmax = zmax
        pnz   = nz

    print(B_R.shape)
    B_R = np.transpose(B_R,(1,0,2))
    B_phi = np.transpose(B_phi,(1,0,2))
    B_z = np.transpose(B_z,(1,0,2))
    print(B_R.shape)

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("b_rmin",   (1,), data=Rmin,    dtype="f8")
        g.create_dataset("b_rmax",   (1,), data=Rmax,    dtype="f8")
        g.create_dataset("b_nr",     (1,), data=nR,      dtype="i4")
        g.create_dataset("b_phimin", (1,), data=phimin,  dtype="f8")
        g.create_dataset("b_phimax", (1,), data=phimax,  dtype="f8")
        g.create_dataset("b_nphi",   (1,), data=nphi,    dtype="i4")
        g.create_dataset("b_zmin",   (1,), data=zmin,    dtype="f8")
        g.create_dataset("b_zmax",   (1,), data=zmax,    dtype="f8")
        g.create_dataset("b_nz",     (1,), data=nz,      dtype="i4")
        g.create_dataset("psi_rmin", (1,), data=pRmin,   dtype="f8")
        g.create_dataset("psi_rmax", (1,), data=pRmax,   dtype="f8")
        g.create_dataset("psi_nr",   (1,), data=pnR,     dtype="i4")
        g.create_dataset("psi_zmin", (1,), data=pzmin,   dtype="f8")
        g.create_dataset("psi_zmax", (1,), data=pzmax,   dtype="f8")
        g.create_dataset("psi_nz",   (1,), data=pnz,     dtype="i4")
        g.create_dataset("axisr",    (1,), data=axisR,   dtype="f8")
        g.create_dataset("axisz",    (1,), data=axisz,   dtype="f8")
        g.create_dataset("psi0",     (1,), data=psiaxis, dtype="f8")
        g.create_dataset("psi1",     (1,), data=psisepx, dtype="f8")

        g.create_dataset("psi",  (pnz, pnR),     data=psiRz, dtype="f8")
        g.create_dataset("br",   (nphi, nz, nR), data=B_R,   dtype="f8")
        g.create_dataset("bphi", (nphi, nz, nR), data=B_phi, dtype="f8")
        g.create_dataset("bz",   (nphi, nz, nR), data=B_z,   dtype="f8")


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

    return out


class B_3DS(AscotData):
    """
    Object representing B_3DS data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
