"""
Analytic tokamak field HDF5 IO.

File: B_GS.py
"""
import numpy as np
import h5py
import random
import datetime
import a5py.preprocessing.analyticequilibrium as psifun
import a5py.ascot5io.B_2DS as B_2DS
import a5py.ascot5io.B_3DS as B_3DS

from . ascot5file import add_group

from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, r0, z0, bphi0, psimult, psicoeff,
               nripple=0, a0=2, alpha0=2, delta0=0.05, desc=None):
    """
    Write analytical tokamak magnetic field input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        r0 : float <br>
            Magnetic axis R coordinate [m].
        z0 : float <br>
            Magnetic axis z coordinate [m].
        bphi0 : float <br>
            Toroidal field at axis [T].
        psimult : float
            Scaling factor for psi.
        coefficients : array_like (13,1) <br>
            Coefficients defining psi [c0, c1, ..., c11, A].
        nripple : float, optional <br>
            Number of TF coils.
        a0 : float, optional <br>
            Minor radius. [m]
        alpha0 : float, optional <br>
            Ripple penetration.
        delta0 : float, optional <br>
            Ripple (maximum) strength.
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "bfield"
    group  = "B_GS"

    c = psi_coeff # For shorter notation.

    # Search for magnetic axis psi
    x = psifun.find_axis(r0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                         c[7], c[8], c[9], c[10], c[11], c[12])
    psi0 = psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                       c[5], c[6], c[7], c[8], c[9], c[10], c[11],
                       c[12]) * psi_mult

    psi1 = 0 # Always true at the separatrix

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("r0",           (1,),   data=r0,       dtype='f8')
        g.create_dataset("z0",           (1,),   data=z0,       dtype='f8')
        g.create_dataset("bphi0",        (1,),   data=bphi0,    dtype='f8')
        g.create_dataset("psi0",         (1,),   data=psi0,     dtype='f8')
        g.create_dataset("psi1",         (1,),   data=psi1,     dtype='f8')
        g.create_dataset("psimult",      (1,),   data=psimult,  dtype='f8')
        g.create_dataset("coefficients", (13,1), data=psicoeff, dtype='f8')
        g.create_dataset("nripple",      (1,),   data=nripple,  dtype='i8')
        g.create_dataset("a0",           (1,),   data=a0,       dtype='f8')
        g.create_dataset("alpha0",       (1,),   data=alpha0,   dtype='f8')
        g.create_dataset("delta0",       (1,),   data=delta0,   dtype='f8')

    return g.name


def read_hdf5(fn, qid):
    """
    Read analytic magnetic field input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "bfield/B_GS_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    return out


def write_hdf5_B_2D(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
                    Rmin, Rmax, nR, zmin, zmax, nz, psi0=None, desc=None):
    """
    Write analytical tokamak magnetic field as a 2D field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    R0, z0 : real
        Magnetic axis Rz-coordinates.
    B_phi0 : real
        Toroidal field at axis.
    psi_mult : real
        Scaling factor for psi-
    psi_coeff : real 13 x 1 numpy array
        Coefficients defining psi.
    Rlim, Rmax, zmin, zmax : real
        Edges of the uniform Rz-grid.
    nR, nz : int
        Number of Rz-grid points.
    desc : str, optional <br>
        Input description.
    """

    rgrid = np.linspace(Rmin, Rmax, nR)
    zgrid = np.linspace(zmin, zmax, nz)

    zg, Rg = np.meshgrid(zgrid, rgrid);
    Rg = np.transpose(Rg)
    zg = np.transpose(zg)

    c = psi_coeff # For shorter notation.
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])

    Br = np.zeros((nz,nR))
    Bz = np.zeros((nz,nR))
    Bphi = (R0/Rg)*B_phi0

    # search for magnetic axis if not given
    if psi0 == None:
        x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
            c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.

    psi1 = 0

    B_2DS.write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz,
                    R0, z0, psiRz, psi0, psi1,
                    Br, Bphi, Bz, desc=desc)


def write_hdf5_B_3D(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
                    Nripple, a0, alpha0, delta0,
                    Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
                    psi0=None, desc=None):
    """
    Write analytical tokamak magnetic field as a 3D field input in HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    R0, z0 : real
        Magnetic axis Rz-coordinates.
    B_phi0 : real
        Toroidal field at axis.
    psi_mult : real
        Scaling factor for psi-
    psi_coeff : real 13 x 1 numpy array
        Coefficients defining psi.
    Nripple : real, optional
        Ripple period, default is 0
    a0 : real, optional
        Minor radius, default is 2 [m]
    alpha0 : real, optional
        Ripple penetration, default is 2
    delta0 : real, optional
        Ripple strength, default is 0.05
    Rlim, Rmax, zmin, zmax, phimin, phimax : real
        Edges of the uniform Rphiz-grid.
    nR, nz, nphi : int
        Number of Rphiz-grid points.
    """

    rgrid = np.linspace(Rmin, Rmax, nR)
    zgrid = np.linspace(zmin, zmax, nz)

    zg, Rg = np.meshgrid(zgrid, rgrid);
    Rg = np.transpose(Rg)
    zg = np.transpose(zg)

    c = psi_coeff # For shorter notation.
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])

    Br = np.zeros((nz, nR, nphi))
    Bz = np.zeros((nz, nR, nphi))
    Bphi = np.zeros((nz, nR, nphi))

    # search for magnetic axis if not given
    if psi0 == None:
        x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
            c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.

    psi1 = 0

    axisRz = np.array([R0, z0])
    axispsi = psi_mult*psifun.psi0(R0/R0,z0/R0,c[0],c[1],c[2],c[3],c[4],c[5],
                                   c[6],c[7],c[8],c[9],c[10],c[11],c[12])
    psivals = np.array([axispsi, 0.0])

    radius = np.sqrt( ( Rg - R0 ) * ( Rg - R0 ) + ( zg - z0 ) * ( zg - z0 ))
    theta = np.arctan2( zg - z0, Rg - R0 )
    delta = delta0 * np.exp(-0.5*theta*theta) * np.power( radius / a0, alpha0 )

    phigrid = np.linspace(phimin,phimax,nphi+1)
    phigrid = phigrid[0:-1]

    for i in range(0,nphi):
        Bphi[:,:,i] = ((R0/Rg)*B_phi0 * ( 1 + delta * np.cos(Nripple * phigrid[i]) ))

    Br = np.transpose(Br,(0,2,1))
    Bphi = np.transpose(Bphi,(0,2,1))
    Bz = np.transpose(Bz,(0,2,1))

    B_3DS.write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
                    R0, z0, psiRz, psi0, psi1,
                    Br, Bphi, Bz, desc=desc)


class B_GS(AscotData):
    """
    Object representing B_GS data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
