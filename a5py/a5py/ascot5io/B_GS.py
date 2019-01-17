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

def write_hdf5(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
               Nripple=0, a0=2, alpha0=2, delta0=0.05, psi0=None, desc=None):
    """
    Write analytical tokamak magnetic field input in HDF5 file.

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

    Notes
    -------

    For details of 2D field, see a5py.preprocessing.analyticfield.

    3D field consists of toroidal field ripple with Nripple toroidal
    field coils, If Nripple = 0, the field is axisymmetric.

    The ripple has the form

        Bphi = Bphi_2D * ( 1 + delta * cos(Nripple * phi) ),

    where delta = delta0 * y1(theta) * y2(r) with r being the distance
    from magnetic axis, i.e. minor radius coordinate, and theta is
    poloidal angle. Functions y1 and y2 are

        y1 = exp(-theta^2 / 2),
        y2 = ( r/a0 )^alpha0.

    Since ripple does not have B_R and B_z components, the 3D field is
    not divergence free and one cannot see ripple in Poincare plots.
    """

    parent = "bfield"
    group  = "B_GS"

    c = psi_coeff # For shorter notation.

    # Search for magnetic axis if not given.
    if psi0 == None:
        x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                             c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                           c[5], c[6], c[7], c[8], c[9], c[10], c[11],
                           c[12]) * psi_mult

    psi1 = 0 # Always true at the separatrix

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("R0",        (1,), data=R0,        dtype='f8')
        g.create_dataset("z0",        (1,), data=z0,        dtype='f8')
        g.create_dataset("B_phi0",    (1,), data=B_phi0,    dtype='f8')
        g.create_dataset("psi0",      (1,), data=psi0,      dtype='f8')
        g.create_dataset("psi1",      (1,), data=psi1,      dtype='f8')
        g.create_dataset("psi_mult",  (1,), data=psi_mult,  dtype='f8')
        g.create_dataset("psi_coeff",       data=psi_coeff, dtype='f8')
        g.create_dataset("Nripple",   (1,), data=Nripple,   dtype='i8')
        g.create_dataset("a0",        (1,), data=a0,        dtype='f8')
        g.create_dataset("alpha0",    (1,), data=alpha0,    dtype='f8')
        g.create_dataset("delta0",    (1,), data=delta0,    dtype='f8')


def read_hdf5(fn, qid):
    """
    Read analytic magnetic field input from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the bfield to be read.

    Returns
    -------

    Dictionary containing magnetic field data.
    """

    path = "bfield" + "/B_GS-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        out["R0"]        = f[path]["R0"][:]
        out["z0"]        = f[path]["z0"][:]
        out["B_phi0"]    = f[path]["B_phi0"][:]
        out["psi0"]      = f[path]["psi0"][:]
        out["psi1"]      = f[path]["psi1"][:]
        out["psi_mult"]  = f[path]["psi_mult"][:]
        out["psi_coeff"] = f[path]["psi_coeff"][:]

        out["Nripple"] = f[path]["Nripple"][:]
        out["a0"]      = f[path]["a0"][:]
        out["alpha0"]  = f[path]["alpha0"][:]
        out["delta0"]  = f[path]["delta0"][:]

    return out


def write_hdf5_B_2D(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
                    Rmin, Rmax, nR, zmin, zmax, nz, psi0=None, desc=None):
    """
    Write analytical tokamak magnetic fieldd as a 2D field input in HDF5 file.

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
    """

    rgrid = np.linspace(Rmin, Rmax, nR)
    zgrid = np.linspace(zmin, zmax, nz)

    zg, Rg = np.meshgrid(zgrid, rgrid);
    Rg = np.transpose(Rg)
    zg = np.transpose(zg)

    c = psi_coeff # For shorter notation.
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])

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
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])

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
    axispsi = psi_mult*psifun.psi0(R0/R0,z0/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12])
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

    def read(self):
        return read_hdf5(self._file, self.get_qid())
