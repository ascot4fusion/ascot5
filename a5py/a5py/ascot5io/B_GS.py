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
import a5py.ascot5io.B_STS as B_STS

from . ascot5file import add_group

from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, r0, z0, bphi0, psimult, coefficients, psi0=None, psi1=0,
               raxis=None, zaxis=None, nripple=0, a0=2, alpha0=2, delta0=0.05,
               desc=None):
    """
    Write analytical tokamak magnetic field input in HDF5 file.

    If psi0 is None, magnetic axis location is found numerically and psi0 value
    is evaluated at that point. The analytical field is defined by assuming
    midplane is at z=0, but it is moved to z=z0 here.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        r0 : float <br>
            Major radius R coordinate [m].
        z0 : float <br>
            Distance by which midplane is moved from z=0 plane [m].
        bphi0 : float <br>
            Toroidal field at axis [T].
        psimult : float
            Scaling factor for psi.
        coefficients : array_like (13,1) <br>
            Coefficients defining psi [c0, c1, ..., c11, A].
        psi0 : float, optional <br>
            Poloidal flux at magnetic axis.
        psi1 : float, optional <br>
            Poloidal flux at the separatrix.
        raxis : float, optional <br>
            Magnetic axis R coordinate [m].
        zaxis : float, optional <br>
            Magnetic axis z coordinate [m].
        nripple : float, optional, default is 0. <br>
            Number of TF coils.
        a0 : float, optional, default is 2 [m]. <br>
            Minor radius. [m]
        alpha0 : float, optional, default is 2. <br>
            Ripple penetration.
        delta0 : float, optional, default is 0.05. <br>
            Ripple (maximum) strength.
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    assert coefficients.size == 13

    parent = "bfield"
    group  = "B_GS"
    gname  = ""

    c = coefficients # For shorter notation.

    # Search for magnetic axis psi
    if psi0 is None:
        x = psifun.find_axis(r0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                             c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                           c[5], c[6], c[7], c[8], c[9], c[10], c[11],
                           c[12]) * psimult * 1.001
        raxis = x[0]*r0
        zaxis = x[1]*r0

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("r0",           (1,),   data=r0,           dtype='f8')
        g.create_dataset("z0",           (1,),   data=z0,           dtype='f8')
        g.create_dataset("raxis",        (1,),   data=raxis,        dtype='f8')
        g.create_dataset("zaxis",        (1,),   data=zaxis,        dtype='f8')
        g.create_dataset("bphi0",        (1,),   data=bphi0,        dtype='f8')
        g.create_dataset("psi0",         (1,),   data=psi0,         dtype='f8')
        g.create_dataset("psi1",         (1,),   data=psi1,         dtype='f8')
        g.create_dataset("psimult",      (1,),   data=psimult,      dtype='f8')
        g.create_dataset("coefficients", (13,1), data=coefficients, dtype='f8')
        g.create_dataset("nripple",      (1,),   data=nripple,      dtype='i8')
        g.create_dataset("a0",           (1,),   data=a0,           dtype='f8')
        g.create_dataset("alpha0",       (1,),   data=alpha0,       dtype='f8')
        g.create_dataset("delta0",       (1,),   data=delta0,       dtype='f8')

    return gname


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


def write_hdf5_dummy(fn, kind="GS", desc="Dummy"):
    r0      = 6.2
    z0      = 0
    bphi0   = 5.3
    psimult = 200

    # ITER-like but circular equilibrium
    coefficients = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                              6.200e-03, -1.205e-03, -3.701e-05,  0,
                              0,          0,          0,          0,    -0.155])

    if kind == "GS":
        return write_hdf5(fn, r0, z0, bphi0, psimult, coefficients,
                          nripple=1, a0=2, alpha0=2, delta0=0.05, desc=desc)

    if kind == "2DS":
        R = (1,  6, 50)
        z = (-4, 4, 100)
        return write_hdf5_B_2DS(fn, r0, z0, bphi0, psimult, coefficients,
                                R[0], R[1], R[2], z[0], z[1], z[2], psi0=None,
                                desc=desc)

    if kind == "3DS":
        R   = (1,  6, 50)
        z   = (-4, 4, 100)
        phi = (0, 360, 100)
        return write_hdf5_B_3DS(fn, r0, z0, bphi0, psimult, coefficients,
                                1, 2, 2, 0.05,
                                R[0], R[1], R[2], z[0], z[1], z[2],
                                phi[0], phi[1], phi[2], psi0=None, desc=desc)


def write_hdf5_B_2DS(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
                     Rmin, Rmax, nR, zmin, zmax, nz, psi0=None, desc=None):
    """
    Write analytical tokamak magnetic field as a 2D field input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        R0 : float <br>
            Major axis R coordinate [m].
        z0 : float <br>
            Major axis z coordinate [m].
        B_phi0 : float <br>
            Toroidal field at axis.
        psi_mult : real
            Scaling factor for psi.
        psi_coeff : real 13 x 1 numpy array
            Coefficients defining psi.
        Rmin, Rmax, zmin, zmax : real
            Edges of the uniform Rz-grid.
        nR, nz : int
            Number of Rz-grid points.
        psi0 : float, optional <br>
            Poloidal flux at magnetic axis.
        desc : str, optional <br>
            Input description.
    """

    rgrid = np.linspace(Rmin, Rmax, nR)
    zgrid = np.linspace(zmin, zmax, nz)

    zg, Rg = np.meshgrid(zgrid, rgrid);

    c = psi_coeff # For shorter notation.
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])

    Br = np.zeros((nR,nz))
    Bz = np.zeros((nR,nz))
    Bphi = (R0/Rg)*B_phi0

    # search for magnetic axis if not given
    if psi0 == None:
        x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
            c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
        R0 = x[0]*R0
        z0 = x[1]*R0

    psi1 = 0

    return B_2DS.write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz,
                            R0, z0, psiRz, psi0, psi1,
                            Br, Bphi, Bz, desc=desc)


def write_hdf5_B_3DS(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
                     Nripple, a0, alpha0, delta0,
                     Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
                     psi0=None, raxis=None, zaxis=None, desc=None):
    """
    Write analytical tokamak magnetic field as a 3D field input in HDF5 file.

    Args:
        fn : str
            Full path to the HDF5 file.
        R0, z0 : real
            Major axis Rz-coordinates.
        B_phi0 : real
            Toroidal field at axis.
        psi_mult : real
            Scaling factor for psi.
        psi_coeff : real 13 x 1 numpy array
            Coefficients defining psi.
        Nripple : real
            Number of TF coils.
        a0 : real
            Minor radius [m].
        alpha0 : real
            Ripple penetration.
        delta0 : real
            Ripple strength.
        Rmin, Rmax, zmin, zmax, phimin, phimax : real
            Edges of the uniform Rphiz-grid.
        nR, nz, nphi : int
            Number of Rphiz-grid points.
        psi0 : float, optional <br>
            Poloidal flux at magnetic axis.
        raxis : float, optional <br>
            Magnetic axis R coordinate [m].
        zaxis : float, optional <br>
            Magnetic axis z coordinate [m].
        desc : str, optional <br>
            Input description.
    """

    rgrid = np.linspace(Rmin, Rmax, nR)
    zgrid = np.linspace(zmin, zmax, nz)

    zg, Rg = np.meshgrid(zgrid, rgrid);

    c = psi_coeff # For shorter notation.
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])

    Br = np.zeros((nR, nphi, nz))
    Bz = np.zeros((nR, nphi, nz))
    Bphi = np.zeros((nR, nphi, nz))

    # Ripple
    radius = np.sqrt( ( Rg - R0 ) * ( Rg - R0 ) + ( zg - z0 ) * ( zg - z0 ))
    theta = np.arctan2( zg - z0, Rg - R0 )
    delta = delta0 * np.exp(-0.5*theta*theta) * np.power( radius / a0, alpha0 )

    phigrid = np.linspace(phimin,phimax,nphi+1)
    phigrid = phigrid[0:-1]

    for i in range(0,nphi):
        Bphi[:,i,:] = ((R0/Rg)*B_phi0 * ( 1 + delta * np.cos(Nripple * phigrid[i]*2*np.pi/360) ))

    # search for magnetic axis if not given
    if psi0 == None:
        x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
            c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
        raxis = x[0]*R0
        zaxis = x[1]*R0

    psi1 = 0

    return B_3DS.write_hdf5(fn, Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax,
                            nphi, raxis, zaxis, psiRz, psi0, psi1,
                            Br, Bphi, Bz, desc=desc)

def write_hdf5_B_STS(fn, R0, z0, B_phi0, psi_mult, psi_coeff,
                     Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
                     psi0=None, raxis=None, zaxis=None, desc=None):
    """
    Write analytical tokamak magnetic field as a stellarator field input in HDF5 file.

    Args:
        fn : str
            Full path to the HDF5 file.
        R0, z0 : real
            Major axis Rz-coordinates.
        B_phi0 : real
            Toroidal field at axis.
        psi_mult : real
            Scaling factor for psi.
        psi_coeff : real 13 x 1 numpy array
            Coefficients defining psi.
        Rmin, Rmax, zmin, zmax, phimin, phimax : real
            Edges of the uniform Rphiz-grid.
        nR, nz, nphi : int
            Number of Rphiz-grid points.
        psi0 : float, optional <br>
            Poloidal flux at magnetic axis.
        raxis : float, optional <br>
            Magnetic axis R coordinate [m].
        zaxis : float, optional <br>
            Magnetic axis z coordinate [m].
        desc : str, optional <br>
            Input description.
    """

    rgrid = np.linspace(Rmin, Rmax, nR)
    zgrid = np.linspace(zmin, zmax, nz)

    zg, Rg = np.meshgrid(zgrid, rgrid);

    c = psi_coeff # For shorter notation.
    psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])
    psiRz_R = psi_mult*psifun.psiX(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])
    psiRz_z = psi_mult*psifun.psiY(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                 c[7],c[8],c[9],c[10],c[11],c[12])

    Br = np.zeros((nR, nphi, nz))
    Bz = np.zeros((nR, nphi, nz))
    Bphi = np.zeros((nR, nphi, nz))
    psi  = np.zeros((nR, nphi, nz))
    axisr= np.zeros((nphi,))
    axisz= np.zeros((nphi,))

    
    phigrid = np.linspace(phimin,phimax,nphi+1)
    phigrid = phigrid[0:-1]

    # search for magnetic axis if not given
    if psi0 == None:
        x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
            c[7], c[8], c[9], c[10], c[11], c[12])
        psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
        raxis = x[0]*R0
        zaxis = x[1]*R0

    psi1=0 #??
        
    #Normalize psi (It tries to impersonate the VMEC psi)
    #print("Initial psi0: {}".format(psi0))
    #print("Initial psi1: {}".format(psi1))
    psiRz -= psi0 # move axis to 0
    psi1  -= psi0 # ..
    psiRz /= psi1 # scale LCFS to 1
    psi0 = 0.0
    psi1 = 1.0
    
    
    # The magnetic field one gets from psi as:
    #  B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
    #  B_z &=  \frac{1}{R}\frac{\partial\psi}{\partial R}
    #  f
    #  B_\phi = \frac{B_0R_0}{R}
    #  f

    for i in range(0,nphi):
        Bphi[:,i,:] = ((R0/Rg)*B_phi0 )
        psi[ :,i,:] = psiRz
        Br[  :,i,:] = ( -1.0 / Rg / R0 ) * psiRz_z 
        Bz[  :,i,:] = (  1.0 / Rg / R0 ) * psiRz_R 
        axisr[ i  ] = raxis
        axisz[ i  ] = zaxis


    return B_STS.write_hdf5(fn, b_rmin=Rmin, b_rmax=Rmax, b_nr=nR, b_zmin=zmin, b_zmax=zmax, b_nz=nz, b_phimin=phimin, b_phimax=phimax,
                            b_nphi=nphi, axisr=axisr, axisz=axisz, psi=psi, psi0=psi0, psi1=psi1,
                            br=Br, bphi=Bphi, bz=Bz, desc=desc,
                            axis_phimin=phimin, axis_phimax=phimax, axis_nphi=nphi)


class B_GS(AscotData):
    """
    Object representing B_GS data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
