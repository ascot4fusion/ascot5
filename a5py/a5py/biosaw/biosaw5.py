import numpy as np
from scipy import interpolate
import ctypes
from numpy.ctypeslib import ndpointer
from copy import deepcopy
import a5py.ascot5io.B_3DS as b3d

libascot = ctypes.CDLL("libascot.so")
real_p = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")

libascot.biosaw_calc_B.restype = None
libascot.biosaw_calc_B.argtypes = [ctypes.c_int, real_p, real_p, real_p,
                                   ctypes.c_int, real_p, real_p, real_p,
                                   real_p, real_p, real_p]

def calc_B(Rv, phiv, zv, coil_x, coil_y, coil_z):
    coil_n = len(coil_x)
    coil_x = np.array(coil_x, dtype="f8")
    coil_y = np.array(coil_y, dtype="f8")
    coil_z = np.array(coil_z, dtype="f8")

    [R,phi,z] = np.meshgrid(Rv,phiv,zv)

    x = R * np.cos(phi)
    y = R * np.sin(phi)

    b3d = dict()

    b3d["b_rmin"] = Rv[0]
    b3d["b_rmax"] = Rv[-1]
    b3d["b_nr"] = Rv.size
    b3d["b_phimin"] = phiv[0]
    b3d["b_phimax"] = phiv[-1]
    b3d["b_nphi"] = phiv.size
    b3d["b_zmin"] = zv[0]
    b3d["b_zmax"] = zv[-1]
    b3d["b_nz"] = zv.size

    b3d["br"] = np.zeros(R.shape)
    b3d["bphi"] = np.zeros(R.shape)
    b3d["bz"] = np.zeros(R.shape)
    
    n = len(R.ravel())

    Bx = np.zeros(n, dtype="f8")
    By = np.zeros(n, dtype="f8")
    Bz = np.zeros(n, dtype="f8")

    libascot.biosaw_calc_B(n, x.ravel(), y.ravel(), z.ravel(), coil_n, coil_x, coil_y, coil_z, Bx, By, Bz)

    Bx = Bx.reshape((phiv.size,Rv.size,zv.size))
    By = By.reshape((phiv.size,Rv.size,zv.size))
    Bz = Bz.reshape((phiv.size,Rv.size,zv.size))

    b3d["br"] = Bx*np.cos(phi) + By*np.sin(phi)
    b3d["bphi"] = -Bx*np.sin(phi) + By*np.cos(phi)
    b3d["bz"] = Bz

    return b3d

def revolve_coils(b3d_coil, n_coil, B0=None, R0=None, width=None):
    b3d = deepcopy(b3d_coil)

    b3d["br"] = 0*b3d["br"]
    b3d["bphi"] = 0*b3d["bphi"]
    b3d["bz"] = 0*b3d["bz"]

    n_phi = b3d["b_nphi"]
    n_phi_coil = int(n_phi/n_coil)

    if width is None:
        width = 1

    for i in range(n_phi):
        for j in range(n_coil):
            for k in range(width):
                index = (k + i + j*n_phi_coil) % n_phi
                b3d["br"][index,:,:] += b3d_coil["br"][i,:,:]
                b3d["bphi"][index,:,:] += b3d_coil["bphi"][i,:,:]
                b3d["bz"][index,:,:] += b3d_coil["bz"][i,:,:]

    if B0 is not None and R0 is not None:
        Rv = np.linspace(b3d["b_rmin"],b3d["b_rmax"],b3d["b_nr"])
        zv = np.linspace(b3d["b_zmin"],b3d["b_zmax"],b3d["b_nz"])
        f=interpolate.interp2d(Rv, zv, np.mean(b3d["bphi"],0).T, kind='cubic')

        b3d["br"] = b3d["br"] * B0/f(R0,0)
        b3d["bphi"] = b3d["bphi"] * B0/f(R0,0)
        b3d["bz"] = b3d["bz"] * B0/f(R0,0)

    return b3d

def combine(fn, b2d, b3d):
    b3d.write_hdf5(fn,
                   b3d["b_rmin"], b3d["b_rmax"], b3d["b_nr"],
                   b3d["b_zmin"], b3d["b_zmax"], b3d["b_nz"],
                   b3d["b_phimin"], b3d["b_phimax"], b3d["b_nphi"],
                   b2d["axisr"], b2d["axisz"],
                   b2d["psi"], b2d["psi0"], b2d["psi1"],
                   np.transpose(b3d["br"],(1,0,2)),
                   np.transpose(b3d["bphi"],(1,0,2)),
                   np.transpose(b3d["bz"],(1,0,2)),
                   psi_rmin=b2d["rmin"], psi_rmax=b2d["rmax"], psi_nr=b2d["nr"],
                   psi_zmin=b2d["zmin"], psi_zmax=b2d["zmax"], psi_nz=b2d["nz"])
