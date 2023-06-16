import sys
import ctypes
import copy
import numpy as np
import numpy.ctypeslib as npctypes
from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy import Ascotpy

libascot = ctypes.CDLL("libascot.so")

class dist_5D_data(ctypes.Structure):
    _fields_ = [("n_r", ctypes.c_int),
                ("min_r", ctypes.c_double),
                ("max_r", ctypes.c_double),
                ("n_phi", ctypes.c_int),
                ("min_phi", ctypes.c_double),
                ("max_phi", ctypes.c_double),
                ("n_z", ctypes.c_int),
                ("min_z", ctypes.c_double),
                ("max_z", ctypes.c_double),
                ("n_ppara", ctypes.c_int),
                ("min_ppara", ctypes.c_double),
                ("max_ppara", ctypes.c_double),
                ("n_pperp", ctypes.c_int),
                ("min_pperp", ctypes.c_double),
                ("max_pperp", ctypes.c_double),
                ("n_time", ctypes.c_int),
                ("min_time", ctypes.c_double),
                ("max_time", ctypes.c_double),
                ("n_q", ctypes.c_int),
                ("min_q", ctypes.c_double),
                ("max_q", ctypes.c_double),
                ("histogram", ctypes.POINTER(ctypes.c_double))]

    def __init__(self, dist):
        self.n_r = dist["nr"]
        self.min_r = dist["r_edges"][0]
        self.max_r = dist["r_edges"][-1]
        self.n_phi = dist["nphi"]
        self.min_phi = dist["phi_edges"][0]
        self.max_phi = dist["phi_edges"][-1]
        self.n_z = dist["nz"]
        self.min_z = dist["z_edges"][0]
        self.max_z = dist["z_edges"][-1]
        self.n_ppara = dist["nppar"]
        self.min_ppara = dist["ppar_edges"][0]
        self.max_ppara = dist["ppar_edges"][-1]
        self.n_pperp = dist["npperp"]
        self.min_pperp = dist["pperp_edges"][0]
        self.max_pperp = dist["pperp_edges"][-1]
        self.n_time = dist["ntime"]
        self.min_time = dist["time_edges"][0]
        self.max_time = dist["time_edges"][-1]
        self.n_q = dist["ncharge"]
        self.min_q = dist["charge_edges"][0]
        self.max_q = dist["charge_edges"][-1]
        self.histogram = npctypes.as_ctypes(np.ascontiguousarray(dist["histogram"].ravel(), dtype="f8"))

class afsi_thermal_data(ctypes.Structure):
    _fields_ = [("n_r", ctypes.c_int),
                ("min_r", ctypes.c_double),
                ("max_r", ctypes.c_double),
                ("n_phi", ctypes.c_int),
                ("min_phi", ctypes.c_double),
                ("max_phi", ctypes.c_double),
                ("n_z", ctypes.c_int),
                ("min_z", ctypes.c_double),
                ("max_z", ctypes.c_double),
                ("temperature", ctypes.POINTER(ctypes.c_double)),
                ("density", ctypes.POINTER(ctypes.c_double))]

    def __init__(self, Rv, phiv, zv, temperature, density):
        self.n_r = len(Rv)-1
        self.min_r = Rv[0]
        self.max_r = Rv[-1]
        self.n_phi = len(phiv)-1
        self.min_phi = phiv[0]
        self.max_phi = phiv[-1]
        self.n_z = len(zv)-1
        self.min_z = zv[0]
        self.max_z = zv[-1]
        self.temperature = npctypes.as_ctypes(np.ascontiguousarray(temperature.ravel(), dtype="f8"))
        self.density = npctypes.as_ctypes(np.ascontiguousarray(density.ravel(), dtype="f8"))

class afsi_data(ctypes.Structure):
    _fields_ = [("type", ctypes.c_int),
                ("dist_5D", ctypes.POINTER(dist_5D_data)),
                ("dist_thermal", ctypes.POINTER(afsi_thermal_data))]

    def __init__(self, dist_5D=None, dist_thermal=None):
        if dist_5D is not None:
            self.type = 1
            self.dist_5D = ctypes.pointer(dist_5D)
        elif dist_thermal is not None:
            self.type = 2
            self.dist_thermal = ctypes.pointer(dist_thermal)
        else:
            self.type = 0

class afsi_dist_5D(ctypes.Structure):
    _fields_ = [("n_r", ctypes.c_int),
                ("min_r", ctypes.c_double),
                ("max_r", ctypes.c_double),
                ("n_phi", ctypes.c_int),
                ("min_phi", ctypes.c_double),
                ("max_phi", ctypes.c_double),
                ("n_z", ctypes.c_int),
                ("min_z", ctypes.c_double),
                ("max_z", ctypes.c_double),
                ("n_pitch", ctypes.c_int),
                ("min_pitch", ctypes.c_double),
                ("max_pitch", ctypes.c_double),
                ("n_energy", ctypes.c_int),
                ("min_energy", ctypes.c_double),
                ("max_energy", ctypes.c_double),
                ("histogram", ctypes.POINTER(ctypes.c_double))]
                #("histogram", npctypes.ndpointer(dtype=np.float64,ndim=1,flags="C_CONTIGUOUS"))]

class afsi_dist_6D(ctypes.Structure):
    _fields_ = [("n_r", ctypes.c_int),
                ("min_r", ctypes.c_double),
                ("max_r", ctypes.c_double),
                ("n_phi", ctypes.c_int),
                ("min_phi", ctypes.c_double),
                ("max_phi", ctypes.c_double),
                ("n_z", ctypes.c_int),
                ("min_z", ctypes.c_double),
                ("max_z", ctypes.c_double),
                ("n_vr", ctypes.c_int),
                ("min_vr", ctypes.c_double),
                ("max_vr", ctypes.c_double),
                ("n_vphi", ctypes.c_int),
                ("min_vphi", ctypes.c_double),
                ("max_vphi", ctypes.c_double),
                ("n_vz", ctypes.c_int),
                ("min_vz", ctypes.c_double),
                ("max_vz", ctypes.c_double),
                ("histogram", ctypes.POINTER(ctypes.c_double))]

libascot.afsi_mc.restype = None
libascot.afsi_mc.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(afsi_data), ctypes.POINTER(afsi_data), ctypes.POINTER(afsi_dist_5D)]

libascot.afsi_test_dist.restype = None
libascot.afsi_test_dist.argtypes = [ctypes.POINTER(dist_5D_data)]


def afsi_thermal(reaction, fn1, mult, nmc, it, ispecies1=1, ispecies2=1):
    """
    Calculate thermonuclear fusion.

    Args:
        reaction : int <br>
            Fusion reaction index
        fn1 : str <br>
            Path to hdf5 file.
        mult : int <br>
            Multiplier for fusion rate (1/1+dij), 0.5 if identical populations
        nmc : int <br>
            Number of MC samples.
        it : int <br>
            Time index from 5D distribution
        ispecies1 : int <br>
            Species index for first reactant.
        ispecies2 : int <br>
            Species index for second reactant.

    Returns:
        Fusion product distribution.
    """
    d1=Ascot(fn1).active.dist5d.read()
    d1["ntime"]=1
    d1["time_edges"]=d1["time_edges"][it:it+1]
    d1["histogram"]=d1["histogram"][:,:,:,:,:,[it],:]
    dist1 = dist_5D_data(d1)

    apy=Ascotpy(fn1)
    apy.init(bfield=True,plasma=True)
    temp = np.zeros((dist1.n_r, dist1.n_phi, dist1.n_z))
    dens1 = np.zeros((dist1.n_r, dist1.n_phi, dist1.n_z))
    dens2 = np.zeros((dist1.n_r, dist1.n_phi, dist1.n_z))
    for i in range(dist1.n_r):
        for j in range(dist1.n_phi):
            for k in range(dist1.n_z):
                temp[i,j,k] = apy.evaluate(d1["r"][i],0,d1["z"][k],d1["time"][it],"ti1")
                dens1[i,j,k] = apy.evaluate(d1["r"][i],0,d1["z"][k],d1["time"][it],"ni"+str(ispecies1))
                dens2[i,j,k] = apy.evaluate(d1["r"][i],0,d1["z"][k],d1["time"][it],"ni"+str(ispecies2))

    thermal1 = afsi_thermal_data(d1["r_edges"],d1["phi_edges"],d1["z_edges"],temp,dens1)
    thermal2 = afsi_thermal_data(d1["r_edges"],d1["phi_edges"],d1["z_edges"],temp,dens2)

    data_th1 = afsi_data(dist_thermal=thermal1)
    data_th2 = afsi_data(dist_thermal=thermal2)

    fusiondist=init_fusiondist(dist1)

    libascot.afsi_mc(reaction, nmc, data_th1, data_th2, fusiondist)
    fusiondist.histogram=mult*npctypes.as_array(fusiondist.histogram,shape=[dist1.n_r*dist1.n_phi*dist1.n_z]).reshape(dist1.n_phi,dist1.n_z,dist1.n_r)
    return fusiondist


def afsi_beamthermal(reaction, fn1, mult, nmc, it, ispecies=1):
    """
    Calculate beam-thermal fusion.

    Args:
        reaction : int <br>
            Fusion reaction index
        fn1 : str <br>
            Path to hdf5 file.
        mult : int <br>
            Multiplier for fusion rate (1/1+dij), 0.5 if identical populations
        nmc : int <br>
            Number of MC samples.
        it : int <br>
            Time index from 5D distribution
        ispecies : int <br>
            Species index for thermal reactant.

    Returns:
        Fusion product distribution.
    """

    d1=Ascot(fn1).active.dist5d.read()
    d1["ntime"]=1
    d1["time_edges"]=d1["time_edges"][it:it+1]
    d1["histogram"]=d1["histogram"][:,:,:,:,:,[it],:]
    dist1 = dist_5D_data(d1)

    apy=Ascotpy(fn1)
    apy.init(bfield=True,plasma=True)
    temp = np.zeros((dist1.n_r, dist1.n_phi, dist1.n_z))
    dens = np.zeros((dist1.n_r, dist1.n_phi, dist1.n_z))
    for i in range(dist1.n_r):
        for j in range(dist1.n_phi):
            for k in range(dist1.n_z):
                temp[i,j,k] = apy.evaluate(d1["r"][i],0,d1["z"][k],d1["time"][it],"ti1")
                dens[i,j,k] = apy.evaluate(d1["r"][i],0,d1["z"][k],d1["time"][it],"ni"+ispecies)
    thermal = afsi_thermal_data(d1["r_edges"],d1["phi_edges"],d1["z_edges"],temp,dens,dens)

    data_fast1 = afsi_data(dist_5D=dist1)
    data_th = afsi_data(dist_thermal=thermal)

    fusiondist=init_fusiondist(dist1)

    libascot.afsi_mc(reaction, nmc, data_th, data_fast1, fusiondist)
    fusiondist.histogram=mult*npctypes.as_array(fusiondist.histogram,shape=[dist1.n_r*dist1.n_phi*dist1.n_z]).reshape(dist1.n_phi,dist1.n_z,dist1.n_r)
    return fusiondist

def afsi_beambeam(reaction, fn1, mult, nmc, it, fn2=None):
    """
    Calculate beam-beam fusion.

    Args:
        reaction : int <br>
            Fusion reaction index
        fn1 : str <br>
            Path to hdf5 file.
        mult : int <br>
            Multiplier for fusion rate (1/1+dij), 0.5 if identical populations
        nmc : int <br>
            Number of MC samples.
        it : int <br>
            Time index from 5D distribution
        fn2 : str <br>
            Path to second hdf5 file.

    Returns:
        Fusion product distribution.
    """
    d1=Ascot(fn1).active.dist5d.read()
    d1["ntime"]=1
    d1["time_edges"]=d1["time_edges"][it:it+1]
    d1["histogram"]=d1["histogram"][:,:,:,:,:,[it],:]
    dist1 = dist_5D_data(d1)

    if fn2 is not None:
        d2=Ascot(fn2).active.dist5d.read()
        d2["ntime"]=1
        d2["time_edges"]=d2["time_edges"][it:it+1]
        d2["histogram"]=d2["histogram"][:,:,:,:,:,[it],:]
        dist2 = dist_5D_data(d2)

    data_fast1 = afsi_data(dist_5D=dist1)
    if fn2 is not None:
        data_fast2 = afsi_data(dist_5D=dist2)

    fusiondist=init_fusiondist(dist1)

    if fn2 is not None:
        libascot.afsi_mc(reaction, nmc, data_fast1, data_fast2, fusiondist)
    else:
        libascot.afsi_mc(reaction, nmc, data_fast1, data_fast1, fusiondist)

    fusiondist.histogram=mult*npctypes.as_array(fusiondist.histogram,shape=[dist1.n_r*dist1.n_phi*dist1.n_z]).reshape(dist1.n_phi,dist1.n_z,dist1.n_r)
    return fusiondist

def init_fusiondist(dist1):
    fusiondist = afsi_dist_5D()
    fusiondist.n_r = dist1.n_r
    fusiondist.min_r = dist1.min_r
    fusiondist.max_r = dist1.max_r
    fusiondist.n_phi = dist1.n_phi
    fusiondist.min_phi = dist1.min_phi
    fusiondist.max_phi = dist1.max_phi
    fusiondist.n_z = dist1.n_z
    fusiondist.min_z = dist1.min_z
    fusiondist.max_z = dist1.max_z
    fusiondist.n_pitch = 1
    fusiondist.min_pitch = 0.0
    fusiondist.max_pitch = 1.0
    fusiondist.n_energy = 1
    fusiondist.min_energy = 0.0
    fusiondist.max_energy = 20e6
    fusiondist.histogram = npctypes.as_ctypes(np.ascontiguousarray(np.zeros(dist1.n_r*dist1.n_phi*dist1.n_z), dtype="f8"))
