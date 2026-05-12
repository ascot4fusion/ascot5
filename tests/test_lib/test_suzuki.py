import ctypes
import unyt
import numpy as np
from a5py.libascot import LIBASCOT, init_fun


init_fun(
    "suzuki_sigmav",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    restype=ctypes.c_uint32
)

def test_suzuki():
    n = 100
    import matplotlib.pyplot as plt
    sigma = np.zeros(n)
    ekin = np.logspace(1, 3, n)*1e3
    Zeff = 1

    for i in range(n):
        sigmav = np.zeros(1) * unyt.m**2
        h_anum = 1
        anum = np.array([1, 1, 1, 1, 1])
        znum = np.array([1, 2, 6, 8, 26], dtype=np.int32)
        nion = anum.size

        Ekin = ekin[i] * unyt.eV
        mass = 1.007 * unyt.amu
        vnorm = 1#np.sqrt(2 * Ekin / mass).to("m/s")
        EperAmu = (Ekin).to("J") / h_anum

        ne = 1e20
        Te = (20e3 * unyt.eV).to("J") if EperAmu >= 200e3 * unyt.eV else EperAmu / 10
        nfe = ne * (Zeff - 0.99999) / (  (5/0.05) * (2**2 - 2) + (1.5/0.05) * (6**2 - 6) + (0.5/0.05) * (8**2 -8) + (26**2 - 26)  )
        n0 = ne - ((5/0.05) * 2 + (1.5/0.05) * 6 + (0.5/0.05) * 8 + 26) * nfe
        ni = np.array([n0, (5/0.05)*nfe, (1.5/0.05)*nfe, (0.5/0.05)*nfe, nfe])
        #print(ni)
        #print(np.sum(ni*znum**2)/ne)

        err = LIBASCOT.suzuki_sigmav(
            sigmav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            EperAmu,
            vnorm,
            ne,
            Te,
            nion,
            ni.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            anum.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            znum.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if err:
            print("fail")
        sigma[i] = sigmav.to("cm**2")

    plt.loglog(ekin/1e3, sigma/1e-16)

    Zeff = 3

    for i in range(n):
        sigmav = np.zeros(1) * unyt.m**2
        h_anum = 1
        anum = np.array([1, 1, 1, 1, 1])
        znum = np.array([1, 2, 6, 8, 26], dtype=np.int32)
        nion = anum.size

        Ekin = ekin[i] * unyt.eV
        mass = 1.007 * unyt.amu
        vnorm = 1#np.sqrt(2 * Ekin / mass).to("m/s")
        EperAmu = (Ekin).to("J") / h_anum

        ne = 1e20
        Te = (20e3 * unyt.eV).to("J") if EperAmu >= 200e3 * unyt.eV else EperAmu / 10
        nfe = ne * (Zeff - 0.99999) / (  (5/0.05) * (2**2 - 2) + (1.5/0.05) * (6**2 - 6) + (0.5/0.05) * (8**2 -8) + (26**2 - 26)  )
        n0 = ne - ((5/0.05) * 2 + (1.5/0.05) * 6 + (0.5/0.05) * 8 + 26) * nfe
        ni = np.array([n0, (5/0.05)*nfe, (1.5/0.05)*nfe, (0.5/0.05)*nfe, nfe])
        #print(ni)
        #print(np.sum(ni*znum**2)/ne)

        err = LIBASCOT.suzuki_sigmav(
            sigmav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            EperAmu,
            vnorm,
            ne,
            Te,
            nion,
            ni.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            anum.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            znum.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if err:
            print("fail")
        sigma[i] = sigmav.to("cm**2")

    plt.loglog(ekin/1e3, sigma/1e-16)
    plt.xlabel("ekin [keV]")
    plt.ylabel("sigma [10^-16 cm^2]")
    plt.show()
    assert False
