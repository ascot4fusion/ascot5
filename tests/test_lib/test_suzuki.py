import ctypes
import unyt
import pytest
import numpy as np
from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "suzuki_sigma",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    restype=ctypes.c_uint32,
)

ekin_tabulated = np.array(
    [
        10.0,
        16.0,
        26.0,
        42.0,
        69.0,
        111.0,
        180.0,
        291.0,
        471.0,
        763.0,
        1236.0,
        2000.0,
    ]
)
sigma_zeff1 = np.array(
    [
        10.5,
        8.0,
        6.1,
        4.4,
        3.0,
        2.0,
        1.4,
        1.0,
        0.7,
        0.5,
        0.3,
        0.2,
    ]
)
sigma_zeff3 = np.array(
    [
        10.2,
        8.4,
        6.6,
        5.2,
        4.1,
        3.2,
        2.4,
        1.8,
        1.4,
        1.0,
        0.7,
        0.5,
    ]
)


@pytest.mark.parametrize(
    "Zeff, ekin, expected",
    [
        (1.000001, ekin_tabulated, sigma_zeff1),
        (3.0, ekin_tabulated, sigma_zeff3),
    ],
)
def test_suzuki_sigma(Zeff, ekin, expected):
    ne = 1e20
    nfe = (
        ne
        * (Zeff - 1.0)
        / (
            (5 / 0.05) * (2**2 - 2)
            + (1.5 / 0.05) * (6**2 - 6)
            + (0.5 / 0.05) * (8**2 - 8)
            + (26**2 - 26)
        )
    )
    n0 = ne - ((5 / 0.05) * 2 + (1.5 / 0.05) * 6 + (0.5 / 0.05) * 8 + 26) * nfe
    ni = np.array([n0, (5 / 0.05) * nfe, (1.5 / 0.05) * nfe, (0.5 / 0.05) * nfe, nfe])
    sigma = np.zeros(ekin.size) * unyt.m**2
    for i, e in enumerate(ekin * 1e3 * unyt.eV):
        H0_anum = 1
        anum = np.array([1, 1, 1, 1, 1])
        znum = np.array([1, 2, 6, 8, 26], dtype=np.int32)
        nion = anum.size

        EperAmu = e.to("J") / H0_anum
        Te = (20e3 * unyt.eV).to("J") if EperAmu >= 200e3 * unyt.eV else EperAmu / 10
        err = LIBASCOT.suzuki_sigma(
            sigma[i : i + 1].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            EperAmu,
            ne,
            Te,
            nion,
            ni.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            anum.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            znum.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        if err:
            assert False

    assert np.allclose(
        sigma, (expected * 1e-16 * unyt.cm**2).to("m**2"), atol=0.0, rtol=0.1
    )
