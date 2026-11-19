import ctypes as ct
import pytest

import unyt
import numpy as np
from a5py.libascot import LIBASCOT, init_fun
from a5py.data.bfield import Bfield
from a5py.data.efield import Efield

init_fun(
    "Efield_eval_e",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Efield),
    ct.POINTER(Bfield),
    restype=ct.c_int,
)


class Helpers:

    @staticmethod
    def efield_eval_e(efield, bfield, r, phi, z, t):
        obj = Efield()
        obj.use(efield)

        bobj = Bfield()
        bobj.use(bfield)

        n = r.size
        e = np.zeros((n, 3), dtype=np.float64)
        for i in range(n):
            err = LIBASCOT.Efield_eval_e(
                e[i, :].ctypes.data_as(ct.POINTER(ct.c_double)),
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
                ct.byref(bobj),
            )
            if err:
                assert np.allclose(e[i, :], [0, 0, 0])
                raise ValueError(err)

        return e.T


@pytest.fixture
def helpers():
    return Helpers


@pytest.fixture
def bdata():
    rgrid = np.linspace(1.0, 3.0, 100) * unyt.m
    zgrid = np.linspace(-1.0, 1.0, 100) * unyt.m
    axisrz = np.array([2.0, 0.0]) * unyt.m

    R, Z = np.meshgrid(rgrid, zgrid)
    psi = ((R - axisrz[0]) ** 2 + (Z - axisrz[1]) ** 2).v * unyt.Wb/unyt.rad
    return {
        "rgrid": rgrid,
        "zgrid": zgrid,
        "psi": psi,
        "bphi": np.full(psi.shape, 1.0) * unyt.T,
        "axisrz": axisrz,
        "psilimits": (0.0, 1.0) * unyt.Wb/unyt.rad,
    }
