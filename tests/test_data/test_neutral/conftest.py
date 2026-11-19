import ctypes as ct
import pytest

import numpy as np
from a5py.libascot import LIBASCOT, init_fun
from a5py.data.neutral import Neutral

init_fun(
    "Neutral_eval_density",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_size_t,
    ct.POINTER(Neutral),
    restype=ct.c_int,
)

init_fun(
    "Neutral_eval_temperature",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_size_t,
    ct.POINTER(Neutral),
    restype=ct.c_int,
)

init_fun(
    "Neutral_eval_nT",
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Neutral),
    restype=ct.c_int,
)

init_fun(
    "Neutral_get_n_species",
    ct.POINTER(Neutral),
    restype=ct.c_size_t,
)

init_fun(
    "Neutral_get_species_znum",
    ct.POINTER(Neutral),
    restype=ct.POINTER(ct.c_int),
)

init_fun(
    "Neutral_get_species_anum",
    ct.POINTER(Neutral),
    restype=ct.POINTER(ct.c_int),
)


class Helpers:

    @staticmethod
    def neutral_eval_nT(neutral, rho, r, phi, z, t):
        obj = Neutral()
        obj.use(neutral)

        nspecies = LIBASCOT.Neutral_get_n_species(ct.byref(obj))

        npnt = r.size
        n = np.zeros((npnt, nspecies), dtype=np.float64)
        T = np.zeros((npnt, nspecies), dtype=np.float64)
        for i in range(npnt):
            err = 0
            err0 = LIBASCOT.Neutral_eval_nT(
                n[i, :].ctypes.data_as(ct.POINTER(ct.c_double)),
                T[i, :].ctypes.data_as(ct.POINTER(ct.c_double)),
                rho[i],
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
            )
            if err0:
                assert np.allclose(n[i, :], 1)
                assert np.allclose(T[i, :], 1)
                err = err0

            for j in range(nspecies):
                ni, Ti = np.array(
                    [
                        0.0,
                    ]
                ), np.array(
                    [
                        0.0,
                    ]
                )
                err0 = LIBASCOT.Neutral_eval_density(
                    ni.ctypes.data_as(ct.POINTER(ct.c_double)),
                    rho[i],
                    r[i],
                    phi[i],
                    z[i],
                    t[i],
                    j,
                    ct.byref(obj),
                )
                if err0:
                    assert np.allclose(ni, 1)
                    err = err0

                err0 = LIBASCOT.Neutral_eval_temperature(
                    Ti.ctypes.data_as(ct.POINTER(ct.c_double)),
                    rho[i],
                    r[i],
                    phi[i],
                    z[i],
                    t[i],
                    j,
                    ct.byref(obj),
                )
                if err0:
                    assert np.allclose(Ti, 1)
                    err = err0

                assert np.allclose(n[i, j], ni)
                assert np.allclose(T[i, j], Ti)

            if err:
                raise ValueError(err)

        return n, T

    @staticmethod
    def assert_species(neutral):
        obj = Neutral()
        obj.use(neutral)

        nspecies = LIBASCOT.Neutral_get_n_species(ct.byref(obj))
        anum = np.ctypeslib.as_array(
            LIBASCOT.Neutral_get_species_anum(ct.byref(obj)), shape=(nspecies,)
        )
        znum = np.ctypeslib.as_array(
            LIBASCOT.Neutral_get_species_znum(ct.byref(obj)), shape=(nspecies,)
        )

        return anum, znum


@pytest.fixture
def helpers():
    return Helpers
