import ctypes as ct
import pytest

import unyt
import numpy as np
from a5py.libascot import LIBASCOT, init_fun
from a5py.data.plasma import Plasma

init_fun(
    "Plasma_eval_density",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_size_t,
    ct.POINTER(Plasma),
    restype=ct.c_int,
)

init_fun(
    "Plasma_eval_temperature",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_size_t,
    ct.POINTER(Plasma),
    restype=ct.c_int,
)

init_fun(
    "Plasma_eval_nT",
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Plasma),
    restype=ct.c_int,
)

init_fun(
    "Plasma_eval_flow",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Plasma),
    restype=ct.c_int,
)

init_fun(
    "Plasma_get_n_species",
    ct.POINTER(Plasma),
    restype=ct.c_size_t,
)

init_fun(
    "Plasma_get_species_znum",
    ct.POINTER(Plasma),
    restype=ct.POINTER(ct.c_int),
)

init_fun(
    "Plasma_get_species_anum",
    ct.POINTER(Plasma),
    restype=ct.POINTER(ct.c_int),
)

init_fun(
    "Plasma_get_species_mass",
    ct.POINTER(Plasma),
    restype=ct.POINTER(ct.c_double),
)

init_fun(
    "Plasma_get_species_charge",
    ct.POINTER(Plasma),
    restype=ct.POINTER(ct.c_int),
)


class Helpers:

    @staticmethod
    def plasma_eval_nT(plasma, rho, r, phi, z, t):
        obj = Plasma()
        obj.use(plasma)

        nspecies = LIBASCOT.Plasma_get_n_species(ct.byref(obj))

        npnt = r.size
        n = np.zeros((npnt, nspecies), dtype=np.float64) * unyt.m**-3
        T = np.zeros((npnt, nspecies), dtype=np.float64) * unyt.J
        for i in range(npnt):
            err = 0
            err0 = LIBASCOT.Plasma_eval_nT(
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
                assert np.allclose(n[i, :], 1e20)
                assert np.allclose(T[i, :], (1e3*unyt.eV).to("J"))
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
                err0 = LIBASCOT.Plasma_eval_density(
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
                    assert np.allclose(ni, 1e20)
                    err = err0

                err0 = LIBASCOT.Plasma_eval_temperature(
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
                    assert np.allclose(Ti, (1e3*unyt.eV).to("J"))
                    err = err0

                assert np.allclose(n[i, j], ni)
                assert np.allclose(T[i, j], Ti)

            if err:
                raise ValueError(err)

        return n[:,0], T[:,0], n[:,1:], T[:,1]

    @staticmethod
    def plasma_eval_rotation(plasma, rho, r, phi, z, t):
        obj = Plasma()
        obj.use(plasma)

        npnt = r.size
        flow = np.zeros((npnt,), dtype=np.float64)
        for i in range(npnt):
            err = LIBASCOT.Plasma_eval_flow(
                flow[i:i+1].ctypes.data_as(ct.POINTER(ct.c_double)),
                rho[i],
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
            )
            if err:
                assert np.allclose(flow[i,], 0)
                raise ValueError(err)

        return flow/r.v

    @staticmethod
    def assert_species(plasma):
        obj = Plasma()
        obj.use(plasma)

        nspecies = LIBASCOT.Plasma_get_n_species(ct.byref(obj))
        anum = np.ctypeslib.as_array(
            LIBASCOT.Plasma_get_species_anum(ct.byref(obj)), shape=(nspecies-1,)
        )
        znum = np.ctypeslib.as_array(
            LIBASCOT.Plasma_get_species_znum(ct.byref(obj)), shape=(nspecies-1,)
        )
        mass = np.ctypeslib.as_array(
            LIBASCOT.Plasma_get_species_mass(ct.byref(obj)), shape=(nspecies,)
        )
        charge = np.ctypeslib.as_array(
            LIBASCOT.Plasma_get_species_charge(ct.byref(obj)), shape=(nspecies,)
        )

        return anum, znum, charge, mass


@pytest.fixture
def helpers():
    return Helpers
