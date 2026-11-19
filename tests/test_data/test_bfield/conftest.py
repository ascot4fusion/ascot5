import ctypes as ct
import pytest

import numpy as np
from a5py.libascot import LIBASCOT, init_fun
from a5py.data.bfield import Bfield

init_fun(
    "Bfield_eval_psi",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Bfield),
    restype=ct.c_int,
)

init_fun(
    "Bfield_eval_psi_dpsi",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Bfield),
    restype=ct.c_int,
)

init_fun(
    "Bfield_eval_rho",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.POINTER(Bfield),
    restype=ct.c_int,
)

init_fun(
    "Bfield_eval_rho_drho",
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double),
    ct.POINTER(Bfield),
    restype=ct.c_int,
)

init_fun(
    "Bfield_eval_b",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Bfield),
    restype=ct.c_int,
)

init_fun(
    "Bfield_eval_b_db",
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.c_double,
    ct.POINTER(Bfield),
    restype=ct.c_int,
)

init_fun(
    "Bfield_eval_axis_rz",
    ct.POINTER(ct.c_double),
    ct.POINTER(Bfield),
    ct.c_double,
    restype=ct.c_int,
)


class Helpers:

    @staticmethod
    def bfield_eval_psi_dpsi_and_rho_drho(bfield, r, phi, z, t):
        obj = Bfield()
        obj.use(bfield)

        n = r.size
        psidpsi, rhodrho = np.zeros((n, 4)), np.zeros((n, 4))
        for i in range(n):
            err = 0
            psi, rho = np.zeros(1), np.zeros(2)
            err0 = LIBASCOT.Bfield_eval_psi(
                psi.ctypes.data_as(ct.POINTER(ct.c_double)),
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
            )
            if err0:
                assert np.allclose(psi, 1)
                err = ValueError(err0)

            err0 = LIBASCOT.Bfield_eval_rho(
                rho.ctypes.data_as(ct.POINTER(ct.c_double)),
                psi[0],
                ct.byref(obj),
            )
            if err0:
                assert np.allclose(rho, [1, 0])
                err = ValueError(err0)

            err0 = LIBASCOT.Bfield_eval_psi_dpsi(
                psidpsi[i, :].ctypes.data_as(ct.POINTER(ct.c_double)),
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
            )
            if err or err0:
                assert np.allclose(psidpsi[i, :], [1, 0, 0, 0])
                err = ValueError(err0)

            err0 = LIBASCOT.Bfield_eval_rho_drho(
                rhodrho[i, :].ctypes.data_as(ct.POINTER(ct.c_double)),
                psidpsi[i, :].ctypes.data_as(ct.POINTER(ct.c_double)),
                ct.byref(obj),
            )
            if err0:
                assert np.allclose(rhodrho[i, :], [1, 0, 0, 0])
                err = ValueError(err0)

            if err:
                raise err

            assert np.isclose(psi[0], psidpsi[i, 0])
            assert np.isclose(rho[0], rhodrho[i, 0])
            assert np.isclose(rhodrho[i, 1], rho[1] * psidpsi[i, 1])
            assert np.isclose(rhodrho[i, 2], rho[1] * psidpsi[i, 2])
            assert np.isclose(rhodrho[i, 3], rho[1] * psidpsi[i, 3])

        return psidpsi.T, rhodrho.T[::4, :].ravel()

    @staticmethod
    def bfield_eval_b_db(bfield, r, phi, z, t):
        obj = Bfield()
        obj.use(bfield)

        n = r.size
        b, jac = np.zeros((3, n)), np.zeros((12, n))
        for i in range(n):
            err = 0
            bvec = np.zeros(3)
            err0 = LIBASCOT.Bfield_eval_b(
                bvec.ctypes.data_as(ct.POINTER(ct.c_double)),
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
            )
            if err or err0:
                assert np.allclose(bvec, [1., 1., 1.])
                err = ValueError(err0)

            bdb = np.zeros(15)
            err0 = LIBASCOT.Bfield_eval_b_db(
                bdb.ctypes.data_as(ct.POINTER(ct.c_double)),
                r[i],
                phi[i],
                z[i],
                t[i],
                ct.byref(obj),
            )
            if err or err0:
                assert np.allclose(bdb[:3], [1., 1., 1.])
                assert np.allclose(bdb[3:], np.zeros(12))
                err = ValueError(err0)

            if err:
                raise err

            b[:, i] = bdb[:3]
            jac[:, i] = bdb[3:]
            assert np.allclose(bvec, b[:, i])

        return b, jac

    @staticmethod
    def bfield_eval_axis_rz(bfield, phi):
        obj = Bfield()
        obj.use(bfield)

        n = phi.size
        r, z = np.zeros((1, n)), np.zeros((1, n))
        for i in range(n):
            rz = np.zeros(2)
            err = LIBASCOT.Bfield_eval_axis_rz(
                rz.ctypes.data_as(ct.POINTER(ct.c_double)),
                ct.byref(obj),
                phi[i],
            )
            if err:
                assert np.allclose(rz, [1., 0.])
                raise ValueError(err)

            r[0, i] = rz[0]
            z[0, i] = rz[1]

        return r, z

    @staticmethod
    def assert_gradient_of_psi(bfield, r, dr, phi, dphi, z, dz):
        t = np.ones(r.size)
        psidpsi0, _ = Helpers.bfield_eval_psi_dpsi_and_rho_drho(
            bfield,
            r,
            phi.to("rad").v,
            z,
            t,
        )

        psi, _ = Helpers.bfield_eval_psi_dpsi_and_rho_drho(
            bfield,
            r + dr,
            phi.to("rad").v,
            z,
            t,
        )
        assert np.allclose(psi[0, :] - psidpsi0[0, :], psidpsi0[1, :] * dr, atol=1e-3)

        psi, _ = Helpers.bfield_eval_psi_dpsi_and_rho_drho(
            bfield,
            r,
            (phi + dphi).to("rad").v,
            z,
            t,
        )
        assert np.allclose(psi[0, :] - psidpsi0[0, :], psidpsi0[2, :] * dphi, atol=1e-3)

        psi, _ = Helpers.bfield_eval_psi_dpsi_and_rho_drho(
            bfield,
            r,
            phi.to("rad").v,
            z + dz,
            t,
        )
        assert np.allclose(psi[0, :] - psidpsi0[0, :], psidpsi0[3, :] * dz, atol=1e-3)

    @staticmethod
    def assert_gradient_of_b(bfield, r, dr, phi, dphi, z, dz):
        t = np.ones(r.size)
        b, jac = Helpers.bfield_eval_b_db(
            bfield,
            r,
            phi.to("rad").v,
            z,
            t,
        )
        br0, bphi0, bz0 = b[0, :], b[1, :], b[2, :]
        jac = jac[:9, :].reshape((3, 3, r.size))

        b, _ = Helpers.bfield_eval_b_db(
            bfield,
            r + dr,
            phi.to("rad").v,
            z,
            t,
        )
        br, bphi, bz = b[0, :], b[1, :], b[2, :]
        assert np.allclose(br - br0, jac[0, 0] * dr, atol=1e-3)
        assert np.allclose(bphi - bphi0, jac[1, 0] * dr, atol=1e-3)
        assert np.allclose(bz - bz0, jac[2, 0] * dr, atol=1e-3)

        b, _ = Helpers.bfield_eval_b_db(
            bfield,
            r,
            (phi + dphi).to("rad").v,
            z,
            t,
        )
        br, bphi, bz = b[0, :], b[1, :], b[2, :]
        dphi = dphi.to("rad").v
        assert np.allclose(br - br0, jac[0, 1] * dphi, atol=1e-3)
        assert np.allclose(bphi - bphi0, jac[1, 1] * dphi, atol=1e-3)
        assert np.allclose(bz - bz0, jac[2, 1] * dphi, atol=1e-3)

        b, _ = Helpers.bfield_eval_b_db(
            bfield,
            r,
            phi.to("rad").v,
            z + dz,
            t,
        )
        br, bphi, bz = b[0, :], b[1, :], b[2, :]
        assert np.allclose(br - br0, jac[0, 2] * dz, atol=1e-3)
        assert np.allclose(bphi - bphi0, jac[1, 2] * dz, atol=1e-3)
        assert np.allclose(bz - bz0, jac[2, 2] * dz, atol=1e-3)


@pytest.fixture
def helpers():
    return Helpers
