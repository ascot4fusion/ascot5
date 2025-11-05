import ctypes

import pytest
import numpy as np

from a5py.libascot import LIBASCOT, init_fun, Spline1D, Spline2D, Spline3D

init_fun(
    "solve_compact_cubic_spline",
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
)

init_fun(
    "Spline1D_init",
    ctypes.POINTER(Spline1D),
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    restype=ctypes.c_int32,
)

init_fun(
    "Spline2D_init",
    ctypes.POINTER(Spline2D),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    restype=ctypes.c_int32,
)

init_fun(
    "Spline3D_init",
    ctypes.POINTER(Spline3D),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    restype=ctypes.c_int32,
)

init_fun(
    "Spline1D_eval_f_df",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(Spline1D),
    ctypes.c_double,
    restype=ctypes.c_int32
)

init_fun(
    "Spline2D_eval_f_df",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(Spline2D),
    ctypes.c_double,
    ctypes.c_double,
    restype=ctypes.c_int32
)

init_fun(
    "Spline3D_eval_f_df",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(Spline3D),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    restype=ctypes.c_int32
)


def test_solve_compact_cubic_spline():
    err = []
    ngrid = np.array([10, 20, 40, 80, 160])
    for n in ngrid:
        x = np.linspace(-np.pi, np.pi, n+1)

        y = np.sin(x)
        c = np.full((x.size,2), 0.0, dtype=np.float64)

        LIBASCOT.solve_compact_cubic_spline(
            x.size,
            c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            0
            )

        dx = x[1] - x[0]
        c[:,1] /= dx**2
        d2ydx2 = -y

        err.append(np.max(np.abs(c[:,1] - d2ydx2)))

    err = np.array(err)
    a, b = err[:-1] / err[1:], (ngrid[1:] / ngrid[:-1])**2
    assert np.isclose(a, b, rtol=0.2).all()

@pytest.mark.parametrize(
    "periodic", [True, False]
)
def test_spline_1d(periodic):
    nq = 1000
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    xq = xlim[0] + np.random.rand(nq) * (xlim[1] - xlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])
    spline = Spline1D()

    fx0 = np.sin(xq)
    fx1 = np.cos(xq)
    fx2 = -np.sin(xq)
    err_f_df = np.zeros((3, ngrid.size))
    for i, n in enumerate(ngrid):
        if periodic:
            x = np.linspace(xlim[0], xlim[1], n+1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], n)
        y = np.sin(x)

        LIBASCOT.Spline1D_init(
            ctypes.byref(spline),
            n,
            1 if periodic else 0,
            xlim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )

        for j in range(nq):
            f_df = np.zeros(3)
            LIBASCOT.Spline1D_eval_f_df(
                f_df.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(spline),
                xq[j],
            )
            err_f_df[0,i] = max(err_f_df[0,i], abs(f_df[0] - fx0[j]))
            err_f_df[1,i] = max(err_f_df[1,i], abs(f_df[1] - fx1[j]))
            err_f_df[2,i] = max(err_f_df[2,i], abs(f_df[2] - fx2[j]))

    a = ngrid[1:] / ngrid[:-1]
    rate = err_f_df[:, :-1] / err_f_df[:, 1:]

    assert np.all(rate[0,:] > a**3), "f convergence too slow"
    assert np.all(rate[1,:] > a**2), "fx convergence too slow"
    assert np.all(rate[2,:] > a**1), "fxx convergence too slow"


@pytest.mark.parametrize(
    "periodic", [True, False]
)
def test_spline_2d(periodic):
    nq = 1000
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi/2, np.pi/2], dtype=np.float64)
    xq = xlim[0] + np.random.rand(nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + np.random.rand(nq) * (ylim[1] - ylim[0])
    ngrid = np.array([10, 20, 40, 80, 160])
    spline = Spline2D()

    fx0, fy0 = np.sin(xq), np.sin(2*yq)
    fx1, fy1 = np.cos(xq), 2 * np.cos(2*yq)
    fx2, fy2 = -np.sin(xq), -4 * np.sin(2*yq)
    err_f_df = np.zeros((6, ngrid.size))
    for i, n in enumerate(ngrid):
        nx, ny = n, n
        if periodic:
            x = np.linspace(xlim[0], xlim[1], nx+1)[:-1]
            y = np.linspace(ylim[0], ylim[1], ny+1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], nx)
            y = np.linspace(ylim[0], ylim[1], ny)

        X, Y = np.meshgrid(x, y, indexing="ij")
        f = np.sin(X) * np.sin(2 * Y)

        LIBASCOT.Spline2D_init(
            ctypes.byref(spline),
            nx, ny,
            1 if periodic else 0,
            1 if periodic else 0,
            xlim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ylim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            f.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        )

        for j in range(nq):
            f_df = np.zeros(6)
            LIBASCOT.Spline2D_eval_f_df(
                f_df.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(spline),
                xq[j], yq[j]
            )
            err_f_df[0,i] = max(err_f_df[0,i], abs(f_df[0] - (fx0*fy0)[j]))
            err_f_df[1,i] = max(err_f_df[1,i], abs(f_df[1] - (fx1*fy0)[j]))
            err_f_df[2,i] = max(err_f_df[2,i], abs(f_df[2] - (fx0*fy1)[j]))
            err_f_df[3,i] = max(err_f_df[3,i], abs(f_df[3] - (fx2*fy0)[j]))
            err_f_df[4,i] = max(err_f_df[4,i], abs(f_df[4] - (fx0*fy2)[j]))
            err_f_df[5,i] = max(err_f_df[5,i], abs(f_df[5] - (fx1*fy1)[j]))

    a = ngrid[1:] / ngrid[:-1]
    rate = err_f_df[:, :-1] / err_f_df[:, 1:]

    assert np.all(rate[0,:] > a**3), "f value convergence too slow"
    assert np.all(rate[1,:] > a**2), "fx convergence too slow"
    assert np.all(rate[2,:] > a**2), "fy convergence too slow"
    assert np.all(rate[3,:] > a**1), "fxx convergence too slow"
    assert np.all(rate[4,:] > a**1), "fyy convergence too slow"
    assert np.all(rate[5,:] > a**1), "fxy convergence too slow"


@pytest.mark.parametrize(
    "periodic", [True, False]
)
def test_spline_2d(periodic):
    nq = 1000
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi/2, np.pi/2], dtype=np.float64)
    zlim = np.array([-np.pi/3, np.pi/3], dtype=np.float64)
    xq = xlim[0] + np.random.rand(nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + np.random.rand(nq) * (ylim[1] - ylim[0])
    zq = zlim[0] + np.random.rand(nq) * (zlim[1] - zlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])
    spline = Spline3D()

    fx0, fy0, fz0 = np.sin(xq), np.sin(2*yq), np.sin(3*zq)
    fx1, fy1, fz1 = np.cos(xq), 2 * np.cos(2*yq), 3 * np.cos(3*zq)
    fx2, fy2, fz2 = -np.sin(xq), -4 * np.sin(2*yq), -9 * np.sin(3*zq)
    err_f_df = np.zeros((10, ngrid.size))
    for i, n in enumerate(ngrid):
        nx, ny, nz = n, n, n
        if periodic:
            x = np.linspace(xlim[0], xlim[1], nx+1)[:-1]
            y = np.linspace(ylim[0], ylim[1], ny+1)[:-1]
            z = np.linspace(zlim[0], zlim[1], nz+1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], nx)
            y = np.linspace(ylim[0], ylim[1], ny)
            z = np.linspace(zlim[0], zlim[1], nz)

        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        f = (np.sin(X) * np.sin(2 * Y) * np.sin(3 * Z))

        LIBASCOT.Spline3D_init(
            ctypes.byref(spline),
            nx, ny, nz,
            1 if periodic else 0,
            1 if periodic else 0,
            1 if periodic else 0,
            xlim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ylim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            zlim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            f.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        )

        for j in range(nq):
            f_df = np.zeros(10)
            LIBASCOT.Spline3D_eval_f_df(
                f_df.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(spline),
                xq[j], yq[j], zq[j]
            )
            err_f_df[0,i] = max(err_f_df[0,i], abs(f_df[0] - (fx0*fy0*fz0)[j]))
            err_f_df[1,i] = max(err_f_df[1,i], abs(f_df[1] - (fx1*fy0*fz0)[j]))
            err_f_df[2,i] = max(err_f_df[2,i], abs(f_df[2] - (fx0*fy1*fz0)[j]))
            err_f_df[3,i] = max(err_f_df[3,i], abs(f_df[3] - (fx0*fy0*fz1)[j]))
            err_f_df[4,i] = max(err_f_df[4,i], abs(f_df[4] - (fx2*fy0*fz0)[j]))
            err_f_df[5,i] = max(err_f_df[5,i], abs(f_df[5] - (fx0*fy2*fz0)[j]))
            err_f_df[6,i] = max(err_f_df[6,i], abs(f_df[6] - (fx0*fy0*fz2)[j]))
            err_f_df[7,i] = max(err_f_df[7,i], abs(f_df[7] - (fx1*fy1*fz0)[j]))
            err_f_df[8,i] = max(err_f_df[8,i], abs(f_df[8] - (fx1*fy0*fz1)[j]))
            err_f_df[9,i] = max(err_f_df[9,i], abs(f_df[9] - (fx0*fy1*fz1)[j]))

    a = ngrid[1:] / ngrid[:-1]
    rate = err_f_df[:, :-1] / err_f_df[:, 1:]

    assert np.all(rate[0,:] > a**3), "f value convergence too slow"
    assert np.all(rate[1,:] > a**2), "fx convergence too slow"
    assert np.all(rate[2,:] > a**2), "fy convergence too slow"
    assert np.all(rate[3,:] > a**2), "fz convergence too slow"
    assert np.all(rate[4,:] > a**1), "fxx convergence too slow"
    assert np.all(rate[5,:] > a**1), "fyy convergence too slow"
    assert np.all(rate[6,:] > a**1), "fzz convergence too slow"
    assert np.all(rate[7,:] > a**1), "fxy convergence too slow"
    assert np.all(rate[8,:] > a**1), "fxz convergence too slow"
    assert np.all(rate[9,:] > a**1), "fyz convergence too slow"
