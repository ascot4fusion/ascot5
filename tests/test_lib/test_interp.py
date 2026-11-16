import ctypes

import pytest
import numpy as np

from a5py.libascot import (
    LIBASCOT, init_fun, Linear2D, Linear3D, Spline1D, Spline2D, Spline3D,
)

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
    "Spline1D_eval_f_df",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(Spline1D),
    ctypes.c_double,
    restype=ctypes.c_int32
)

for (scheme, ty2D, ty3D) in [("Linear", Linear2D, Linear3D), ("Spline", Spline2D, Spline3D)]:
    init_fun(
        f"{scheme}2D_init",
        ctypes.POINTER(ty2D),
        ctypes.c_size_t,
        ctypes.c_size_t,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        # The linear init functions return void
        restype=None if scheme == "Linear" else ctypes.c_int32,
    )

    init_fun(
        f"{scheme}2D_eval_f_df",
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ty2D),
        ctypes.c_double,
        ctypes.c_double,
        restype=ctypes.c_int32
    )

    init_fun(
        f"{scheme}3D_init",
        ctypes.POINTER(ty3D),
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
        # The linear init functions return void
        restype=None if scheme == "Linear" else ctypes.c_int32,
    )

    init_fun(
        f"{scheme}3D_eval_f_df",
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ty3D),
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
    "periodic", [True, False], ids=["periodic", "natural"]
)
def test_spline_1d(periodic):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    xq = xlim[0] + rng.random(nq) * (xlim[1] - xlim[0])
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
    "periodic", [True, False], ids=["periodic", "natural"]
)
@pytest.mark.parametrize(
    "interpolator", [
        pytest.param(
            {
                "struct": Linear2D,
                "init": "Linear2D_init",
                "eval_f_df": "Linear2D_eval_f_df",
                "min_order": 1,
            },
            id="linear",
        ),
        pytest.param(
            {
                "struct": Spline2D,
                "init": "Spline2D_init",
                "eval_f_df": "Spline2D_eval_f_df",
                "min_order": 3,
            },
            id="spline",
        ),
    ],
)
def test_interpolate_2d(periodic, interpolator):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi/2, np.pi/2], dtype=np.float64)
    xq = xlim[0] + rng.random(nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + rng.random(nq) * (ylim[1] - ylim[0])
    ngrid = np.array([10, 20, 40, 80, 160])
    interp = interpolator["struct"]()
    init = getattr(LIBASCOT, interpolator["init"])
    eval_f_df = getattr(LIBASCOT, interpolator["eval_f_df"])

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

        init(
            ctypes.byref(interp),
            nx, ny,
            1 if periodic else 0,
            1 if periodic else 0,
            xlim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ylim.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            f.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        )

        for j in range(nq):
            f_df = np.zeros(6)
            eval_f_df(
                f_df.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(interp),
                xq[j], yq[j]
            )
            err_f_df[0,i] = max(err_f_df[0,i], abs(f_df[0] - (fx0*fy0)[j]))
            err_f_df[1,i] = max(err_f_df[1,i], abs(f_df[1] - (fx1*fy0)[j]))
            err_f_df[2,i] = max(err_f_df[2,i], abs(f_df[2] - (fx0*fy1)[j]))
            err_f_df[3,i] = max(err_f_df[3,i], abs(f_df[3] - (fx2*fy0)[j]))
            err_f_df[4,i] = max(err_f_df[4,i], abs(f_df[4] - (fx0*fy2)[j]))
            err_f_df[5,i] = max(err_f_df[5,i], abs(f_df[5] - (fx1*fy1)[j]))

    a = ngrid[1:] / ngrid[:-1]
    rate_df = err_f_df[:, :-1] / err_f_df[:, 1:]

    labels = ["f value", "fx", "fy", "fxx", "fyy", "fxy"]
    # Each derivative lowers the expected order of convergence by one
    min_orders = interpolator["min_order"] - np.array([0, 1, 1, 2, 2, 2])
    for k in range(6):
        assert np.all(rate_df[k,:] > a**min_orders[k]), \
            f"f_df: {labels[k]} convergence too slow"


@pytest.mark.parametrize(
    "periodic", [True, False], ids=["periodic", "natural"]
)
@pytest.mark.parametrize(
    "interpolator", [
        pytest.param(
            {
                "struct": Linear3D,
                "init": "Linear3D_init",
                "eval_f_df": "Linear3D_eval_f_df",
                "min_order": 1,
            },
            id="linear",
        ),
        pytest.param(
            {
                "struct": Spline3D,
                "init": "Spline3D_init",
                "eval_f_df": "Spline3D_eval_f_df",
                "min_order": 3,
            },
            id="spline",
        ),
    ],
)
def test_interpolate_3d(periodic, interpolator):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi/2, np.pi/2], dtype=np.float64)
    zlim = np.array([-np.pi/3, np.pi/3], dtype=np.float64)
    xq = xlim[0] + rng.random(nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + rng.random(nq) * (ylim[1] - ylim[0])
    zq = zlim[0] + rng.random(nq) * (zlim[1] - zlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])
    interp = interpolator["struct"]()
    init = getattr(LIBASCOT, interpolator["init"])
    eval_f_df = getattr(LIBASCOT, interpolator["eval_f_df"])

    fx0, fy0, fz0 = np.sin(xq), np.sin(2*yq), np.sin(3*zq)
    fx1, fy1, fz1 = np.cos(xq), 2 * np.cos(2*yq), 3 * np.cos(3*zq)
    err_f_df = np.zeros((4, ngrid.size))
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

        init(
            ctypes.byref(interp),
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
            f_df = np.zeros(4)
            eval_f_df(
                f_df.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(interp),
                xq[j], yq[j], zq[j]
            )
            err_f_df[0,i] = max(err_f_df[0,i], abs(f_df[0] - (fx0*fy0*fz0)[j]))
            err_f_df[1,i] = max(err_f_df[1,i], abs(f_df[1] - (fx1*fy0*fz0)[j]))
            err_f_df[2,i] = max(err_f_df[2,i], abs(f_df[2] - (fx0*fy1*fz0)[j]))
            err_f_df[3,i] = max(err_f_df[3,i], abs(f_df[3] - (fx0*fy0*fz1)[j]))

    a = ngrid[1:] / ngrid[:-1]
    rate_df = err_f_df[:, :-1] / err_f_df[:, 1:]

    labels = ["f value", "fx", "fy", "fz"]
    # Each derivative lowers the expected order of convergence by one
    min_orders = interpolator["min_order"] - np.array([0, 1, 1, 1])
    for k in range(4):
        assert np.all(rate_df[k,:] > a**min_orders[k]), \
            f"f_df: {labels[k]} convergence too slow"
