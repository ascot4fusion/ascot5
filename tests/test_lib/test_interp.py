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


def _make_interpolator(scheme, dim, periodic, lims, data):
    """Initialize a linear (2D or 3D) or spline (1D, 2D or 3D) interpolator.

    Returns a function that evaluates the interpolant and its first derivatives
    at a point, asserts that the evaluation succeeded (or failed, with
    expect_error=True), and returns the output array.
    """
    struct = {
        ("Spline", 1): Spline1D,
        ("Linear", 2): Linear2D, ("Linear", 3): Linear3D,
        ("Spline", 2): Spline2D, ("Spline", 3): Spline3D,
    }[scheme, dim]()
    init = getattr(LIBASCOT, f"{scheme}{dim}D_init")
    eval_f_df = getattr(LIBASCOT, f"{scheme}{dim}D_eval_f_df")
    ptr = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    data = np.ascontiguousarray(data, dtype=np.float64)
    lims = [np.asarray(l, dtype=np.float64) for l in lims]
    err = init(
        ctypes.byref(struct), *data.shape, *(dim*[int(periodic)]),
        *(ptr(l) for l in lims), ptr(data),
    )
    assert not err

    def evaluate(*x, expect_error=False):
        out = np.zeros({1: 3, 2: 6, 3: 4}[dim])
        err = eval_f_df(ptr(out), ctypes.byref(struct), *x)
        assert bool(err) == expect_error, (x, err)
        return out

    # The linear interpolators don't copy the data, so it must stay alive
    evaluate.keepalive = (struct, data, lims)
    return evaluate


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
    shift = 2 if periodic else 0
    xq = xlim[0] + rng.uniform(-shift, 1 + shift, nq) * (xlim[1] - xlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])

    fx0 = np.sin(xq)
    fx1 = np.cos(xq)
    fx2 = -np.sin(xq)
    err_f_df = np.zeros((3, ngrid.size))
    for i, n in enumerate(ngrid):
        if periodic:
            x = np.linspace(xlim[0], xlim[1], n+1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], n)
        evaluate = _make_interpolator("Spline", 1, periodic, [xlim], np.sin(x))

        for j in range(nq):
            f_df = evaluate(xq[j])
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
    "scheme, min_order", [("Linear", 1), ("Spline", 3)],
    ids=["linear", "spline"],
)
def test_interpolate_2d(periodic, scheme, min_order):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi/2, np.pi/2], dtype=np.float64)
    shift = 2 if periodic else 0
    xq = xlim[0] + rng.uniform(-shift, 1 + shift, nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + rng.uniform(-shift, 1 + shift, nq) * (ylim[1] - ylim[0])
    ngrid = np.array([10, 20, 40, 80, 160])

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

        evaluate = _make_interpolator(scheme, 2, periodic, [xlim, ylim], f)

        for j in range(nq):
            f_df = evaluate(xq[j], yq[j])
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
    min_orders = min_order - np.array([0, 1, 1, 2, 2, 2])
    for k in range(6):
        assert np.all(rate_df[k,:] > a**min_orders[k]), \
            f"f_df: {labels[k]} convergence too slow"


@pytest.mark.parametrize(
    "periodic", [True, False], ids=["periodic", "natural"]
)
@pytest.mark.parametrize(
    "scheme, min_order", [("Linear", 1), ("Spline", 3)],
    ids=["linear", "spline"],
)
def test_interpolate_3d(periodic, scheme, min_order):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi/2, np.pi/2], dtype=np.float64)
    zlim = np.array([-np.pi/3, np.pi/3], dtype=np.float64)
    shift = 2 if periodic else 0
    xq = xlim[0] + rng.uniform(-shift, 1 + shift, nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + rng.uniform(-shift, 1 + shift, nq) * (ylim[1] - ylim[0])
    zq = zlim[0] + rng.uniform(-shift, 1 + shift, nq) * (zlim[1] - zlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])

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

        evaluate = _make_interpolator(
            scheme, 3, periodic, [xlim, ylim, zlim], f)

        for j in range(nq):
            f_df = evaluate(xq[j], yq[j], zq[j])
            err_f_df[0,i] = max(err_f_df[0,i], abs(f_df[0] - (fx0*fy0*fz0)[j]))
            err_f_df[1,i] = max(err_f_df[1,i], abs(f_df[1] - (fx1*fy0*fz0)[j]))
            err_f_df[2,i] = max(err_f_df[2,i], abs(f_df[2] - (fx0*fy1*fz0)[j]))
            err_f_df[3,i] = max(err_f_df[3,i], abs(f_df[3] - (fx0*fy0*fz1)[j]))

    a = ngrid[1:] / ngrid[:-1]
    rate_df = err_f_df[:, :-1] / err_f_df[:, 1:]

    labels = ["f value", "fx", "fy", "fz"]
    # Each derivative lowers the expected order of convergence by one
    min_orders = min_order - np.array([0, 1, 1, 1])
    for k in range(4):
        assert np.all(rate_df[k,:] > a**min_orders[k]), \
            f"f_df: {labels[k]} convergence too slow"


@pytest.mark.parametrize("dim", [2, 3], ids=["2d", "3d"])
@pytest.mark.parametrize("scheme", ["Linear", "Spline"])
def test_grid_edges(scheme, dim):
    """Test the interpolators on the edges of a grid with the natural boundary
    condition.

    The grids are chosen so that the naive cell index (x - xmin) / dx rounds
    wrong at the upper edge: on the first axis x = xlim[1] lands two cells
    back, and on the second axis one ulp below xlim[1] lands one cell past the
    last one.
    """
    lims = [[4.0, 8.5], [-4.0, 4.0], [-3.0, -1.0]][:dim]
    ngrid = [8, 3, 5][:dim]
    axes = [np.linspace(*l, n) for l, n in zip(lims, ngrid)]
    f = np.random.default_rng(6).random(ngrid)
    evaluate = _make_interpolator(scheme, dim, False, lims, f)

    # Both schemes interpolate, so the grid points recover the data
    for idx in np.ndindex(*ngrid):
        out = evaluate(*(a[i] for a, i in zip(axes, idx)))
        assert np.isclose(out[0], f[idx]), idx

    # One ulp below the upper corner along each axis
    for k in range(dim):
        point = [l[1] for l in lims]
        point[k] = np.nextafter(point[k], -np.inf)
        out = evaluate(*point)
        assert np.isclose(out[0], f[(-1,)*dim]), k

    # Just outside each edge, and NaN
    for k in range(dim):
        for x in [lims[k][0] - 1e-9, lims[k][1] + 1e-9, np.nan]:
            point = np.mean(lims, axis=1)
            point[k] = x
            evaluate(*point, expect_error=True)
