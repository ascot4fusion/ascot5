import ctypes

import pytest
import numpy as np

from a5py import libascot
from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "solve_compact_cubic_spline",
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
)

schemes_and_ranks = [
    ("Linear1D", ("f")),
    ("Linear2D", ("f", "f_df")),
    ("Linear3D", ("f", "f_df")),
    ("Spline1D", ("f", "f_ddf")),
    ("Spline2D", ("f", "f_ddf")),
    ("Spline3D", ("f", "f_df", "f_ddf")),
]

for scheme, ranks in schemes_and_ranks:

    for dim in range(1,4):
        if f"{dim}D" in scheme:
            break

    init_fun(
        f"{scheme}_init",
        ctypes.POINTER(getattr(libascot, scheme)),
        *(dim * [ctypes.c_size_t]),
        *(dim * [ctypes.c_int]),
        *(dim * [ctypes.POINTER(ctypes.c_double)]),
        ctypes.POINTER(ctypes.c_double),
        # The linear init functions return void
        restype=None if "Linear" in scheme else ctypes.c_int32,
    )

    for rank in ranks:
        init_fun(
            f"{scheme}_eval_{rank}",
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(getattr(libascot, scheme)),
            *(dim * [ctypes.c_double]),
            restype=ctypes.c_int32,
        )


def make_interpolator(scheme, periodic, lims, data, rank):
    """Initialize an interpolator and return its evaluation function.

    Supports linear interpolation in 2D and 3D, and spline interpolation
    in 1D, 2D and 3D.

    The returned function evaluates the interpolant and its first
    derivatives at a point. Set ``expect_error=True`` to assert that
    evaluation fails.
    """
    struct = getattr(libascot, scheme)()
    for dim in range(1,4):
        if f"{dim}D" in scheme:
            break

    init = getattr(LIBASCOT, f"{scheme}_init")
    evaluate = getattr(LIBASCOT, f"{scheme}_eval_{rank}")

    def as_ptr(array):
        return array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    data = np.ascontiguousarray(data, dtype=np.float64)
    lims = [np.asarray(limit, dtype=np.float64) for limit in lims]
    err = init(
        ctypes.byref(struct),
        *data.shape,
        *([int(periodic)] * dim),
        *(as_ptr(limit) for limit in lims),
        as_ptr(data),
    )
    assert not err, f"Interpolator initialization failed: {err}"

    output_size = 1
    if rank == "f_df":
        output_size += dim
    elif rank == "f_ddf":
        output_size += dim * 2
        if dim == 2:
            output_size += 1
        elif dim == 3:
            output_size += 3

    def interpolate(*x, expect_error=False):
        out = np.zeros(output_size)
        err = evaluate(as_ptr(out), ctypes.byref(struct), *x)

        assert bool(err) == expect_error, (x, err)
        return out

    # Ensure linear interpolators retain pointers to the input data.
    interpolate.keepalive = (struct, data, lims)

    return interpolate


def test_solve_compact_cubic_spline():
    err = []
    ngrid = np.array([10, 20, 40, 80, 160])
    for n in ngrid:
        x = np.linspace(-np.pi, np.pi, n + 1)

        y = np.sin(x)
        c = np.full((x.size, 2), 0.0, dtype=np.float64)

        LIBASCOT.solve_compact_cubic_spline(
            x.size,
            c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            0,
        )

        dx = x[1] - x[0]
        c[:, 1] /= dx**2
        d2ydx2 = -y

        err.append(np.max(np.abs(c[:, 1] - d2ydx2)))

    err = np.array(err)
    a, b = err[:-1] / err[1:], (ngrid[1:] / ngrid[:-1]) ** 2
    assert np.isclose(a, b, rtol=0.2).all()


@pytest.mark.parametrize("periodic", [True, False], ids=["periodic", "natural"])
@pytest.mark.parametrize(
    "scheme, rank",
    [
        ("Spline1D", "f"),
        ("Spline1D", "f_ddf"),
        ("Linear1D", "f"),
    ],
)
def test_spline_1d(scheme, rank, periodic):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    shift = 2 if periodic else 0
    xq = xlim[0] + rng.uniform(-shift, 1 + shift, nq) * (xlim[1] - xlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])

    min_order = 3 if "Spline" in scheme else 1

    labels = ["f value", "fx", "fxx"]
    f0 = (np.sin(xq), np.cos(xq), -np.sin(xq))
    if rank == "f":
        err_f = np.zeros((1, ngrid.size))
    elif rank == "f_ddf":
        err_f = np.zeros((3, ngrid.size))

    for i, n in enumerate(ngrid):
        if periodic:
            x = np.linspace(xlim[0], xlim[1], n + 1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], n)
        evaluate = make_interpolator(scheme, periodic, [xlim], np.sin(x), rank)

        for j in range(nq):
            f = evaluate(xq[j])
            for k in range(err_f.shape[0]):
                err_f[k, i] = max(err_f[k, i], abs(f[k] - f0[k][j]))

    a = ngrid[1:] / ngrid[:-1]
    rate = err_f[:, :-1] / err_f[:, 1:]

    min_orders = min_order - np.array([0, 1, 2])
    for k in range(err_f.shape[0]):
        assert np.all(
            rate[0, :] > a ** min_orders[k]
        ), f"{labels[k]} convergence too slow"


@pytest.mark.parametrize("periodic", [True, False], ids=["periodic", "natural"])
@pytest.mark.parametrize(
    "scheme, rank",
    [("Linear2D", "f"), ("Linear2D", "f_df"), ("Spline2D", "f"), ("Spline2D", "f_ddf")],
)
def test_interpolate_2d(scheme, rank, periodic):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi / 2, np.pi / 2], dtype=np.float64)
    shift = 2 if periodic else 0
    xq = xlim[0] + rng.uniform(-shift, 1 + shift, nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + rng.uniform(-shift, 1 + shift, nq) * (ylim[1] - ylim[0])
    ngrid = np.array([10, 20, 40, 80, 160])

    fx0, fy0 = np.sin(xq), np.sin(2 * yq)
    fx1, fy1 = np.cos(xq), 2 * np.cos(2 * yq)
    fx2, fy2 = -np.sin(xq), -4 * np.sin(2 * yq)

    min_order = 3 if "Spline" in scheme else 1

    # Each derivative lowers the expected order of convergence by one
    min_orders = min_order - np.array([0, 1, 1, 2, 2, 2])
    labels = ["f value", "fx", "fy", "fxx", "fyy", "fxy"]
    f0 = (fx0 * fy0, fx1 * fy0, fx0 * fy1, fx2 * fy0, fx0 * fy2, fx1 * fy1)

    if rank == "f":
        err_f = np.zeros((1, ngrid.size))
    elif rank == "f_df":
        err_f = np.zeros((3, ngrid.size))
    elif rank == "f_ddf":
        err_f = np.zeros((6, ngrid.size))

    for i, n in enumerate(ngrid):
        nx, ny = n, n
        if periodic:
            x = np.linspace(xlim[0], xlim[1], nx + 1)[:-1]
            y = np.linspace(ylim[0], ylim[1], ny + 1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], nx)
            y = np.linspace(ylim[0], ylim[1], ny)

        X, Y = np.meshgrid(x, y, indexing="ij")
        f = np.sin(X) * np.sin(2 * Y)

        evaluate = make_interpolator(scheme, periodic, [xlim, ylim], f, rank)

        for j in range(nq):
            f = evaluate(xq[j], yq[j])
            for k in range(err_f.shape[0]):
                err_f[k, i] = max(err_f[k, i], abs(f[k] - f0[k][j]))

    a = ngrid[1:] / ngrid[:-1]
    rate = err_f[:, :-1] / err_f[:, 1:]

    for k in range(err_f.shape[0]):
        assert np.all(
            rate[0, :] > a ** min_orders[k]
        ), f"{labels[k]} convergence too slow"


@pytest.mark.parametrize("periodic", [True, False], ids=["periodic", "natural"])
@pytest.mark.parametrize(
    "scheme, rank",
    [("Linear3D", "f"), ("Linear3D", "f_df"), ("Spline3D", "f"), ("Spline3D", "f_df"), ("Spline3D", "f_ddf")],
)
def test_interpolate_3d(scheme, rank, periodic):
    nq = 1000
    rng = np.random.default_rng(123)
    xlim = np.array([-np.pi, np.pi], dtype=np.float64)
    ylim = np.array([-np.pi / 2, np.pi / 2], dtype=np.float64)
    zlim = np.array([-np.pi / 3, np.pi / 3], dtype=np.float64)
    shift = 2 if periodic else 0
    xq = xlim[0] + rng.uniform(-shift, 1 + shift, nq) * (xlim[1] - xlim[0])
    yq = ylim[0] + rng.uniform(-shift, 1 + shift, nq) * (ylim[1] - ylim[0])
    zq = zlim[0] + rng.uniform(-shift, 1 + shift, nq) * (zlim[1] - zlim[0])
    ngrid = np.array([10, 20, 40, 80, 160])

    fx0, fy0, fz0 = np.sin(xq), np.sin(2 * yq), np.sin(3 * zq)
    fx1, fy1, fz1 = np.cos(xq), 2 * np.cos(2 * yq), 3 * np.cos(3 * zq)
    fx2, fy2, fz2 = -np.sin(xq), -4 * np.sin(2 * yq), -9 * np.sin(3 * zq)

    min_order = 3 if "Spline" in scheme else 1

    # Each derivative lowers the expected order of convergence by one
    min_orders = min_order - np.array([0, 1, 1, 1, 2, 2, 2, 2, 2, 2])
    labels = ["f value", "fx", "fy", "fz", "fxx", "fyy", "fzz", "fxy", "fxz", "fyz"]
    f0 = (fx0 * fy0 * fz0, fx1 * fy0 * fz0, fx0 * fy1 * fz0, fx0 * fy0 * fz1,
          fx2 * fy0 * fz0, fx0 * fy2 * fz0, fx0 * fy0 * fz2,
          fx1 * fy1 * fz0, fx1 * fy0 * fz1, fx0 * fy1 * fz1,
          )

    if rank == "f":
        err_f = np.zeros((1, ngrid.size))
    elif rank == "f_df":
        err_f = np.zeros((4, ngrid.size))
    elif rank == "f_ddf":
        err_f = np.zeros((10, ngrid.size))

    for i, n in enumerate(ngrid):
        nx, ny, nz = n, n, n
        if periodic:
            x = np.linspace(xlim[0], xlim[1], nx + 1)[:-1]
            y = np.linspace(ylim[0], ylim[1], ny + 1)[:-1]
            z = np.linspace(zlim[0], zlim[1], nz + 1)[:-1]
        else:
            x = np.linspace(xlim[0], xlim[1], nx)
            y = np.linspace(ylim[0], ylim[1], ny)
            z = np.linspace(zlim[0], zlim[1], nz)

        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        f = np.sin(X) * np.sin(2 * Y) * np.sin(3 * Z)

        evaluate = make_interpolator(scheme, periodic, [xlim, ylim, zlim], f, rank)
        for j in range(nq):
            f = evaluate(xq[j], yq[j], zq[j])
            for k in range(err_f.shape[0]):
                err_f[k, i] = max(err_f[k, i], abs(f[k] - f0[k][j]))


    a = ngrid[1:] / ngrid[:-1]
    rate_df = err_f[:, :-1] / err_f[:, 1:]

    for k in range(err_f.shape[0]):
        assert np.all(
            rate_df[k, :] > a ** min_orders[k]
        ), f"{labels[k]} convergence too slow"


@pytest.mark.parametrize(
        "scheme, rank",
        [
            ("Linear1D", "f"),
            ("Linear2D", "f"),
            ("Linear2D", "f_df"),
            ("Linear3D", "f"),
            ("Linear3D", "f_df"),
            ("Spline1D", "f"),
            ("Spline1D", "f_ddf"),
            ("Spline2D", "f"),
            ("Spline2D", "f_ddf"),
            ("Spline3D", "f"),
            ("Spline3D", "f_df"),
            ("Spline3D", "f_ddf"),
            ])
def test_grid_edges(scheme, rank):
    """Test the interpolators on the edges of a grid with the natural boundary
    condition.

    The grids are chosen so that the naive cell index (x - xmin) / dx rounds
    wrong at the upper edge: on the first axis x = xlim[1] lands two cells
    back, and on the second axis one ulp below xlim[1] lands one cell past the
    last one.
    """
    for dim in range(1,4):
        if f"{dim}D" in scheme:
            break

    lims = [[4.0, 8.5], [-4.0, 4.0], [-3.0, -1.0]][:dim]
    ngrid = [8, 3, 5][:dim]
    axes = [np.linspace(*l, n) for l, n in zip(lims, ngrid)]
    f = np.random.default_rng(6).random(ngrid)
    evaluate = make_interpolator(scheme, False, lims, f, rank)

    # Both schemes interpolate, so the grid points recover the data
    for idx in np.ndindex(*ngrid):
        out = evaluate(*(a[i] for a, i in zip(axes, idx)))
        assert np.isclose(out[0], f[idx]), idx

    # One ulp below the upper corner along each axis
    for k in range(dim):
        point = [l[1] for l in lims]
        point[k] = np.nextafter(point[k], -np.inf)
        out = evaluate(*point)
        assert np.isclose(out[0], f[(-1,) * dim]), k

    # Just outside each edge, and NaN
    for k in range(dim):
        for x in [lims[k][0] - 1e-9, lims[k][1] + 1e-9, np.nan]:
            point = np.mean(lims, axis=1)
            point[k] = x
            evaluate(*point, expect_error=True)
