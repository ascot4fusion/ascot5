import ctypes
import pytest

import numpy as np

from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "math_jac_rpz2xyz",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_double,
)

init_fun(
    "math_jac_xyz2rpz",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_double,
)

init_fun(
    "math_test_eval_vector_operations",
    *(10 * [ctypes.POINTER(ctypes.c_double)]),
)

init_fun(
    "math_test_eval_vector_transformations",
    *(8 * [ctypes.POINTER(ctypes.c_double)]),
)

init_fun(
    "math_test_eval_bin_index",
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_size_t),
)

init_fun(
    "math_test_eval_fmod",
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
)

init_fun(
    "math_test_eval_iabs",
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
)

init_fun(
    "math_crossed_plane",
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    restype=ctypes.c_double,
)

init_fun(
    "math_point_in_polygon",
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_double,
    restype=ctypes.c_int32,
)


@pytest.mark.parametrize(
    "a,expected",
    [
        (1, 1),
        (-1, 1),
        (0, 0),
    ],
)
def test_math_iabs(a, expected):
    b = np.empty(1, dtype=np.int32)
    LIBASCOT.math_test_eval_iabs(a, b.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    assert b[0] == expected


@pytest.mark.parametrize(
    "a,b,expected",
    [
        (5.0, 2.0, 1.0),
        (10.0, 3.0, 1.0),
        (7.5, 2.5, 0.0),
        (1.0, 2.0, 1.0),
        (0.0, 5.0, 0.0),
        (-5.0, 2.0, 1.0),
        (5.0, -2.0, -1.0),
        (-5.0, -2.0, -1.0),
    ],
)
def test_mathlib_eval_fmod(a, b, expected):
    out = np.empty(1, dtype=np.float64)
    LIBASCOT.math_test_eval_fmod(
        a, b, out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    assert out[0] == pytest.approx(expected)
    assert np.isclose(out[0], np.mod(a, b))


@pytest.mark.parametrize(
    "alpha, beta, gamma, k",
    [
        (0.0, 1.0, 0.1, 0.1),
        (0.0, 1.0, 1.1, -1.0),
        (-0.5, 0.5, 0.1, 0.6),
        (0.5, -0.5, 0.1, 0.4),
    ],
)
def test_mathlib_crossed_plane(alpha, beta, gamma, k):
    k0 = LIBASCOT.math_crossed_plane(alpha, beta, gamma)
    assert k0 == pytest.approx(k)


@pytest.mark.parametrize(
    "x, y, expected",
    [
        (0.9, 0.9, 1),
        (0.6, 0.5, 1),
        (0.4, 0.5, 0),
        (1.1, 0.5, 0),
        (1.0, 0.5, 0),
    ],
)
def test_mathlib_point_in_polygon(x, y, expected):
    vx = np.array([0.0, 1.0, 1.0, 0.0, 0.5])
    vy = np.array([0.0, 0.0, 1.0, 1.0, 0.5])
    out = LIBASCOT.math_point_in_polygon(
        vx.size,
        vx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        vy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        x,
        y,
    )
    assert out == expected

    vx = np.array([0.0, 0.5, 0.0, 1.0, 1.0])
    vy = np.array([0.0, 0.5, 1.0, 1.0, 0.0])
    out = LIBASCOT.math_point_in_polygon(
        vx.size,
        vx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        vy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        x,
        y,
    )
    assert out == expected


@pytest.mark.parametrize(
    "x, expected",
    [
        (0.3, 1),
        (0.0, 0),
        (1.0, 4),
    ],
)
def test_mathlib_eval_bin_index(x, expected):
    nx, xmin, xmax = 5, 0.0, 1.0
    val = np.array([0], dtype=np.uint64)
    LIBASCOT.math_test_eval_bin_index(
        nx,
        xmin,
        xmax,
        x,
        val.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t)),
    )
    assert val[0] == expected


def test_mathlib_vector_operations():
    rng = np.random.default_rng(0)

    for _ in range(100):
        a = rng.standard_normal(3)
        b = rng.standard_normal(3)
        c = rng.standard_normal(3)

        dot = np.empty(1, dtype=np.float64)
        cross = np.empty(3, dtype=np.float64)
        triple = np.empty(1, dtype=np.float64)
        det = np.empty(1, dtype=np.float64)
        norm = np.empty(1, dtype=np.float64)
        normc = np.empty(1, dtype=np.float64)
        unit = np.empty(3, dtype=np.float64)

        LIBASCOT.math_test_eval_vector_operations(
            a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            dot.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cross.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            triple.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            det.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            norm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            normc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            unit.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )

        assert np.isclose(dot[0], np.dot(a, b))
        assert np.allclose(cross, np.cross(a, b))
        assert np.isclose(triple[0], np.dot(np.cross(a, b), c))
        assert np.isclose(det[0], np.linalg.det(np.column_stack((a, b, c))))
        assert np.isclose(norm[0], np.linalg.norm(a))
        assert np.isclose(normc[0], norm[0])
        assert np.allclose(unit, a / norm[0])


def test_mathlib_vector_transformations():
    rng = np.random.default_rng(0)

    for _ in range(100):
        xyz = rng.standard_normal(3)
        rpz = np.array(
            [rng.uniform(0.1, 10.0), rng.uniform(-np.pi, np.pi), rng.standard_normal()]
        )
        vxyz = rng.standard_normal(3)
        vrpz = rng.standard_normal(3)

        phi = rpz[1]
        C = np.array(
            [[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]]
        )

        xyz_out = np.empty(3, dtype=np.float64)
        rpz_out = np.empty(3, dtype=np.float64)
        vxyz_out = np.empty(3, dtype=np.float64)
        vrpz_out = np.empty(3, dtype=np.float64)
        LIBASCOT.math_test_eval_vector_transformations(
            xyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            rpz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            vxyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            vrpz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            xyz_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            rpz_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            vxyz_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            vrpz_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )

        assert np.allclose(
            xyz_out,
            np.array(
                [
                    rpz[0] * np.cos(rpz[1]),
                    rpz[0] * np.sin(rpz[1]),
                    rpz[2],
                ]
            ),
        )
        assert np.allclose(
            rpz_out,
            np.array(
                [
                    np.hypot(xyz[0], xyz[1]),
                    np.arctan2(xyz[1], xyz[0]),
                    xyz[2],
                ]
            ),
        )
        assert np.allclose(vxyz_out, C @ vrpz)
        assert np.allclose(vrpz_out, C.T @ vxyz)


def test_mathlib_jac_rpz2xyz():
    np.random.seed(0)
    jacrpz = np.random.rand(12)
    jacxyz = np.zeros(12, dtype=np.float64)
    rpz0 = np.random.rand(3)
    rpz1 = rpz0 + 1e-3
    xyz1 = np.array([rpz0[0] * np.cos(rpz0[1]), rpz0[0] * np.sin(rpz0[1]), rpz0[2]])
    xyz2 = np.array([rpz1[0] * np.cos(rpz1[1]), rpz1[0] * np.sin(rpz1[1]), rpz1[2]])

    A0 = np.asarray(jacrpz[:3], dtype=float)
    J = np.asarray(jacrpz[3:], dtype=float).reshape(3, 3)

    delta = rpz0 - rpz1
    Arpz1 = A0 + J @ delta
    C = np.array(
        [
            [np.cos(rpz0[1]), np.sin(rpz0[1]), 0],
            [-np.sin(rpz0[1]), np.cos(rpz0[1]), 0],
            [0, 0, 1],
        ]
    )
    Axyz0 = A0 @ C
    Axyz1 = Arpz1 @ C

    LIBASCOT.math_jac_rpz2xyz(
        jacxyz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        jacrpz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        rpz0[0],
        rpz0[1],
    )

    A0 = np.asarray(jacxyz[:3], dtype=float)
    J = np.asarray(jacxyz[3:], dtype=float).reshape(3, 3)
    delta = xyz1 - xyz2
    A1 = A0 + J @ delta

    assert not np.allclose(A0, A1, atol=1e-3, rtol=0)
    assert np.allclose(A1, Axyz1, atol=1e-3, rtol=0)
    assert np.allclose(A0, Axyz0, atol=1e-8, rtol=0)


def test_mathlib_jac_xyz2rpz():

    np.random.seed(1)
    jacxyz0 = np.random.rand(12)
    jacxyz1 = np.zeros(12, dtype=np.float64)
    jacrpz = np.zeros(12, dtype=np.float64)
    xyz = np.random.rand(3)
    rpz = np.array(
        [np.sqrt(xyz[0] ** 2 + xyz[1] ** 2), np.arctan2(xyz[1], xyz[0]), xyz[2]]
    )

    LIBASCOT.math_jac_xyz2rpz(
        jacrpz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        jacxyz0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        rpz[0],
        rpz[1],
    )

    LIBASCOT.math_jac_rpz2xyz(
        jacxyz1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        jacrpz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        rpz[0],
        rpz[1],
    )

    assert np.allclose(jacxyz1, jacxyz0, atol=1e-8, rtol=0)


@pytest.mark.parametrize("positive_x", [0, 1], ids=("x<0", "x>0"))
@pytest.mark.parametrize("pi", [-2 * np.pi, 0.0, 2 * np.pi], ids=("-2pi", "0", "2pi"))
@pytest.mark.parametrize(
    "alpha1, alpha2, beta, k",
    [
        (np.pi / 4, 3 * np.pi / 4, 2 * np.pi / 4, 0.5),
        (3 * np.pi / 4, np.pi / 4, 2 * np.pi / 4, 0.5),
        (np.pi / 4, 2 * np.pi / 4, 3 * np.pi / 4, -1.0),
        (2 * np.pi / 4, np.pi / 4, 3 * np.pi / 4, -1.0),
        (7 * np.pi / 4, np.pi / 4, 0.0, 0.5),
        (np.pi / 4, 7 * np.pi / 4, 0.0, 0.5),
    ],
    ids=(
        "cross-positive",
        "cross-negative",
        "no-cross-positive",
        "no-cross-negative",
        "cross-at-zero-positive",
        "cross-at-zero-negative",
    ),
)
def test_math_crossed_plane(alpha1, alpha2, beta, k, pi, positive_x):
    if not positive_x:
        alpha1 += np.pi
        alpha2 += np.pi
        beta += np.pi
    kout = LIBASCOT.math_crossed_plane(alpha1 + pi, alpha2 + pi, beta)
    assert np.isclose(k, kout)
