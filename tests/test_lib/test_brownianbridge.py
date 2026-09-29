import numpy as np
import pytest
import ctypes

from a5py.libascot import LIBASCOT, init_fun


class BrownianBridge(ctypes.Structure):
    _fields_ = [
        ("nmrk", ctypes.c_size_t),
        ("ndim", ctypes.c_size_t),
        ("ntime", ctypes.c_size_t),
        ("time", ctypes.POINTER(ctypes.c_double)),
        ("wiener", ctypes.POINTER(ctypes.c_double)),
    ]


init_fun(
    "BrownianBridge_init",
    ctypes.POINTER(BrownianBridge),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    restype=ctypes.c_int32,
)

init_fun(
    "BrownianBridge_free",
    ctypes.POINTER(BrownianBridge),
)

init_fun(
    "BrownianBridge_generate0th",
    ctypes.POINTER(BrownianBridge),
    ctypes.c_size_t,
    ctypes.c_double,
)

init_fun(
    "BrownianBridge_generate5",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(BrownianBridge),
    restype=ctypes.c_int32,
)

init_fun(
    "BrownianBridge_clear",
    ctypes.POINTER(BrownianBridge),
    ctypes.c_double,
    ctypes.c_size_t,
)


def test_generate0th():
    bbridge = BrownianBridge()
    LIBASCOT.BrownianBridge_init(ctypes.byref(bbridge), 2, 5, 3)
    LIBASCOT.BrownianBridge_generate0th(ctypes.byref(bbridge), 0, 0.5)
    LIBASCOT.BrownianBridge_generate0th(ctypes.byref(bbridge), 1, 1.5)

    assert bbridge.time[0] == 0.5
    assert bbridge.wiener[0] == 0
    assert bbridge.wiener[1] == 0
    assert bbridge.time[5 + 0] == 1.5
    assert bbridge.wiener[5 * 3 + 0] == 0
    assert bbridge.wiener[5 * 3 + 1] == 0
    for i in range(1, 5):
        assert bbridge.time[i] == -1
        assert bbridge.time[5 + i] == -1

    LIBASCOT.BrownianBridge_free(ctypes.byref(bbridge))


@pytest.mark.parametrize(
    "sample, expected_time, expected_w, expected_err",
    [
        (
            [
                (0, 0.5),
            ],
            (0.1, -1, -1, 0.5, 0.1, -1, -1, -1),
            (0., 0., 0., 0.63, 0., 0., 0., 0.),
            False,
        ),
        (
            [
                (0, 0.5),
                (0, 0.5),
            ],
            (0.1, -1, -1, 0.5, 0.1, -1, -1, -1),
            (0., 0., 0., 0.63, 0., 0., 0., 0.),
            False,
        ),
        (
            [
                (1, 0.5),
            ],
            (0.1, -1, -1, -1, 0.1, -1, -1, 0.5),
            (0., 0., 0., 0., 0., 0., 0., 0.63),
            False,
        ),
        (
            [(0, 0.5), (0, 0.4), (0, 0.6)],
            (0.1, 0.6, 0.4, 0.5, 0.1, -1, -1, -1),
            (0., 1.58, 1.02, 0.63, 0., 0., 0., 0.),
            False,
            ),
        (
            [(0, 0.5), (0, 0.4), (0, 0.6), (0, 0.3)],
            (0.1, 0.6, 0.4, 0.5, 0.1, -1, -1, -1),
            (0., 1.58, 1.02, 0.63, 0., 0., 0., 0.),
            True,
        ),
    ],
)
def test_generate5(sample, expected_time, expected_w, expected_err):
    bbridge = BrownianBridge()
    LIBASCOT.BrownianBridge_init(ctypes.byref(bbridge), 2, 4, 5)
    for i in range(4*5*2):
        bbridge.wiener[i] = 0
    LIBASCOT.BrownianBridge_generate0th(ctypes.byref(bbridge), 0, 0.1)
    LIBASCOT.BrownianBridge_generate0th(ctypes.byref(bbridge), 1, 0.1)


    for i, (mrk, t) in enumerate(sample):
        wiener = np.zeros(5) + np.nan
        wiener[0] = 1 + i
        err = LIBASCOT.BrownianBridge_generate5(
            wiener.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            t,
            mrk,
            ctypes.byref(bbridge),
        )
    assert np.allclose(bbridge.time[:8], expected_time)
    assert np.allclose(bbridge.wiener[:8*5:5], expected_w, atol=1e-2)
    assert err == expected_err

    #real var = (t - t0) * (1 + bridge * ((t1 - t) / (t1 - t0) - 1));
    #real temp = bridge * (t - t0) / (t1 - t0);
    #wiener[i] = W0[i] + temp * (W1[i] - W0[i]) + sqrt(var) * wiener[i];

    LIBASCOT.BrownianBridge_free(ctypes.byref(bbridge))


def test_clear():
    bbridge = BrownianBridge()
    LIBASCOT.BrownianBridge_init(ctypes.byref(bbridge), 1, 4, 5)
    LIBASCOT.BrownianBridge_generate0th(ctypes.byref(bbridge), 0, 0.1)
    for mrk, t in [(0, 0.5), (0, 0.4), (0, 0.6)]:
        wiener = np.zeros(5) + np.nan
        LIBASCOT.BrownianBridge_generate5(
            wiener.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            t,
            mrk,
            ctypes.byref(bbridge),
        )

    LIBASCOT.BrownianBridge_clear(ctypes.byref(bbridge), 0.4, 0)
    assert np.allclose(bbridge.time[:4], (0.4, 0.6, -1, 0.5))
