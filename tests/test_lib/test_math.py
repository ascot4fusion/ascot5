import ctypes
import pytest

import numpy as np

from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "test_matrix_multiplication",
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t * 3,
)

def test_matrix_multiplication():
    LIBASCOT.test_matrix_multiplication(
        np.ones(9).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.ones(9).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.ones(9).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        np.array([3, 3, 3], dtype=np.int64),
    )

@pytest.mark.parametrize("positive_x", [0, 1], ids=("x<0", "x>0"))
@pytest.mark.parametrize(
        "pi",
        [-2*np.pi, 0., 2*np.pi],
        ids=("-2pi", "0", "2pi")
        )
@pytest.mark.parametrize(
        "alpha1, alpha2, beta, k",
        [
            (np.pi/4, 3*np.pi/4, 2*np.pi/4, 0.5),
            (3*np.pi/4, np.pi/4, 2*np.pi/4, 0.5),
            (np.pi/4, 2*np.pi/4, 3*np.pi/4, -1.),
            (2*np.pi/4, np.pi/4, 3*np.pi/4, -1.),
            (7*np.pi/4, np.pi/4, 0.0, 0.5),
            (np.pi/4, 7*np.pi/4, 0.0, 0.5),
        ],
        ids=("cross-positive", "cross-negative", "no-cross-positive",
             "no-cross-negative", "cross-at-zero-positive", "cross-at-zero-negative")
        )
def test_math_crossed_plane(alpha1, alpha2, beta, k, pi, positive_x):
    if not positive_x:
        alpha1 += np.pi
        alpha2 += np.pi
        beta += np.pi
    kout = LIBASCOT.math_crossed_plane(alpha1 + pi, alpha2 + pi, beta)
    assert np.isclose(k, kout)
