import ctypes
import numpy as np
from scipy import stats

from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "random_test_init",
    ctypes.c_void_p,
    ctypes.c_size_t,
)

init_fun(
    "random_test_uniform_normal",
    ctypes.c_void_p,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
)

init_fun(
    "random_test_uniform_normal_simd",
    ctypes.c_void_p,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
)


def test_random_init():
    n = 2
    temp = np.empty(n, dtype=np.float64)

    LIBASCOT.random_test_init(None, 0)
    uniform1 = np.empty(n, dtype=np.float64)
    LIBASCOT.random_test_uniform_normal(
        None,
        n,
        uniform1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    uniform2 = np.empty(n, dtype=np.float64)
    LIBASCOT.random_test_uniform_normal(
        None,
        n,
        uniform2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    LIBASCOT.random_test_init(None, 0)
    uniform3 = np.empty(n, dtype=np.float64)
    LIBASCOT.random_test_uniform_normal(
        None,
        n,
        uniform3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    assert not np.allclose(uniform1, uniform2)
    assert np.allclose(uniform1, uniform3)


def test_random_scalar():
    n = 10000
    uniform = np.empty(n, dtype=np.float64)
    normal = np.empty(n, dtype=np.float64)
    LIBASCOT.random_test_init(None, 0)
    LIBASCOT.random_test_uniform_normal(
        None,
        n,
        uniform.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        normal.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    _, p_value = stats.kstest(normal, "norm")
    assert p_value > 0.01

    _, p_value = stats.kstest(uniform, "uniform")
    assert p_value > 0.01


def test_random_vector():

    n = 10000
    uniform = np.empty(n, dtype=np.float64)
    normal = np.empty(n, dtype=np.float64)
    LIBASCOT.random_test_init(None, 0)
    LIBASCOT.random_test_uniform_normal_simd(
        None,
        n,
        uniform.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        normal.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    _, p_value = stats.kstest(normal, "norm")
    assert p_value > 0.01

    _, p_value = stats.kstest(uniform, "uniform")
    assert p_value > 0.01
