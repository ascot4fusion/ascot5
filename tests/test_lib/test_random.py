import ctypes
import numpy as np
from scipy import stats
import pytest

from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "random_seed",
    ctypes.c_size_t,
)

init_fun(
    "random_uniform",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_size_t),
    ctypes.POINTER(ctypes.c_size_t),
    ctypes.c_size_t,
)

init_fun(
    "random_norm",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_size_t),
    ctypes.POINTER(ctypes.c_size_t),
    ctypes.c_size_t,
)

def generate_rnd(n, marker_id, rng_counter, kind):
    rnd = np.zeros(n, dtype=np.float64)
    marker_id = np.array([marker_id], dtype=np.int64)
    rng_counter = np.array([rng_counter], dtype=np.int64)
    if kind == "uniform":
        LIBASCOT.random_uniform(
            rnd.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
            marker_id.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t)),
            rng_counter.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t)),
            0,
        )
    elif kind == "norm":
        LIBASCOT.random_norm(
            rnd.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
            marker_id.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t)),
            rng_counter.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t)),
            0,
        )
    return rnd, rng_counter


def test_random_seed():
    LIBASCOT.random_seed(0)
    rnd1, _ = generate_rnd(1, 1, 0, "norm")
    LIBASCOT.random_seed(1)
    rnd2, _ = generate_rnd(1, 1, 0, "norm")
    assert rnd1[0] != rnd2[0]

def test_random_markerid():
    LIBASCOT.random_seed(0)
    rnd1, _ = generate_rnd(1, 1, 0, "norm")
    rnd2, _ = generate_rnd(1, 2, 0, "norm")
    rnd3, _ = generate_rnd(1, 1, 0, "norm")
    assert rnd1[0] != rnd2[0]
    assert rnd1[0] == rnd3[0]

def test_random_counter():
    LIBASCOT.random_seed(0)
    rnd1, _ = generate_rnd(1, 1, 0, "norm")
    rnd2, _ = generate_rnd(1, 1, 1, "norm")
    rnd3, _ = generate_rnd(1, 1, 0, "norm")
    assert rnd1[0] != rnd2[0]
    assert rnd1[0] == rnd3[0]

@pytest.mark.parametrize(
    "nrng, expected_counter", [(1, 1), (4, 1), (5, 2),]
)
def test_random_output(nrng, expected_counter):
    LIBASCOT.random_seed(0)

    rnd, rng_counter = generate_rnd(nrng, 1, 0, "uniform")
    assert (rnd > 0).all()
    assert rng_counter[0] == expected_counter

    rng_counter[0] = 0
    rnd[:] = 0
    rnd, rng_counter = generate_rnd(nrng, 1, 0, "norm")
    assert (rnd != 0).all()
    assert rng_counter[0] == expected_counter

def test_random_uniform():
    LIBASCOT.random_seed(0)
    rnd, _ = generate_rnd(100000, 1, 0, "norm")
    _, p_value = stats.kstest(rnd, "norm")
    assert p_value > 0.01

def test_random_normal():
    LIBASCOT.random_seed(0)
    rnd, _ = generate_rnd(100000, 1, 0, "uniform")
    _, p_value = stats.kstest(rnd, "uniform")
    assert p_value > 0.01
