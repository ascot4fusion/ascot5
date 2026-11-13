import ctypes

import unyt
import numpy as np
import pytest

from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "boschhale_sigma",
    ctypes.c_int,
    ctypes.c_double,
    restype=ctypes.c_double,
)

init_fun(
    "boschhale_sigmav",
    ctypes.c_int,
    ctypes.c_double,
    restype=ctypes.c_double,
)

@pytest.mark.parametrize(
    "reaction",
    [0, 1, 2, 3],
    ids=("DT_He4n", "DHe3_He4p", "DD_Tp", "DD_He3n"),
)
def test_cross_section(reaction):

    # This is from the Bosch-Hale paper. The firs column is energy in keV
    # and the rest are cross sections of different reactions in mbarn.
    cross_section_data = np.array([
        [3,   9.808e-3, 1.119e-11, 2.513e-4, 2.445e-4],
        [4,   1.073e-1, 1.718e-9,  2.146e-3, 2.093e-3],
        [5,   5.383e-1, 5.199e-8,  9.038e-3, 8.834e-3],
        [6,   1.749e0,  6.336e-7,  2.569e-2, 2.517e-2],
        [7,   4.335e0,  4.373e-6,  5.720e-2, 5.616e-2],
        [8,   8.968e0,  2.058e-5,  1.081e-1, 1.064e-1],
        [9,   1.632e1,  7.374e-5,  1.820e-1, 1.794e-1],
        [10,  2.702e1,  2.160e-4,  2.812e-1, 2.779e-1],
        [12,  6.065e1,  1.206e-3,  5.607e-1, 5.563e-1],
        [15,  1.479e2,  7.944e-3,  1.180e0,  1.178e0],
        [20,  4.077e2,  6.568e-2,  2.670e0,  2.691e0],
        [30,  1.381e3,  7.704e-1,  6.681e0,  6.857e0],
        [40,  2.817e3,  3.266e0,  1.116e1,  1.165e1],
        [50,  4.219e3,  8.688e0,  1.557e1,  1.649e1],
        [60,  4.984e3,  1.788e1,  1.971e1,  2.115e1],
        [70,  4.987e3,  3.144e1,  2.351e1,  2.554e1],
        [80,  4.545e3,  4.978e1,  2.698e1,  2.964e1],
        [100, 3.427e3,  1.021e2,  3.304e1,  3.701e1],
        [120, 2.576e3,  1.768e2,  3.812e1,  4.343e1],
        [140, 2.004e3,  2.740e2,  4.245e1,  4.904e1],
        [150, 1.792e3,  3.303e2,  4.438e1,  5.160e1],
        [200, 1.138e3,  6.378e2,  5.234e1,  6.239e1],
        [250, 8.135e2,  8.122e2,  5.831e1,  7.071e1],
        [300, 6.234e2,  7.731e2,  6.302e1,  7.732e1],
        [400, 4.126e2,  5.304e2,  7.005e1,  8.702e1],
    ])

    for ekin, sigma in zip(
        (cross_section_data[:, 0] * 1e3 * unyt.eV).to("J"),
        cross_section_data[:, 1 + reaction] * 1e-31 * unyt.m**2,
    ):
        sigma0 = LIBASCOT.boschhale_sigma(reaction, ekin) * unyt.m**2
        assert np.isclose(sigma0, sigma, atol=0., rtol=1e-3)


@pytest.mark.parametrize(
    "reaction",
    [0, 1, 2, 3],
    ids=("DT_He4n", "DHe3_He4p", "DD_Tp", "DD_He3n"),
)
def test_reactivity(reaction):

    # This is from the Bosch-Hale paper. The firs column is temperature in keV
    # and the rest are reactivities of different reactions in mbarn. Note that
    # the original Bosch & Hale contained errors and even errata probably has
    # an error at T = 1.3 keV, DD_He3n.
    reactivity_data = np.array(
        [
            [0.2, 1.254e-26, 1.414e-35, 4.640e-28, 4.482e-28],
            [0.3, 7.292e-25, 1.033e-32, 2.071e-26, 2.004e-26],
            [0.4, 9.344e-24, 6.537e-31, 2.237e-25, 2.168e-25],
            [0.5, 5.697e-23, 1.241e-29, 1.204e-24, 1.169e-24],
            [0.6, 2.253e-22, 1.166e-28, 4.321e-24, 4.200e-24],
            [0.7, 6.740e-22, 6.960e-28, 1.193e-23, 1.162e-23],
            [0.8, 1.662e-21, 3.032e-27, 2.751e-23, 2.681e-23],
            [1.0, 6.857e-21, 3.057e-26, 1.017e-22, 9.933e-23],
            [1.3, 3.174e-20, 3.708e-25, 4.143e-22, 4.062e-22],
            [1.5, 6.923e-20, 1.317e-24, 8.431e-22, 8.284e-22],
            [1.8, 1.773e-19, 6.053e-24, 1.976e-21, 1.948e-21],
            [2.0, 2.977e-19, 1.399e-23, 3.150e-21, 3.110e-21],
            [2.5, 8.425e-19, 7.477e-23, 7.969e-21, 7.905e-21],
            [3.0, 1.867e-18, 2.676e-22, 1.608e-20, 1.602e-20],
            [4.0, 5.974e-18, 1.710e-21, 4.428e-20, 4.447e-20],
            [5.0, 1.366e-17, 6.377e-21, 9.024e-20, 9.128e-20],
            [6.0, 2.554e-17, 1.739e-20, 1.545e-19, 1.573e-19],
            [8.0, 6.222e-17, 7.504e-20, 3.354e-19, 3.457e-19],
            [10.0, 1.136e-16, 2.126e-19, 5.781e-19, 6.023e-19],
            [12.0, 1.747e-16, 4.715e-19, 8.723e-19, 9.175e-19],
            [15.0, 2.740e-16, 1.175e-18, 1.390e-18, 1.481e-18],
            [20.0, 4.330e-16, 3.482e-18, 2.399e-18, 2.603e-18],
            [30.0, 6.681e-16, 1.363e-17, 4.728e-18, 5.271e-18],
            [40.0, 7.998e-16, 3.160e-17, 7.249e-18, 8.235e-18],
            [50.0, 8.649e-16, 5.554e-17, 9.838e-18, 1.133e-17],
        ]
    )

    for Ti, sigmav in zip(
        (reactivity_data[:, 0] * 1e3 * unyt.eV).to("J"),
        reactivity_data[:, 1 + reaction] * 1e-6 * unyt.m**3 / unyt.s,
    ):
        sigmav0 = LIBASCOT.boschhale_sigmav(reaction, Ti) * unyt.m**3 / unyt.s
        assert np.isclose(sigmav0, sigmav, atol=0., rtol=1e-3)
