import ctypes
import numpy as np
import pytest

from a5py.data.atomic import *
from a5py.templates.import_atomic import ImportAtomicData, calculate_maxwellian_average
from a5py import Ascot

from a5py.libascot import LIBASCOT, init_fun

init_fun(
    "Atomic_eval_bms",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.POINTER(Struct),
    restype=ctypes.c_uint32,
)


init_fun(
    "Atomic_eval_cx",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.POINTER(Struct),
    restype=ctypes.c_uint32,
)


init_fun(
    "Atomic_eval_ii",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.POINTER(Struct),
    restype=ctypes.c_uint32,
)

init_fun(
    "Atomic_eval_ei",
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,
    ctypes.POINTER(Struct),
    restype=ctypes.c_uint32,
)


def test_init():
    ascot = Ascot()

    d = BMSdata(
        "H",
        0,
        "H",
        1,
        np.logspace(15, 21, 7),
        np.linspace(1e3, 3e4, 10),
        np.logspace(3, 6, 20),
        np.full((20, 10, 7), 1.0),
    )

    #leaf = ascot.data.create_atomic([d])
    leaf = ascot.data.create_atomic([], [])

    ratecoeff = np.zeros((1,))
    LIBASCOT.Atomic_eval_bms(
        ratecoeff.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1,
        1,
        (1e5*unyt.eV).to("J").v,
        ((1.0*unyt.amu).to("kg")).v,
        1,
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        (1e3*unyt.eV).to("J").v,
        1e19,
        np.full((1,), 1e19, dtype=np.float64).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        0,
        ctypes.byref(leaf._cdata),
    )
    print(ratecoeff)


def test_import_bms():
    ascot = Ascot()
    template = ImportAtomicData(ascot)
    template.import_openadas_bms(
        "bmsdata.txt", ("H", 0), ("H", 1)
    )
    leaf = template.create_input()

    ratecoeff = np.zeros((1,))
    LIBASCOT.Atomic_eval_bms(
        ratecoeff.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1,
        1,
        (1e5*unyt.eV).to("J").v,
        ((1.0*unyt.amu).to("kg")).v,
        1,
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        (1e3*unyt.eV).to("J").v,
        1e19,
        np.full((1,), 1e19, dtype=np.float64).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        0,
        ctypes.byref(leaf._cdata),
    )
    print(ratecoeff)


def test_import_cx():
    ascot = Ascot()
    template = ImportAtomicData(ascot)
    template.import_aladdin_cx(
        "cxdata.txt", ("H", 1), "H"
    )
    leaf = template.create_input()

    ratecoeff = np.zeros((1,))
    err = LIBASCOT.Atomic_eval_cx(
        ratecoeff.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1,
        (1.0e2*unyt.eV).to("J").v,
        ((1.0*unyt.amu).to("kg")).v,
        1,
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        (9.999e1*unyt.eV).to("J").v,
        np.full((1,), 1e19, dtype=np.float64).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        0,
        ctypes.byref(leaf._cdata),
    )
    print(err, ((ratecoeff/1e19)*unyt.m**3/unyt.s).to("cm**3/s").v)


def test_import_ii():
    ascot = Ascot()
    template = ImportAtomicData(ascot)
    #template.import_aladdin_ii(
    #   "iidata.txt", ("He4", 1), ("H", 1), "target"
    #)
    leaf = template.create_input()

    ratecoeff = np.zeros((1,))
    err = LIBASCOT.Atomic_eval_ii(
        ratecoeff.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        1,
        (80.0e3*unyt.eV*4).to("J").v,
        ((4.0*unyt.amu).to("kg")).v,
        1,
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        np.array([1], dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        (1e3*unyt.eV).to("J").v,
        np.full((1,), 1e19, dtype=np.float64).ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        0,
        ctypes.byref(leaf._cdata),
    )
    print(err, ((ratecoeff/1e19)*unyt.m**3/unyt.s).to("cm**3/s").v)
    # 6e-10 cm**3/s


def test_import_ei():
    ascot = Ascot()
    template = ImportAtomicData(ascot)
    #template.import_aladdin_ei(
    #    "eidata.txt", ("He4", 1),
    #)
    leaf = template.create_input()

    ratecoeff = np.zeros((1,))
    err = LIBASCOT.Atomic_eval_ei(
        ratecoeff.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        2,
        4,
        1,
        (1.0e2*unyt.eV).to("J").v,
        1e19,
        0,
        ctypes.byref(leaf._cdata),
    )
    print(err, ((ratecoeff/1e19)*unyt.m**3/unyt.s).to("cm**3/s").v)


def test_maxwellian_average():
    mass = (1.0 * unyt.amu).to("kg").v
    energy = (np.logspace(-10, 10, 100) * unyt.eV).to("J").v
    sigma = np.ones_like(energy)
    T = (1e2 * unyt.eV).to("J").v
    sigma_v = calculate_maxwellian_average(T*1000, T, mass, mass, energy, sigma)
    print(sigma_v)
