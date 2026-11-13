import ctypes as ct
import pytest
import unyt
import numpy as np

from a5py.engine.solve import SimData

from a5py.libascot import LIBASCOT, init_fun
from a5py.data.hist import Struct as HistStruct, Hist
from a5py import Ascot, SimulationOptions

real = ct.c_double


class MarkerGyroOrbit(ct.Structure):
    _fields_ = [
        ("size", ct.c_size_t),
        ("znum", ct.c_int),
        ("anum", ct.c_int),
        ("mass", real),
        ("charge", ct.POINTER(ct.c_int)),
        ("r", ct.POINTER(real)),
        ("phi", ct.POINTER(real)),
        ("z", ct.POINTER(real)),
        ("p_r", ct.POINTER(real)),
        ("p_phi", ct.POINTER(real)),
        ("p_z", ct.POINTER(real)),
        ("time", ct.POINTER(real)),
        ("br", ct.POINTER(real)),
        ("bphi", ct.POINTER(real)),
        ("bz", ct.POINTER(real)),
        ("dbrdr", ct.POINTER(real)),
        ("dbphidr", ct.POINTER(real)),
        ("dbzdr", ct.POINTER(real)),
        ("dbrdphi", ct.POINTER(real)),
        ("dbphidphi", ct.POINTER(real)),
        ("dbzdphi", ct.POINTER(real)),
        ("dbrdz", ct.POINTER(real)),
        ("dbphidz", ct.POINTER(real)),
        ("dbzdz", ct.POINTER(real)),
        ("bounces", ct.POINTER(ct.c_int)),
        ("weight", ct.POINTER(real)),
        ("cputime", ct.POINTER(real)),
        ("rho", ct.POINTER(real)),
        ("theta", ct.POINTER(real)),
        ("id", ct.POINTER(ct.c_size_t)),
        ("endcond", ct.POINTER(ct.c_uint)),
        ("walltile", ct.POINTER(ct.c_size_t)),
        ("mileage", ct.POINTER(real)),
        ("running", ct.POINTER(ct.c_int)),
        ("err", ct.POINTER(ct.c_uint)),
        ("index", ct.POINTER(ct.c_size_t)),
        ("rng_count", ct.POINTER(ct.c_size_t)),
    ]


init_fun(
    "MarkerGyroOrbit_allocate",
    ct.POINTER(MarkerGyroOrbit),
    ct.c_size_t,
    restype=ct.c_int,
)

init_fun(
    "MarkerGyroOrbit_deallocate",
    ct.POINTER(MarkerGyroOrbit),
)

# fusion_source_compute_product(
#   Reaction reaction, real Q, real vol, MarkerGyroOrbit *react1,
#   MarkerGyroOrbit *react2, MarkerGyroOrbit *prod1, MarkerGyroOrbit *prod2)

init_fun(
    "fusion_source_sample_thermal",
    ct.POINTER(SimData),
    ct.c_size_t,
    ct.POINTER(MarkerGyroOrbit),
    ct.POINTER(ct.c_double),
)

init_fun(
    "fusion_source_sample_beam",
    ct.POINTER(HistStruct),
    ct.c_size_t,
    ct.c_double,
    ct.POINTER(MarkerGyroOrbit),
    ct.POINTER(ct.c_double),
)

init_fun(
    "fusion_source_compute_product",
    ct.c_int,
    ct.c_double,
    ct.POINTER(MarkerGyroOrbit),
    ct.POINTER(MarkerGyroOrbit),
    ct.POINTER(MarkerGyroOrbit),
    ct.POINTER(MarkerGyroOrbit),
)

def test_fusion_source_compute_product():
    ascot = Ascot()
    ascot.data.create_bfieldcartesian(
        bxyz=np.array([0, 1, 1]), jacobian=np.zeros((3, 3)), axisrz=(1, 0), rhoval=0.5
    )
    ascot.data.create_plasmalinear1d(
        ["H"], np.linspace(0, 1, 5), 1e20 * np.ones((5, 1)), 1e3 * np.ones(5)
    )
    nsample = 100
    sim = SimData()
    sim.bfield.use(ascot.data.bfield.active)
    sim.plasma.use(ascot.data.plasma.active)
    react1, react2, prod1, prod2 = (MarkerGyroOrbit(), MarkerGyroOrbit(), MarkerGyroOrbit(), MarkerGyroOrbit())
    LIBASCOT.MarkerGyroOrbit_allocate(ct.byref(react1), nsample)
    LIBASCOT.MarkerGyroOrbit_allocate(ct.byref(react2), nsample)
    LIBASCOT.MarkerGyroOrbit_allocate(ct.byref(prod1), nsample)
    LIBASCOT.MarkerGyroOrbit_allocate(ct.byref(prod2), nsample)
    for j, mrk in enumerate([react1, react2, prod1, prod2]):
        for i in range(nsample):
            mrk.id[i] = 1 + i
            mrk.r[i] = 1.0
            mrk.z[i] = 0.0
            mrk.phi[i] = 0.0
            mrk.rho[i] = 0.5
            mrk.rng_count[i] = 0
            if j == 0:
                react1.mass = (2.0 * unyt.amu).to("kg").v
                react1.charge[i] = 1
            elif j == 1:
                react2.mass = (3.0 * unyt.amu).to("kg").v
                react2.charge[i] = 1
            elif j == 2:
                prod1.mass = (4.0 * unyt.amu).to("kg").v
                prod1.charge[i] = 2
            elif j == 3:
                prod2.mass = (1.0 * unyt.amu).to("kg").v
                prod2.charge[i] = 0

    density = np.zeros(1, dtype="f8")
    LIBASCOT.fusion_source_sample_thermal(
        ct.byref(sim),
        0,
        ct.byref(react1),
        density.ctypes.data_as(ct.POINTER(ct.c_double)),
    )
    LIBASCOT.fusion_source_sample_thermal(
        ct.byref(sim),
        0,
        ct.byref(react2),
        density.ctypes.data_as(ct.POINTER(ct.c_double)),
    )

    LIBASCOT.fusion_source_compute_product(
        0,
        1.0,
        ct.byref(react1),
        ct.byref(react2),
        ct.byref(prod1),
        ct.byref(prod2),
    )
    pr = np.array(react1.p_r[:nsample]) + np.array(react2.p_r[:nsample]) - np.array(prod1.p_r[:nsample]) - np.array(prod2.p_r[:nsample])
    pz = np.array(react1.p_z[:nsample]) + np.array(react2.p_z[:nsample]) - np.array(prod1.p_z[:nsample]) - np.array(prod2.p_z[:nsample])
    pphi = np.array(react1.p_phi[:nsample]) + np.array(react2.p_phi[:nsample]) - np.array(prod1.p_phi[:nsample]) - np.array(prod2.p_phi[:nsample])
    #print(np.amax(np.abs(pz)), np.amax(np.abs(pphi)), np.amax(np.abs(pr)), np.amax(np.array(prod1.p_phi[:nsample])))
    ekinr1 = 0.5*(np.array(react1.p_r[:nsample])**2 + np.array(react1.p_z[:nsample])**2 + np.array(react1.p_phi[:nsample])**2)/react1.mass
    ekinr2 = 0.5*(np.array(react2.p_r[:nsample])**2 + np.array(react2.p_z[:nsample])**2 + np.array(react2.p_phi[:nsample])**2)/react2.mass
    ekinp1 = 0.5*(np.array(prod1.p_r[:nsample])**2 + np.array(prod1.p_z[:nsample])**2 + np.array(prod1.p_phi[:nsample])**2)/prod1.mass
    ekinp2 = 0.5*(np.array(prod2.p_r[:nsample])**2 + np.array(prod2.p_z[:nsample])**2 + np.array(prod2.p_phi[:nsample])**2)/prod2.mass
    assert np.amax(np.abs(pz)) < 1e-34
    assert np.amax(np.abs(pphi)) < 1e-34
    assert np.amax(np.abs(pr)) < 1e-34

    Q = (17.6e6*unyt.eV).to("J").v
    #print(ekinr1+ekinr2+Q- ekinp1-ekinp2, Q)
    assert np.amax(np.abs(ekinr1+ekinr2+Q- ekinp1-ekinp2)) < 1e-20

    LIBASCOT.MarkerGyroOrbit_deallocate(ct.byref(react1))
    LIBASCOT.MarkerGyroOrbit_deallocate(ct.byref(react2))
    LIBASCOT.MarkerGyroOrbit_deallocate(ct.byref(prod1))
    LIBASCOT.MarkerGyroOrbit_deallocate(ct.byref(prod2))


def test_fusion_source_sample_beam():
    opt = SimulationOptions.from_dict(
        histograms=[
            {
                "dimensions": [
                    ("r", 0., 1., 2),
                    ("z", -1., 1., 2),
                    ("pr", -1.e-21, 1.e-21, 3),
                    ("pphi", -1.e-21, 1.e-21, 3),
                    ("pz", -1.e-21, 1.e-21, 3),
                ],
                "charge_interval": (1, 1),
            }
        ]
    )
    hist = Hist.from_params(opt.histograms[0])
    hist._cdata.values[13] = 1.0
    hist._cdata.values[14] = 1.0
    print(hist._cdata.values_ref.shape)

    ascot = Ascot()
    ascot.data.create_bfieldcartesian(
        bxyz=np.array([0, 1, 1]), jacobian=np.zeros((3, 3)), axisrz=(1, 0), rhoval=0.5
    )
    ascot.data.create_plasmalinear1d(
        ["H"], np.linspace(0, 1, 5), 1e20 * np.ones((5, 1)), 1e3 * np.ones(5)
    )
    nsample = 100000
    sim = SimData()
    sim.bfield.use(ascot.data.bfield.active)
    sim.plasma.use(ascot.data.plasma.active)
    mrk = MarkerGyroOrbit()
    LIBASCOT.MarkerGyroOrbit_allocate(ct.byref(mrk), nsample)
    for i in range(nsample):
        mrk.id[i] = 1 + i
        mrk.r[i] = 1.0
        mrk.z[i] = 0.0
        mrk.phi[i] = 0.0
        mrk.rho[i] = 0.5
        mrk.rng_count[i] = 0
        mrk.mass = (1.0 * unyt.amu).to("kg").v
        mrk.charge[i] = 1

    density = np.zeros(1, dtype="f8")
    LIBASCOT.fusion_source_sample_beam(
        ct.byref(hist._cdata),
        0,
        1.0,
        ct.byref(mrk),
        density.ctypes.data_as(ct.POINTER(ct.c_double)),
    )
    import matplotlib.pyplot as plt

    plt.hist(mrk.p_z[:nsample], np.linspace(-1,1,10)*1e-21)
    plt.show()
    LIBASCOT.MarkerGyroOrbit_deallocate(ct.byref(mrk))


def test_fusion_source_sample_thermal():
    ascot = Ascot()
    ascot.data.create_bfieldcartesian(
        bxyz=np.array([0, 1, 1]), jacobian=np.zeros((3, 3)), axisrz=(1, 0), rhoval=0.5
    )
    ascot.data.create_plasmalinear1d(
        ["H"], np.linspace(0, 1, 5), 1e20 * np.ones((5, 1)), 1e3 * np.ones(5)
    )
    nsample = 100000
    sim = SimData()
    sim.bfield.use(ascot.data.bfield.active)
    sim.plasma.use(ascot.data.plasma.active)
    mrk = MarkerGyroOrbit()
    LIBASCOT.MarkerGyroOrbit_allocate(ct.byref(mrk), nsample)
    for i in range(nsample):
        mrk.id[i] = 1 + i
        mrk.r[i] = 1.0
        mrk.z[i] = 0.0
        mrk.phi[i] = 0.0
        mrk.rho[i] = 0.5
        mrk.rng_count[i] = 0
        mrk.mass = (1.0 * unyt.amu).to("kg").v
        mrk.charge[i] = 1

    density = np.zeros(1, dtype="f8")
    LIBASCOT.fusion_source_sample_thermal(
        ct.byref(sim),
        0,
        ct.byref(mrk),
        density.ctypes.data_as(ct.POINTER(ct.c_double)),
    )
    import matplotlib.pyplot as plt

    plt.hist(mrk.p_r[:nsample], 20)
    plt.show()
    LIBASCOT.MarkerGyroOrbit_deallocate(ct.byref(mrk))
