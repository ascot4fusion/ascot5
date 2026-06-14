import pytest
import unyt
import numpy as np
import os

from a5py import Ascot
from a5py.data.marker.state import MarkerState
from a5py.templates import PremadeMagneticField

from a5py.data.orbit import Orbit
from a5py.data.options import SimulationOptions
from a5py.data.access.hdf5io import TreeFileManager

NMRK = 5
NPOINT = 10
NPOINT_I = np.array([NPOINT, NPOINT, NPOINT-1, 1, 0])


def generate_data(ids, quantity):
    """Generate fixed data for a given quantity and marker."""
    i = ids - 1
    r0, z0 = 6.2, 0.0
    radii = np.linspace(1., 2., NMRK)
    theta = np.linspace(0, 2*np.pi, NPOINT)
    match quantity:
        case "id":
            return np.full(NPOINT, i+1)
        case "mileage":
            return np.linspace(0., 1., NPOINT)
        case "r":
            return radii[i] * np.cos(theta) + r0
        case "z":
            return radii[i] * np.cos(theta) + z0
        case _:
            return np.full(NPOINT, 0)
        
@pytest.fixture()
def file():
    datamanager = TreeFileManager(
        root="root",
        filename="testascot.h5",
        file_exists=False,
        input_categories=[],
        )
    datamanager.set_note(name="orbit", node="root", note="") # Creates node
    yield datamanager.access_data(name="orbit", node="root")

    os.unlink("testascot.h5")

@pytest.fixture
def orbit():
    """Create an orbit object.

    For test cases the orbits are cocentric circles with varying radii. Some
    have all points, some have less, and some have none. The data is stored in
    random order.
    """
    opt = SimulationOptions.from_dict(
        orbit={"buffer_size": NPOINT, "interval": 0.01, "collect": "interval"},
        )
    orbit = Orbit.from_params(NMRK, opt.orbit)
    for i in np.arange(NMRK):
        idx = slice(NPOINT*i, NPOINT*(i+1))
        theta = np.linspace(0, 2*np.pi, NPOINT)
        orbit._cdata.r_ref[idx] = generate_data(i+1, "r")
        orbit._cdata.z_ref[idx] = generate_data(i+1, "z")

        orbit._cdata.mileage_ref[idx] = generate_data(i+1, "mileage")
        orbit._cdata.id_ref[NPOINT*i:NPOINT*i+NPOINT_I[i]] = 1 + i

    idx = np.random.permutation(NMRK*NPOINT)
    for field in ["r", "z", "phi", "p1", "p2", "p3", "mileage",
                  "id", "charge", "poincare", "simmode"]:
        arr = getattr(orbit._cdata, field + "_ref")
        arr[:] = arr[idx]
    yield orbit


def test_init_from_params():
    opt = SimulationOptions.from_dict(
        orbit={"buffer_size": NPOINT, "interval": 0.01, "collect": "interval"},
        )
    Orbit.from_params(NMRK, opt.orbit)


def test_all_returned_data_have_correct_number_of_points(orbit):
    ids = orbit.extract_fields("id")["id"]
    for i in range(NMRK):
        assert np.sum(ids == i + 1) == NPOINT_I[i]

    assert np.sum(ids == 0) == 0
    assert ids.size == np.sum(NPOINT_I)


def test_returned_data_is_sorted_by_id_and_stamp(orbit):
    ids = orbit.extract_fields("id")["id"]
    assert np.all(np.diff(ids) >= 0)

    stamp = orbit.extract_fields("mileage")["mileage"]
    for i in range(NMRK):
        assert np.all(np.diff(stamp[ids == i + 1]) >= 0)


@pytest.mark.parametrize(
    "fields",
    [
        ["id"],
        ["id", "mileage"],
    ],
    ids=["single-query", "multiple-queries"]
)
def test_return_all_queried_fields(orbit, fields):
    out = orbit.extract_fields(*fields)
    assert len(out) == len(fields)
    for field in fields:
        assert field in out

@pytest.mark.parametrize(
    "ids, expected",
    [
        (None, np.arange(1, NMRK)), ([1, 2], [1, 2]), ([NMRK+1], []),
        ([2, NMRK+1], [2]),
    ],
    ids=[
        "no-filter", "filter", "filter-nonexistent-id",
        "filter-existent-and-nonexistent-id",
        ],
)
def test_filter_by_id(orbit, ids, expected):
    out = orbit.extract_fields("id", filter=ids)
    assert np.isclose(np.unique(out["id"]), expected).all()


def test_all_quantities_have_correct_values(orbit):
    out = orbit.extract_fields(
        "r", "z", "phi", "p1", "p2", "p3", "mileage", "id",
        "charge", "poincare", "simmode", filter=1,
        )
    for qnt, val in out.items():
        assert np.allclose(val, generate_data(1, qnt))


def test_store_and_read_hdf5(orbit, file):
    r = orbit.r
    orbit.save(file)
    orbit2 = Orbit()
    orbit2._file = file
    assert all(orbit2.r == r)


def test_returned_arrays_are_copies(orbit):
    some_value = 0.999
    out = orbit.extract_fields("r")
    out["r"][:] = some_value
    assert np.all(out["r"] == some_value)

    out = orbit.extract_fields("r")
    assert not np.any(out["r"] == some_value)


def test_evaluation_gc(orbit):
    a5 = Ascot()
    template = PremadeMagneticField(a5, field="iter-baseline")
    template.create_input()
    mrk = a5.data.create_guidingcentermarker("alpha", 6.2*unyt.m, 0.*unyt.m, 10.*unyt.eV, 0.9, preview=True)#, ids=[1, 2, 3, 4, 5])
    state = MarkerState.from_markerinput(mrk, np.full((mrk.n, 12), 1.0))

    quantities = ["r", "z", "phi", "bnorm", "mileage", "id", "charge", "time", "ppar", "mu", "zeta", "ptor",]
    for qnt in quantities:
        orbit.evaluate_quantity("guidingcenter", state, a5.data.bfield.active, qnt, filter=[1])
    assert False
