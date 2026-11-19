"""Test bfieldanalytical.py."""

import pytest
import numpy as np
import unyt

from a5py.utils.arrayop import assert_input_validation
from a5py.physlib import aeq


@pytest.fixture
def data(request):
    """Test data.

    Data is for an ITER-like equilibrium. Including ripple is optional.
    """
    is3d = getattr(request, "param", False)

    out = {
        "rmajor": unyt.unyt_array(6.2, "m"),
        "axisrz": unyt.unyt_array([6.51738394, 0.0], "m"),
        "bphi": unyt.unyt_array(5.3, "T"),
        "psiscaling": unyt.unyt_array(200.0, "Wb/rad"),
        "coefficients": np.array(
            [
                0.08649128,
                0.3236476,
                0.52270472,
                -0.23197358,
                0.38073753,
                -0.35733467,
                -0.01487402,
                0.14801494,
                0.74018674,
                -0.43977189,
                -0.10713086,
                0.01278622,
                -0.155,
            ]
        ),
        "psilimits": unyt.unyt_array([-7.17652448, 0.0], "Wb/rad"),
    }
    if is3d:
        out.update(
            {
                "nripple": 16,
                "rminor": unyt.unyt_array(1.8, "m"),
                "rippledamping": unyt.unyt_array(0.02, "1"),
                "ripplescaling": unyt.unyt_array(1.1, "1"),
            }
        )
    return out


@pytest.mark.parametrize(
    "ascot",
    [
        pytest.param(True, id="with-file"),
        pytest.param(False, id="no-file"),
    ],
    indirect=True,
)
@pytest.mark.parametrize(
    "data, is3d",
    [
        pytest.param(False, False, id="2d"),
        pytest.param(True, True, id="3d"),
    ],
    indirect=["data"],
)
def test_bfieldanalytical_init(ascot, data, is3d):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_bfieldanalytical(**data)
    out = obj.export()

    assert set(out.keys()) == {
        "rmajor",
        "axisrz",
        "bphi",
        "psiscaling",
        "coefficients",
        "psilimits",
        "nripple",
        "rminor",
        "rippledamping",
        "ripplescaling",
    }

    if not is3d:
        assert obj.nripple == 0
        assert np.allclose(obj.rminor, unyt.unyt_array(1.0, "m"))
        assert np.allclose(obj.rippledamping, unyt.unyt_array(0.01, "1"))
        assert np.allclose(obj.ripplescaling, unyt.unyt_array(1.0, "1"))

    for k, v in data.items():
        assert np.allclose(getattr(obj, k), v)

        with pytest.raises(AttributeError):
            setattr(obj, k, 0)

        assert np.allclose(out[k], v)
        if k not in ["nripple"]:
            assert not getattr(obj, k).flags.writeable
            assert out[k].flags.writeable


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(True, id="3d"),
    ],
    indirect=["data"],
)
def test_bfieldanalytical_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    assert_input_validation(
        func=ascot.data.create_bfieldanalytical,
        inputs={
            "rmajor": (data["rmajor"], (), "m", None),
            "axisrz": (data["axisrz"], (2,), "m", None),
            "bphi": (data["bphi"], (), "T", None),
            "psiscaling": (data["psiscaling"], (), "Wb/rad", None),
            "coefficients": (data["coefficients"], (13,), None, None),
            "psilimits": (data["psilimits"], (2,), "Wb/rad", None),
            "nripple": (data["nripple"], (), None, "i4"),
            "rminor": (data["rminor"], (), "m", None),
            "rippledamping": (data["rippledamping"], (), "1", None),
            "ripplescaling": (data["ripplescaling"], (), "1", None),
        },
    )


def test_bfieldanalytical_eval_axis_rz(ascot, data, helpers):
    """Test that magnetic axis evaluationr returns correct values."""
    obj = ascot.data.create_bfieldanalytical(**data)
    r, z = helpers.bfield_eval_axis_rz(obj, np.array([0.0]))
    assert r[0] == data["axisrz"][0]
    assert z[0] == data["axisrz"][1]


def test_bfieldanalytical_eval_psi_rho_and_gradients(ascot, data, helpers):
    """Test values and gradients of psi and rho."""
    obj = ascot.data.create_bfieldanalytical(**data)

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n)) * unyt.m
    phi = 360 * np.random.rand(n) * unyt.deg
    t = np.zeros(n) * unyt.s

    psi0, rho0 = helpers.bfield_eval_psi_dpsi_and_rho_drho(obj, r, phi, z, t)

    psi = (
        aeq.psi(
            r / obj.rmajor, z / obj.rmajor, obj.coefficients[:-1], obj.coefficients[-1]
        )
        * obj.psiscaling
    )
    rho = np.sqrt((psi - obj.psilimits[0]) / (obj.psilimits[1] - obj.psilimits[0]))

    assert np.isclose(psi, psi0[::4]).all()
    assert np.isclose(rho, rho0).all()

    dr = 1.0e-3 * unyt.m
    dz = 1.0e-3 * unyt.m
    dphi = 1.0e-4 * unyt.deg
    helpers.assert_gradient_of_psi(obj, r, dr, phi, dphi, z, dz)


@pytest.mark.parametrize(
    "data, is3d",
    [
        pytest.param(False, False, id="2d"),
        pytest.param(True, True, id="3d"),
    ],
    indirect=["data"],
)
def test_bfieldanalytical_value_of_b(ascot, data, helpers, is3d):
    """Test that the value of magnetic field is correct."""
    obj = ascot.data.create_bfieldanalytical(**data)

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n)) * unyt.m
    phi = np.zeros(n) * unyt.deg
    t = np.zeros(n) * unyt.s

    psi, _ = helpers.bfield_eval_psi_dpsi_and_rho_drho(obj, r, phi, z, t)
    b, _ = helpers.bfield_eval_b_db(obj, r, phi, z, t)

    assert np.allclose(b[0, :], -((psi[3, :] * unyt.Wb / unyt.m) / r).to("T").v)
    assert np.allclose(b[2, :], ((psi[1, :] * unyt.Wb / unyt.m) / r).to("T").v)

    bphi = data["rmajor"] * data["bphi"] / r
    if is3d:
        rminor = np.sqrt((r - data["rmajor"]) ** 2 + (z - data["axisrz"][1]) ** 2)
        theta = np.arctan2(z - data["axisrz"][1], r - data["rmajor"])
        delta = (
            data["ripplescaling"]
            * np.exp(-0.5 * theta**2)
            * (rminor / data["rminor"]) ** data["rippledamping"]
        )
        bphi *= 1 + delta * np.cos(data["nripple"] * phi)
    assert np.allclose(b[1, :], bphi)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="2d"),
        pytest.param(True, id="3d"),
    ],
    indirect=["data"],
)
def test_bfieldanalytical_eval_gradient_of_b(ascot, data, helpers):
    """Test evaluation of magnetic field derivatives."""
    obj = ascot.data.create_bfieldanalytical(**data)

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n)) * unyt.m
    phi = 360 * np.random.rand(n) * unyt.deg

    dr = 1.0e-3 * unyt.m
    dz = 1.0e-3 * unyt.m
    dphi = 1.0e-4 * unyt.deg
    helpers.assert_gradient_of_b(obj, r, dr, phi, dphi, z, dz)
