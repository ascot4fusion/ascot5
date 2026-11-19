"""Test bfieldspline2d.py."""

import pytest
import numpy as np
import unyt

from a5py.utils.arrayop import assert_input_validation
from a5py.physlib import aeq
from a5py.templates import TransmuteAnalyticalBfieldToSplines, PremadeMagneticField


@pytest.fixture
def data(request, ascot):
    """Test data generated with templates that produce ITER-like fields."""
    optional = getattr(request, "param", False)

    template = PremadeMagneticField(ascot, field="iter-baseline")
    correct = template.create_input(preview=True)

    template = TransmuteAnalyticalBfieldToSplines(
        ascot,
        field=correct,
        nr=100,
        rlim=unyt.unyt_array([4, 8.2], "m"),
        nz=200,
        zlim=unyt.unyt_array([-4, 4], "m"),
    )
    out = template.create_input(preview=True).export()

    if not optional:
        del out["axisrz"]
        del out["psilimits"]
        del out["br"]
        del out["bz"]
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
    "data, optional",
    [
        pytest.param(False, False, id="required-data"),
        pytest.param(True, True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldspline2d_init(ascot, data, optional):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_bfieldspline2d(**data)
    out = obj.export()
    assert set(out.keys()) == {
        "rgrid",
        "zgrid",
        "axisrz",
        "psilimits",
        "psi",
        "br",
        "bz",
        "bphi",
    }

    if not optional:
        assert np.allclose(obj.axisrz.v, [6.2, 0.0])
        assert np.allclose(obj.psilimits.v, [0, 0])
        assert np.allclose(obj.br.v, 0.0)
        assert np.allclose(obj.bz.v, 0.0)

    for k, v in data.items():
        assert np.allclose(getattr(obj, k), v)

        with pytest.raises(AttributeError):
            setattr(obj, k, 0)

        assert np.allclose(out[k], v)
        assert not getattr(obj, k).flags.writeable
        assert out[k].flags.writeable


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldspline2d_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    assert_input_validation(
        func=ascot.data.create_bfieldspline2d,
        inputs={
            "rgrid": (data["rgrid"], (-1,), "m", None),
            "zgrid": (data["zgrid"], (-1,), "m", None),
            "axisrz": (data["axisrz"], (2,), "m", None),
            "psilimits": (data["psilimits"], (2,), "Wb/rad", None),
            "psi": (data["psi"], (100, 200), "Wb/rad", None),
            "br": (data["br"], (100, 200), "T", None),
            "bz": (data["bz"], (100, 200), "T", None),
            "bphi": (data["bphi"], (100, 200), "T", None),
        },
    )


def test_bfieldspline2d_eval_axis_rz(ascot, data, helpers):
    """Test that magnetic axis evaluationr returns correct values."""
    obj = ascot.data.create_bfieldspline2d(**data)
    r, z = helpers.bfield_eval_axis_rz(obj, np.array([0.0]))
    assert r[0] == obj.axisrz[0]
    assert z[0] == obj.axisrz[1]


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldspline2d_eval_psi_rho_and_gradients(ascot, data, helpers):
    """Test values and gradients of psi and rho."""
    obj = ascot.data.create_bfieldspline2d(**data)

    template = PremadeMagneticField(ascot, field="iter-baseline")
    correct = template.create_input(preview=True)

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n)) * unyt.m
    phi = np.zeros(n) * unyt.deg
    t = np.zeros(n) * unyt.s

    psi0, rho0 = helpers.bfield_eval_psi_dpsi_and_rho_drho(obj, r, phi, z, t)

    psi = (
        aeq.psi(
            r / correct.rmajor,
            z / correct.rmajor,
            correct.coefficients[:-1],
            correct.coefficients[-1],
        )
        * correct.psiscaling
    )
    rho = np.sqrt(
        (psi - correct.psilimits[0]) / (correct.psilimits[1] - correct.psilimits[0])
    )

    assert np.isclose(psi, psi0[::4], atol=1e-3).all()
    assert np.isclose(rho, rho0, atol=1e-3).all()

    dr = 1.0e-3 * unyt.m
    dz = 1.0e-3 * unyt.m
    dphi = 1.0e-4 * unyt.deg
    helpers.assert_gradient_of_psi(obj, r, dr, phi, dphi, z, dz)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldspline2d_value_of_b(ascot, data, helpers):
    """Test that the value of magnetic field is correct."""
    obj = ascot.data.create_bfieldspline2d(**data)

    template = PremadeMagneticField(ascot, field="iter-baseline")
    correct = template.create_input(preview=True)

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n)) * unyt.m
    phi = np.zeros(n) * unyt.deg
    t = np.zeros(n) * unyt.s

    psi, _ = helpers.bfield_eval_psi_dpsi_and_rho_drho(obj, r, phi, z, t)
    b, _ = helpers.bfield_eval_b_db(obj, r, phi, z, t)

    assert np.allclose(b[0, :], -((psi[3, :] * unyt.Wb / unyt.m) / r).to("T").v)
    assert np.allclose(b[2, :], ((psi[1, :] * unyt.Wb / unyt.m) / r).to("T").v)

    bphi = correct.rmajor * correct.bphi / r
    assert np.allclose(b[1, :], bphi)


def test_bfieldspline2d_eval_b_and_bdb(ascot, data, helpers):
    """Test evaluation of magnetic field derivatives."""
    obj = ascot.data.create_bfieldspline2d(**data)

    n = 100
    r = (4 + (8.2 - 4) * np.random.rand(n)) * unyt.m
    z = (-2 + (2 - -2) * np.random.rand(n)) * unyt.m
    phi = 360 * np.random.rand(n) * unyt.deg

    dr = 1.0e-3 * unyt.m
    dz = 1.0e-3 * unyt.m
    dphi = 1.0e-4 * unyt.deg
    helpers.assert_gradient_of_b(obj, r, dr, phi, dphi, z, dz)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldspline2d_eval_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_bfieldspline2d(**data)

    r, phi, z, t = (
        unyt.unyt_array([0.9], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )

    with pytest.raises(ValueError):
        helpers.bfield_eval_psi_dpsi_and_rho_drho(obj, r, phi.to("rad"), z, t)

    with pytest.raises(ValueError):
        helpers.bfield_eval_b_db(obj, r, phi.to("rad"), z, t)
