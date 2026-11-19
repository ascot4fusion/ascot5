"""Test bfieldcartesian.py."""

import pytest

import unyt
import numpy as np

from a5py.utils.arrayop import assert_input_validation
from a5py.physlib.formulas import cart2pol, cart2pol_vec


@pytest.fixture
def data(request):
    """Test data."""
    optional = getattr(request, "param", False)
    out = {
        "bxyz": unyt.unyt_array([1.0, 2.0, 3.0], "T"),
        "jacobian": unyt.unyt_array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ],
            "T/m",
        ),
        "axisrz": unyt.unyt_array([1.0, 0.0], "m"),
        "rhoval": unyt.unyt_array(0.5, "1"),
    }
    if optional:
        out.update(
            {
                "psival": unyt.unyt_array(0.8, "Wb/rad"),
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
    "data, optional",
    [
        pytest.param(False, False, id="required-data"),
        pytest.param(True, True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldcartesian_init(ascot, data, optional):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_bfieldcartesian(**data)
    out = obj.export()

    assert set(out.keys()) == {"bxyz", "jacobian", "axisrz", "rhoval", "psival"}

    for k, v in data.items():
        assert np.allclose(getattr(obj, k), v)
        assert not getattr(obj, k).flags.writeable

        with pytest.raises(AttributeError):
            setattr(obj, k, 0)

        assert np.allclose(out[k], v)
        assert out[k].flags.writeable

    if not optional:
        assert np.allclose(obj.psival, 0.5 * unyt.Wb / unyt.rad)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(True, id="optional-data"),
    ],
    indirect=["data"],
)
def test_bfieldcartesian_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    assert_input_validation(
        func=ascot.data.create_bfieldcartesian,
        inputs={
            "bxyz": (data["bxyz"], (3,), "T", None),
            "jacobian": (data["jacobian"], (3,3), "T/m", None),
            "axisrz": (data["axisrz"], (2,), "m", None),
            "rhoval": (data["rhoval"], (), "1", None),
            "psival": (data["psival"], (), "Wb/rad", None),
        },
    )


def test_bfieldcartesian_eval_axis_rz(ascot, data, helpers):
    """Test evaluation of axis location."""
    obj = ascot.data.create_bfieldcartesian(**data)
    r, z = helpers.bfield_eval_axis_rz(obj, np.array([0.0]))
    assert r[0] == data["axisrz"][0]
    assert z[0] == data["axisrz"][1]


def test_bfieldcartesian_eval_psi_rho_and_gradients(ascot, data, helpers):
    """Test evaluation of psi, rho, and their gradients."""
    obj = ascot.data.create_bfieldcartesian(**data)

    psidpsi, rho = helpers.bfield_eval_psi_dpsi_and_rho_drho(
        obj, np.array([1.0]), np.array([0.0]), np.array([0.0]), np.array([0.0])
    )

    assert psidpsi[0] == obj.psival.v
    assert rho[0] == obj.rhoval.v
    for i in range(3):
        assert psidpsi[i + 1] == 0


def test_bfieldcartesian_value_of_b(ascot, data, helpers):
    """Test that magnetic field returns correct values."""
    obj = ascot.data.create_bfieldcartesian(**data)

    n = 100
    x, y, z, t = (
        unyt.unyt_array(1 - 2 * np.random.rand(n), "m"),
        unyt.unyt_array(1 - 2 * np.random.rand(n), "m"),
        unyt.unyt_array(1 - 2 * np.random.rand(n), "m"),
        np.zeros(n),
    )
    r, phi, z = cart2pol(x, y, z)

    b, _ = helpers.bfield_eval_b_db(
        obj,
        r,
        phi.to("rad").v,
        z,
        t,
    )
    br0, bphi0, bz0 = b[0], b[1], b[2]

    jac = obj.jacobian
    bx = obj.bxyz[0] + jac[0, 0] * x + jac[0, 1] * y + jac[0, 2] * z
    by = obj.bxyz[1] + jac[1, 0] * x + jac[1, 1] * y + jac[1, 2] * z
    bz = obj.bxyz[2] + jac[2, 0] * x + jac[2, 1] * y + jac[2, 2] * z

    br, bphi, bz = cart2pol_vec(x, bx, y, by, z, bz)
    assert np.allclose(br, br0)
    assert np.allclose(bz, bz0)
    assert np.allclose(bphi, bphi0)


def test_bfieldcartesian_eval_b_and_bdb(ascot, data, helpers):
    """Test evaluation of magnetic field derivatives."""
    obj = ascot.data.create_bfieldcartesian(**data)

    n = 100
    r = (1 + 2 * np.random.rand(n)) * unyt.m
    z = (-1 + 2 * np.random.rand(n)) * unyt.m
    phi = 360 * np.random.rand(n) * unyt.deg

    dr = 1.0e-1 * unyt.m
    dz = 1.0e-1 * unyt.m
    dphi = 1.0e-3 * unyt.deg
    helpers.assert_gradient_of_b(obj, r, dr, phi, dphi, z, dz)
