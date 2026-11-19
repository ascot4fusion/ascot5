"""Test efieldcartesian.py."""

import pytest

import unyt
import numpy as np

from a5py.physlib.formulas import cart2pol, cart2pol_vec
from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data():
    """Test data."""
    out = {
        "exyz": unyt.unyt_array([1.0, 2.0, 3.0], "V/m"),
    }
    return out


@pytest.mark.parametrize(
    "ascot",
    [
        pytest.param(True, id="with-file"),
        pytest.param(False, id="no-file"),
    ],
    indirect=True,
)
def test_efieldcartesian_init(
    ascot,
    data,
):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_efieldcartesian(**data)

    assert np.allclose(obj.exyz, data["exyz"])
    with pytest.raises(AttributeError):
        setattr(obj, "exyz", 0)

    assert not getattr(obj, "exyz").flags.writeable

    out = obj.export()
    for k, v in data.items():
        assert np.allclose(out[k], v)
        assert out[k].flags.writeable


def test_efieldcartesian_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    assert_input_validation(
        func=ascot.data.create_efieldcartesian,
        inputs={"exyz": (data["exyz"], (3,), "V/m", None)},
    )


def test_efieldcartesian_eval_e(ascot, data, bdata, helpers):
    """Test that the electric field evaluation returns correct values.

    In this case, the vector is always constant in the cartesian basis.
    """
    obj = ascot.data.create_efieldcartesian(**data)
    bobj = ascot.data.create_bfieldspline2d(**bdata)

    n = 10
    x, y, z, t = (
        unyt.unyt_array(1.0 + np.random.rand(n), "m"),
        unyt.unyt_array(np.random.rand(n), "m"),
        unyt.unyt_array(np.random.rand(n), "m"),
        np.zeros(n) * unyt.s,
    )
    r, phi, z = cart2pol(x, y, z)
    e = helpers.efield_eval_e(obj, bobj, r, phi.to("rad"), z, t)

    er, ephi, ez = cart2pol_vec(
        x, data["exyz"][0], y, data["exyz"][1], z, data["exyz"][2]
    )
    erpz = unyt.unyt_array([er, ephi, ez], "V/m")
    assert np.allclose(e, erpz)
