"""Test efieldpotential1d.py."""

import pytest

import unyt
import numpy as np

from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data():
    out = {
        "rhogrid": np.linspace(0, 1, 50),
        "potential": unyt.unyt_array(np.linspace(0, 1, 50), "V"),
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
def test_efieldpotential1d_init(
    ascot,
    data,
):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_efieldpotential1d(**data)
    out = obj.export()
    for k, v in data.items():
        assert np.allclose(getattr(obj, k), v)
        with pytest.raises(AttributeError):
            setattr(obj, k, 0)

        assert not getattr(obj, k).flags.writeable

        assert np.allclose(out[k], v)
        assert out[k].flags.writeable


def test_efieldpotential1d_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    assert_input_validation(
        func=ascot.data.create_efieldpotential1d,
        inputs={
            "rhogrid": (data["rhogrid"], (-1,), "1", None),
            "potential": (data["potential"], (50,), "V", None),
        },
    )


def test_efieldpotential1d_eval_e(ascot, data, bdata, helpers):
    """Test that the electric field evaluation returns correct values."""
    obj = ascot.data.create_efieldpotential1d(**data)
    bobj = ascot.data.create_bfieldspline2d(**bdata)

    n = 100
    axisrz = bobj.axisrz
    rminor = np.random.rand(n) * unyt.m
    pol = 2 * np.pi * np.random.rand(n)
    r, z = axisrz[0] + rminor * np.cos(pol), axisrz[1] + rminor * np.sin(pol)
    phi = 360 * np.random.rand(n) * unyt.deg
    t = np.zeros(n) * unyt.s

    e = helpers.efield_eval_e(obj, bobj, r, phi.to("rad"), z, t)

    assert np.allclose(e[0], -(r - axisrz[0]) / rminor, atol=1e-3)
    assert np.allclose(e[1], 0)
    assert np.allclose(e[2], -(z - axisrz[1]) / rminor, atol=1e-3)


def test_efieldpotential1d_eval_outside(ascot, data, bdata, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_efieldpotential1d(**data)
    bobj = ascot.data.create_bfieldspline2d(**bdata)

    r, phi, z, t = (
        unyt.unyt_array([0.9], "m"),
        unyt.unyt_array([0.], "deg"),
        unyt.unyt_array([0.], "m"),
        np.zeros(1) * unyt.s,
    )

    with pytest.raises(ValueError):
        helpers.efield_eval_e(obj, bobj, r, phi.to("rad"), z, t)
