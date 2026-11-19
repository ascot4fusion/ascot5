"""Test efieldpotential2d.py."""

import pytest

import unyt
import numpy as np

from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data():
    """Test data.

    Potential is equal to the distance from the axis which makes it easy to
    calculate the gradient analytically for comparison.
    """
    rgrid = np.linspace(1.0, 3.0, 80) * unyt.m
    zgrid = np.linspace(-1.0, 1.0, 100) * unyt.m
    axisrz = np.array([2.0, 0.0]) * unyt.m

    R, Z = np.meshgrid(rgrid, zgrid, indexing="ij")
    potential = np.sqrt((R - axisrz[0]) ** 2 + (Z - axisrz[1]) ** 2).v * unyt.V
    out = {
        "rgrid": rgrid,
        "zgrid": zgrid,
        "potential": potential,
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
def test_efieldpotential2d_init(
    ascot,
    data,
):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_efieldpotential2d(**data)
    out = obj.export()
    for k, v in data.items():
        assert np.allclose(getattr(obj, k), v)
        with pytest.raises(AttributeError):
            setattr(obj, k, 0)

        assert not getattr(obj, k).flags.writeable

        assert np.allclose(out[k], v)
        assert out[k].flags.writeable


def test_efieldpotential2d_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    assert_input_validation(
        func=ascot.data.create_efieldpotential2d,
        inputs={
            "rgrid": (data["rgrid"], (-1,), "m", None),
            "zgrid": (data["zgrid"], (-1,), "m", None),
            "potential": (data["potential"], (80, 100), "V", None),
        },
    )


def test_efieldpotential2d_eval_e(ascot, data, bdata, helpers):
    """Test that the electric field evaluation returns correct values."""
    obj = ascot.data.create_efieldpotential2d(**data)
    bobj = ascot.data.create_bfieldspline2d(**bdata)

    n = 100
    r, phi, z, t = (
        unyt.unyt_array(1.0 + 2 * np.random.rand(n), "m"),
        unyt.unyt_array(360 * np.random.rand(n), "deg"),
        unyt.unyt_array(-1.0 + 2 * np.random.rand(n), "m"),
        np.zeros(n) * unyt.s,
    )
    e = helpers.efield_eval_e(obj, bobj, r, phi.to("rad"), z, t)

    axisrz = bobj.axisrz
    rminor = np.sqrt((r - axisrz[0]) ** 2 + (z - axisrz[1]) ** 2)

    assert np.allclose(e[0, :], -((r - axisrz[0]) / rminor).v, atol=1e-3)
    assert np.allclose(e[1, :], 0.0)
    assert np.allclose(e[2, :], -((z - axisrz[1]) / rminor).v, atol=1e-3)


def test_efieldpotential2d_eval_outside(ascot, data, bdata, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_efieldpotential2d(**data)
    bobj = ascot.data.create_bfieldspline2d(**bdata)

    r, phi, z, t = (
        unyt.unyt_array([0.9], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )

    with pytest.raises(ValueError):
        helpers.efield_eval_e(obj, bobj, r, phi.to("rad"), z, t)
