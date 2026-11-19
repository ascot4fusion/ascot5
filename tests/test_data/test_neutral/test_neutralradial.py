"""Test neutralradial.py."""

import pytest
import numpy as np
import unyt

from a5py.physlib import Species
from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data(request):
    """Test data generated with templates that produce ITER-like fields."""
    nspecies = getattr(request, "param", False)

    if nspecies == 2:
        out = {
            "rhogrid": np.linspace(0.0, 1.0, 20) * unyt.dimensionless,
            "species": ["H1", "He4"],
            "n0": np.ones((20, 2)) * unyt.m ** (-3),
            "T0": np.ones((20, 2)) * unyt.eV,
        }
    else:
        out = {
            "rhogrid": np.linspace(0.0, 1.0, 20) * unyt.dimensionless,
            "species": ["H1"],
            "n0": np.ones((20, 1)) * unyt.m ** (-3),
            "T0": np.ones((20, 1)) * unyt.eV,
        }
    for i in range(len(out["species"])):
        out["species"][i] = Species.from_string(out["species"][i])
        out["n0"][:, i] *= (i + 1) * 1e20 * np.linspace(0.1, 1.0, 20)
        out["T0"][:, i] *= (i + 1) * 1e3 * np.linspace(0.1, 1.0, 20)

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
    "data",
    [
        pytest.param(1, id="single-species"),
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_neutralradial_init(ascot, data):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_neutralradial(**data)
    out = obj.export()
    assert set(out.keys()) == {
        "rhogrid",
        "species",
        "n0",
        "T0",
    }

    for k, v in data.items():

        with pytest.raises(AttributeError):
            setattr(obj, k, 0)

        if k == "species":
            for i, s in enumerate(v):
                assert getattr(obj, k)[i].name == s.name
                assert out[k][i].name == s.name

            assert isinstance(getattr(obj, k), tuple)
            assert isinstance(out[k], list)
        else:
            assert np.allclose(getattr(obj, k), v)
            assert np.allclose(out[k], v)
            assert not getattr(obj, k).flags.writeable
            assert out[k].flags.writeable


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_neutralradial_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    def wrap(*args, **kwargs):
        ascot.data.create_neutralradial(*args, **kwargs, species=data["species"])

    assert_input_validation(
        func=wrap,
        inputs={
            "rhogrid": (data["rhogrid"], (-1,), "1", None),
            "n0": (data["n0"], (20, 2), "m**(-3)", None),
            "T0": (data["T0"], (20, 2), "eV", None),
        },
    )

    data["species"][0] = "Orava"
    with pytest.raises(ValueError, match="Unknown species"):
        ascot.data.create_neutralradial(**data)

    data["species"] = []
    with pytest.raises(ValueError, match="At least one species must be specified."):
        ascot.data.create_neutralradial(**data)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_neutralradial_get_species(ascot, data, helpers):
    """Test that correct species number, and atomic and mass number are stored."""
    anum, znum = helpers.assert_species(ascot.data.create_neutralradial(**data))
    assert anum.size == len(data["species"])
    assert znum.size == len(data["species"])

    for i in range(len(data["species"])):
        assert anum[i] == data["species"][i].anum
        assert znum[i] == data["species"][i].znum


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(1, id="single-species"),
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_neutralradial_eval_nT(ascot, data, helpers):
    """Test evaluation of neutral density and temperature."""
    obj = ascot.data.create_neutralradial(**data)

    nmrk = 100
    rho = np.random.rand(nmrk)
    r = np.zeros(nmrk) * unyt.m
    z = np.zeros(nmrk) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s
    n, T = helpers.neutral_eval_nT(obj, rho, r, phi.to("rad"), z, t)

    expected_density = np.zeros((nmrk, len(data["species"]))) * unyt.m ** (-3)
    expected_temperature = np.zeros((nmrk, len(data["species"]))) * unyt.eV

    for i in range(len(data["species"])):
        expected_density[:, i] = np.interp(rho, data["rhogrid"], data["n0"][:, i])
        expected_temperature[:, i] = np.interp(rho, data["rhogrid"], data["T0"][:, i])

    assert np.allclose(n, expected_density)
    assert np.allclose(T, expected_temperature)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_neutralradial_eval_nT_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_neutralradial(**data)

    rho, r, phi, z, t = (
        unyt.unyt_array([1.1], "1"),
        unyt.unyt_array([0.0], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )
    with pytest.raises(ValueError):
        helpers.neutral_eval_nT(obj, rho, r, phi.to("rad"), z, t)
