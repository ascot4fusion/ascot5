"""Test neutralarbitrary.py."""

import pytest
import numpy as np
import unyt

from a5py.physlib import Species
from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data(request):
    """Test data generated with templates that produce ITER-like fields."""
    nspecies = getattr(request, "param", False)
    nr, nphi, nz = 10, 21, 15
    out = {
        "rgrid": np.linspace(1.0, 2.0, nr) * unyt.m,
        "phigrid": np.linspace(0.0, 20.0, nphi) * unyt.deg,
        "zgrid": np.linspace(-1.0, 1.0, nz) * unyt.m,
    }

    if nspecies == 2:
        out.update({"species": ["H1", "He4"]})
    else:
        out.update({"species": ["H1"]})

    n0 = np.ones((nr, nphi - 1, nz, len(out["species"]))) * 1e20 * unyt.m ** (-3)
    T0 = np.ones((nr, nphi - 1, nz, len(out["species"]))) * 1e3 * unyt.eV
    for i in range(len(out["species"])):
        out["species"][i] = Species.from_string(out["species"][i])

        x = out["rgrid"].copy().v[:, None, None]
        y = out["phigrid"].copy().v[None, :, None]
        z = out["zgrid"].copy().v[None, None, :]
        y_periodic = np.where(y <= 10, y, 20 - y)
        data = 14 + 2 * x + 3 * y_periodic + 5 * z

        n0[:, :, :, i] *= data[:, :-1, :] * (i + 1)
        T0[:, :, :, i] *= data[:, :-1, :] * (i + 1)

    out.update({"n0": n0, "T0": T0})

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
def test_neutralarbitrary_init(ascot, data):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_neutralarbitrary(**data)
    out = obj.export()
    assert set(out.keys()) == {
        "rgrid",
        "phigrid",
        "zgrid",
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
def test_neutralarbitrary_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    def wrap(*args, **kwargs):
        ascot.data.create_neutralarbitrary(*args, **kwargs, species=data["species"])

    assert_input_validation(
        func=wrap,
        inputs={
            "rgrid": (data["rgrid"], (-1,), "m", None),
            "phigrid": (data["phigrid"], (-1,), "deg", None),
            "zgrid": (data["zgrid"], (-1,), "m", None),
            "n0": (data["n0"], (10, 20, 15, 2), "m**(-3)", None),
            "T0": (data["T0"], (10, 20, 15, 2), "eV", None),
        },
    )

    data["species"][0] = "Orava"
    with pytest.raises(ValueError, match="Unknown species"):
        ascot.data.create_neutralarbitrary(**data)

    data["species"] = []
    with pytest.raises(ValueError, match="At least one species must be specified."):
        ascot.data.create_neutralarbitrary(**data)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_neutralarbitrary_get_species(ascot, data, helpers):
    """Test that correct species number, and atomic and mass number are stored."""
    anum, znum = helpers.assert_species(ascot.data.create_neutralarbitrary(**data))
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
def test_neutralarbitrary_eval_nT(ascot, data, helpers):
    """Test evaluation of neutral density and temperature."""
    obj = ascot.data.create_neutralarbitrary(**data)

    nmrk = 100
    rho = np.zeros(nmrk)
    r = (1.0 + np.random.rand(nmrk)) * unyt.m
    z = (-1 + 2 * np.random.rand(nmrk)) * unyt.m
    phi = 20 * np.random.rand(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s
    n, T = helpers.neutral_eval_nT(obj, rho, r, phi.to("rad"), z, t)

    expected_density = np.zeros((nmrk, len(data["species"]))) * unyt.m ** (-3)
    expected_temperature = np.zeros((nmrk, len(data["species"]))) * unyt.eV

    y_periodic = np.where(phi.v <= 10, phi.v, 20 - phi.v)
    val = 14 + 2 * r.v + 3 * y_periodic + 5 * z.v

    for i in range(len(data["species"])):
        expected_density[:, i] = val * 1e20 * (i + 1)
        expected_temperature[:, i] = val * 1e3 * (i + 1)

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
def test_neutralarbitrary_eval_nT_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_neutralarbitrary(**data)

    rho, r, phi, z, t = (
        unyt.unyt_array([0.0], "1"),
        unyt.unyt_array([0.9], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )
    with pytest.raises(ValueError):
        helpers.neutral_eval_nT(obj, rho, r, phi.to("rad"), z, t)


def test_neutralarbitrary_periodicity(ascot, data, helpers):
    """Test that the field is periodic when interpolated."""
    obj = ascot.data.create_neutralarbitrary(**data)

    nmrk = 100
    rho = np.zeros(nmrk)
    r = (1.0 + np.random.rand(nmrk)) * unyt.m
    z = (-1 + 2 * np.random.rand(nmrk)) * unyt.m
    phi = 20 * np.random.rand(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s

    n0, T0 = helpers.neutral_eval_nT(obj, rho, r, phi.to("rad").v, z, t)

    phi += 20 * np.random.randint(low=1, high=18, size=nmrk) * unyt.deg
    n1, T1 = helpers.neutral_eval_nT(obj, rho, r, phi.to("rad").v, z, t)

    assert np.allclose(n0, n1)
    assert np.allclose(T0, T1)
