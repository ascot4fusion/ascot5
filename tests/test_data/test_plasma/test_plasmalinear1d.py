"""Test plasmalinear1d.py."""

import pytest
import numpy as np
import unyt

from a5py.physlib import Species
from a5py.utils.arrayop import assert_input_validation


@pytest.fixture
def data(request):
    """Test data."""
    nspecies = getattr(request, "param", False)

    if nspecies == 2:
        out = {
            "rhogrid": np.array([0.0, 0.2, 0.5, 1.0]) * unyt.dimensionless,
            "species": ["D", "T"],
            "ne": 1e20 * np.linspace(0.1, 1.0, 4) * unyt.m ** (-3),
            "Te": 2e3 * np.linspace(0.1, 1.0, 4) * unyt.eV,
            "ni": np.ones((4, 2)) * unyt.m ** (-3),
            "Ti": 1e3 * np.linspace(0.1, 1.0, 4) * unyt.eV,
            "rotation": 1e5 * np.linspace(0.1, 1.0, 4) * unyt.rad / unyt.s,
            "charge": np.array([1, 1]) * unyt.e,
        }
    else:
        out = {
            "rhogrid": np.array([0.0, 0.2, 0.5, 1.0]) * unyt.dimensionless,
            "species": ["H1"],
            "ne": 1e20 * np.linspace(0.1, 1.0, 4) * unyt.m ** (-3),
            "Te": 2e3 * np.linspace(0.1, 1.0, 4) * unyt.eV,
            "ni": np.ones((4, 1)) * unyt.m ** (-3),
            "Ti": 1e3 * np.linspace(0.1, 1.0, 4) * unyt.eV,
            "rotation": 1e5 * np.linspace(0.1, 1.0, 4) * unyt.rad / unyt.s,
            "charge": np.array(
                [
                    1,
                ]
            )
            * unyt.e,
        }
    for i in range(len(out["species"])):
        out["species"][i] = Species.from_string(out["species"][i])
        out["ni"][:, i] *= 1e20 * np.linspace(0.1, 1.0, 4) * (i + 1) / 3.0

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
def test_plasmalinear1d_init(ascot, data):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_plasmalinear1d(**data)
    out = obj.export()
    assert set(out.keys()) == {
        "rhogrid",
        "species",
        "ne",
        "Te",
        "ni",
        "Ti",
        "rotation",
        "charge",
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
        pytest.param(1, id="single-species"),
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear1d_init_optional(ascot, data):
    for optional in ["ne", "Te", "charge", "rotation"]:
        del data[optional]

    obj = ascot.data.create_plasmalinear1d(**data)
    assert np.allclose(obj.Te, data["Ti"])
    assert np.allclose(obj.charge, [1] * len(data["species"]))
    assert np.allclose(obj.rotation, [0.0] * data["rhogrid"].size)

    ne = np.sum(obj.charge.v * data["ni"], axis=1)
    assert np.allclose(obj.ne, ne)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear1d_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    def wrap(*args, **kwargs):
        ascot.data.create_plasmalinear1d(*args, **kwargs, species=data["species"])

    assert_input_validation(
        func=wrap,
        inputs={
            "rhogrid": (data["rhogrid"], (-1,), "1", None),
            "ne": (data["ne"], (4,), "m**(-3)", None),
            "Te": (data["Te"], (4,), "eV", None),
            "ni": (data["ni"], (4, 2), "m**(-3)", None),
            "Ti": (data["Ti"], (4,), "eV", None),
            "charge": (data["charge"], (2,), "e", "i4"),
            "rotation": (data["rotation"], (4,), "rad/s", None),
        },
    )

    data["species"][0] = "Orava"
    with pytest.raises(ValueError, match="Unknown species"):
        ascot.data.create_plasmalinear1d(**data)

    data["species"] = []
    with pytest.raises(ValueError, match="At least one species must be specified."):
        ascot.data.create_plasmalinear1d(**data)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear1d_get_species(ascot, data, helpers):
    """Test that correct species number, and atomic and mass number are stored."""
    anum, znum, charge, mass = helpers.assert_species(
        ascot.data.create_plasmalinear1d(**data)
    )
    assert anum.size == len(data["species"])
    assert znum.size == len(data["species"])
    assert mass.size == len(data["species"]) + 1
    assert charge.size == len(data["species"]) + 1

    electron = Species.from_string("e")
    assert np.isclose(mass[0], electron.mass.to("kg"))
    assert charge[0] == -1 * unyt.e
    for i in range(len(data["species"])):
        assert anum[i] == data["species"][i].anum
        assert znum[i] == data["species"][i].znum
        assert np.isclose(mass[i + 1], data["species"][i].mass.to("kg"))
        assert charge[i + 1] == data["species"][i].znum


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(1, id="single-species"),
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear1d_eval_nT(ascot, data, helpers):
    """Test evaluation of plasma density and temperature."""
    obj = ascot.data.create_plasmalinear1d(**data)

    nmrk = 100
    rho = np.random.rand(nmrk)
    r = np.zeros(nmrk) * unyt.m
    z = np.zeros(nmrk) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s
    ne, Te, ni, Ti = helpers.plasma_eval_nT(obj, rho, r, phi.to("rad"), z, t)

    expected_ne = np.interp(rho, data["rhogrid"], data["ne"])
    expected_Te = np.interp(rho, data["rhogrid"], data["Te"])
    expected_Ti = np.interp(rho, data["rhogrid"], data["Ti"])
    expected_ni = np.zeros((nmrk, len(data["species"]))) * unyt.m ** (-3)

    for i in range(len(data["species"])):
        expected_ni[:, i] = np.interp(rho, data["rhogrid"], data["ni"][:, i])

    assert np.allclose(ne, expected_ne)
    assert np.allclose(Te, expected_Te)
    assert np.allclose(ni, expected_ni)
    assert np.allclose(Ti, expected_Ti)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear1d_eval_nT_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_plasmalinear1d(**data)

    rho, r, phi, z, t = (
        unyt.unyt_array([1.1], "1"),
        unyt.unyt_array([0.0], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )
    with pytest.raises(ValueError):
        helpers.plasma_eval_nT(obj, rho, r, phi.to("rad"), z, t)


@pytest.mark.parametrize(
    "rho0, error",
    [
        pytest.param(0.0, False, id="inside"),
        pytest.param(1.0, True, id="outside"),
    ],
)
@pytest.mark.parametrize(
    "data",
    [
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear2d_eval_rotation(ascot, data, helpers, rho0, error):
    """Test evaluation of plasma rotation."""
    obj = ascot.data.create_plasmalinear1d(**data)

    nmrk = 10
    rho = rho0 + np.random.rand(nmrk)
    r = np.ones(nmrk) * unyt.m
    z = np.zeros(nmrk) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s
    if error:
        with pytest.raises(ValueError):
            helpers.plasma_eval_rotation(obj, rho, r, phi.to("rad"), z, t)
    else:
        rotation = helpers.plasma_eval_rotation(obj, rho, r, phi.to("rad"), z, t)
        expected_rotation = (
            1e5
            * np.interp(rho, data["rhogrid"], np.linspace(0.1, 1.0, 4))
            * unyt.rad
            / unyt.s
        )
        assert np.allclose(rotation, expected_rotation)
