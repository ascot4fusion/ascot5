"""Test plasmadynamic1d.py."""

import pytest
import numpy as np
import unyt

from a5py.physlib import Species
from a5py.utils.arrayop import assert_input_validation


def profile(a, b, c, r, t):
    """Evaluate a positive linear profile a + b*R + c*t."""
    return a + b * r + c * t


@pytest.fixture
def data(request):
    """Test data generated with positive linear profiles."""
    nspecies = getattr(request, "param", 1)
    species_names = ["D", "T"] if nspecies == 2 else ["H1"]

    rhogrid = np.array([0.0, 0.2, 1.0]) * unyt.dimensionless
    timegrid = np.array([0.1, 0.5, 0.7, 1.1]) * unyt.s

    R, T = np.meshgrid(
        rhogrid.v,
        timegrid.to_value(unyt.s),
        indexing="ij",
    )

    # Electron profiles.
    ne = 1e20 * profile(1.0, 0.1, 0.05, R, T) * unyt.m**-3
    Te = 1e3 * profile(2.0, 0.2, 0.1, R, T) * unyt.eV

    # Ion profiles, with distinct profiles for each species.
    ni = (
        np.stack(
            [
                1e20 * profile(0.8 + 0.2 * i, 0.05 + 0.02 * i, 0.03 + 0.01 * i, R, T)
                for i in range(len(species_names))
            ],
            axis=-1,
        )
        * unyt.m**-3
    )

    Ti = 1e3 * profile(1.0, 0.1, 0.05, R, T) * unyt.eV

    rotation = 1e5 * profile(1.0, 0.1, 0.05, R, T) * unyt.rad / unyt.s

    return {
        "rhogrid": rhogrid,
        "timegrid": timegrid,
        "species": [Species.from_string(name) for name in species_names],
        "ne": ne,
        "Te": Te,
        "ni": ni,
        "Ti": Ti,
        "rotation": rotation,
        "charge": np.ones(len(species_names)) * unyt.e,
    }


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
def test_plasmadynamic1d_init(ascot, data):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_plasmadynamic1d(**data)
    out = obj.export()
    assert set(out.keys()) == {
        "rhogrid",
        "timegrid",
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
def test_plasmadynamic1d_init_optional(ascot, data):
    for optional in ["ne", "Te", "charge", "rotation"]:
        del data[optional]

    obj = ascot.data.create_plasmadynamic1d(**data)
    assert np.allclose(obj.Te, data["Ti"])
    assert np.allclose(obj.charge, [1] * len(data["species"]))
    assert np.allclose(
        obj.rotation, np.full((data["rhogrid"].size, data["timegrid"].size), 0)
    )

    ne = np.sum(obj.charge.v * data["ni"], axis=2)
    assert np.allclose(obj.ne, ne)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmadynamic1d_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    def wrap(*args, **kwargs):
        ascot.data.create_plasmadynamic1d(*args, **kwargs, species=data["species"])

    assert_input_validation(
        func=wrap,
        inputs={
            "rhogrid": (data["rhogrid"], (-1,), "1", None),
            "timegrid": (data["timegrid"], (-1,), "s", None),
            "ne": (data["ne"], (3, 4), "m**(-3)", None),
            "Te": (data["Te"], (3, 4), "eV", None),
            "ni": (data["ni"], (3, 4, 2), "m**(-3)", None),
            "Ti": (data["Ti"], (3, 4), "eV", None),
            "charge": (data["charge"], (2,), "e", "i4"),
            "rotation": (data["rotation"], (3, 4), "rad/s", None),
        },
    )

    data["species"][0] = "Orava"
    with pytest.raises(ValueError, match="Unknown species"):
        ascot.data.create_plasmadynamic1d(**data)

    data["species"] = []
    with pytest.raises(ValueError, match="At least one species must be specified."):
        ascot.data.create_plasmadynamic1d(**data)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmadynamic1d_get_species(ascot, data, helpers):
    """Test that correct species number, and atomic and mass number are stored."""
    anum, znum, charge, mass = helpers.assert_species(
        ascot.data.create_plasmadynamic1d(**data)
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
def test_plasmadynamic1d_eval_nT(ascot, data, helpers):
    """Test evaluation of plasma density and temperature."""
    obj = ascot.data.create_plasmadynamic1d(**data)

    nmrk = 100
    rho = np.random.rand(nmrk) * unyt.dimensionless
    r = np.zeros(nmrk) * unyt.m
    z = np.zeros(nmrk) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = (0.1 + np.random.rand(nmrk)) * unyt.s
    ne, Te, ni, Ti = helpers.plasma_eval_nT(obj, rho, r, phi.to("rad"), z, t)

    expected_ne = 1e20 * profile(1.0, 0.1, 0.05, rho.v, t.v) * unyt.m**-3
    expected_Te = 1e3 * profile(2.0, 0.2, 0.1, rho.v, t.v) * unyt.eV
    expected_Ti = 1e3 * profile(1.0, 0.1, 0.05, rho.v, t.v) * unyt.eV
    expected_ni = np.zeros((nmrk, len(data["species"]))) * unyt.m ** (-3)

    for i in range(len(data["species"])):
        expected_ni[:, i] = 1e20 * profile(
            0.8 + 0.2 * i, 0.05 + 0.02 * i, 0.03 + 0.01 * i, rho.v, t.v
        )

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
def test_plasmadynamic1d_eval_nT_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_plasmadynamic1d(**data)

    rho, r, phi, z, t = (
        unyt.unyt_array([2.0], "1"),
        unyt.unyt_array([1.0], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )
    with pytest.raises(ValueError):
        helpers.plasma_eval_nT(obj, rho, r, phi.to("rad"), z, t)


@pytest.mark.parametrize(
    "rho0, time0, error",
    [
        pytest.param(0.0, 0.1, False, id="inside"),
        pytest.param(2.0, 0.2, True, id="outside_rho"),
        pytest.param(0.1, 2.0, True, id="outside_time"),
    ],
)
@pytest.mark.parametrize(
    "data",
    [
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmadynamic1d_eval_rotation(ascot, data, helpers, rho0, time0, error):
    """Test evaluation of plasma rotation."""
    obj = ascot.data.create_plasmadynamic1d(**data)

    nmrk = 10
    rho = (rho0 + np.random.rand(nmrk)) * unyt.dimensionless
    r = np.ones(nmrk) * unyt.m
    z = np.zeros(nmrk) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = (time0 + np.random.rand(nmrk)) * unyt.s
    if error:
        with pytest.raises(ValueError):
            helpers.plasma_eval_rotation(obj, rho, r, phi.to("rad"), z, t)
    else:
        rotation = helpers.plasma_eval_rotation(obj, rho, r, phi.to("rad"), z, t)
        expected_rotation = 1e5 * profile(1.0, 0.1, 0.05, rho.v, t.v) * unyt.rad / unyt.s
        assert np.allclose(rotation, expected_rotation)
