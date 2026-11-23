"""Test plasmalinear2d.py."""

import pytest
import numpy as np
import unyt

from a5py.physlib import Species
from a5py.utils.arrayop import assert_input_validation


def profile(a, b, c, r, z):
    """Evaluate a positive linear profile a + b*R + c*z."""
    return a + b * r + c * z


@pytest.fixture
def data(request):
    """Test data generated with positive linear profiles."""
    nspecies = getattr(request, "param", 1)
    species_names = ["D", "T"] if nspecies == 2 else ["H1"]

    rgrid = np.array([1.0, 1.5, 2.0]) * unyt.m
    zgrid = np.array([-1.0, -0.5, 0.0, 0.5, 1.0]) * unyt.m

    R, Z = np.meshgrid(
        rgrid.to_value(unyt.m),
        zgrid.to_value(unyt.m),
        indexing="ij",
    )

    # Electron profiles.
    ne = 1e20 * profile(1.0, 0.1, 0.05, R, Z) * unyt.m**-3
    Te = 1e3 * profile(2.0, 0.2, 0.1, R, Z) * unyt.eV

    # Ion profiles, with distinct profiles for each species.
    ni = (
        np.stack(
            [
                1e20 * profile(0.8 + 0.2 * i, 0.05 + 0.02 * i, 0.03 + 0.01 * i, R, Z)
                for i in range(len(species_names))
            ],
            axis=-1,
        )
        * unyt.m**-3
    )

    Ti = 1e3 * profile(1.0, 0.1, 0.05, R, Z) * unyt.eV

    rotation = 1e5 * profile(1.0, 0.1, 0.05, R, Z) * unyt.rad / unyt.s

    return {
        "rgrid": rgrid,
        "zgrid": zgrid,
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
def test_plasmalinear2d_init(ascot, data):
    """Test object initialization.

    Check that properties return correct values which are immutable, and that
    the export method returns correct values.
    """
    obj = ascot.data.create_plasmalinear2d(**data)
    out = obj.export()
    assert set(out.keys()) == {
        "rgrid",
        "zgrid",
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
def test_plasmalinear2d_init_optional(ascot, data):
    for optional in ["ne", "Te", "charge", "rotation"]:
        del data[optional]

    obj = ascot.data.create_plasmalinear2d(**data)
    assert np.allclose(obj.Te, data["Ti"])
    assert np.allclose(obj.charge, [1] * len(data["species"]))
    assert np.allclose(
        obj.rotation, np.full((data["rgrid"].size, data["zgrid"].size), 0)
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
def test_plasmalinear2d_validate_parameters(ascot, data):
    """Check that all input parameters are validated."""

    def wrap(*args, **kwargs):
        ascot.data.create_plasmalinear2d(*args, **kwargs, species=data["species"])

    assert_input_validation(
        func=wrap,
        inputs={
            "rgrid": (data["rgrid"], (-1,), "m", None),
            "zgrid": (data["zgrid"], (-1,), "m", None),
            "ne": (data["ne"], (3, 5), "m**(-3)", None),
            "Te": (data["Te"], (3, 5), "eV", None),
            "ni": (data["ni"], (3, 5, 2), "m**(-3)", None),
            "Ti": (data["Ti"], (3, 5), "eV", None),
            "charge": (data["charge"], (2,), "e", "i4"),
            "rotation": (data["rotation"], (3, 5), "rad/s", None),
        },
    )

    data["species"][0] = "Orava"
    with pytest.raises(ValueError, match="Unknown species"):
        ascot.data.create_plasmalinear2d(**data)

    data["species"] = []
    with pytest.raises(ValueError, match="At least one species must be specified."):
        ascot.data.create_plasmalinear2d(**data)


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(False, id="single-species"),
        pytest.param(True, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear2d_get_species(ascot, data, helpers):
    """Test that correct species number, and atomic and mass number are stored."""
    anum, znum, charge, mass = helpers.assert_species(
        ascot.data.create_plasmalinear2d(**data)
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
def test_plasmalinear2d_eval_nT(ascot, data, helpers):
    """Test evaluation of plasma density and temperature."""
    obj = ascot.data.create_plasmalinear2d(**data)

    nmrk = 100
    rho = np.zeros(nmrk)
    r = (1.0 + np.random.rand(nmrk)) * unyt.m
    z = (-1.0 + 2 * np.random.rand(nmrk)) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s
    ne, Te, ni, Ti = helpers.plasma_eval_nT(obj, rho, r, phi.to("rad"), z, t)

    expected_ne = 1e20 * profile(1.0, 0.1, 0.05, r.v, z.v) * unyt.m**-3
    expected_Te = 1e3 * profile(2.0, 0.2, 0.1, r.v, z.v) * unyt.eV
    expected_Ti = 1e3 * profile(1.0, 0.1, 0.05, r.v, z.v) * unyt.eV
    expected_ni = np.zeros((nmrk, len(data["species"]))) * unyt.m ** (-3)

    for i in range(len(data["species"])):
        expected_ni[:, i] = 1e20 * profile(
            0.8 + 0.2 * i, 0.05 + 0.02 * i, 0.03 + 0.01 * i, r.v, z.v
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
def test_plasmalinear2d_eval_nT_outside(ascot, data, helpers):
    """Test that evaluationg outside the data range raises an error."""
    obj = ascot.data.create_plasmalinear2d(**data)

    rho, r, phi, z, t = (
        unyt.unyt_array([0.0], "1"),
        unyt.unyt_array([0.1], "m"),
        unyt.unyt_array([0.0], "deg"),
        unyt.unyt_array([0.0], "m"),
        np.zeros(1) * unyt.s,
    )
    with pytest.raises(ValueError):
        helpers.plasma_eval_nT(obj, rho, r, phi.to("rad"), z, t)


@pytest.mark.parametrize(
    "r0, error",
    [
        pytest.param(1.0, False, id="inside"),
        pytest.param(2.0, True, id="outside"),
    ],
)
@pytest.mark.parametrize(
    "data",
    [
        pytest.param(2, id="multi-species"),
    ],
    indirect=["data"],
)
def test_plasmalinear2d_eval_rotation(ascot, data, helpers, r0, error):
    """Test evaluation of plasma rotation."""
    obj = ascot.data.create_plasmalinear2d(**data)

    nmrk = 10
    rho = np.zeros(nmrk)
    r = (r0 + np.random.rand(nmrk)) * unyt.m
    z = (-1.0 + 2 * np.random.rand(nmrk)) * unyt.m
    phi = np.zeros(nmrk) * unyt.deg
    t = np.zeros(nmrk) * unyt.s
    if error:
        with pytest.raises(ValueError):
            helpers.plasma_eval_rotation(obj, rho, r, phi.to("rad"), z, t)
    else:
        rotation = helpers.plasma_eval_rotation(obj, rho, r, phi.to("rad"), z, t)
        expected_rotation = 1e5 * profile(1.0, 0.1, 0.05, r.v, z.v) * unyt.rad / unyt.s
        assert np.allclose(rotation, expected_rotation)
