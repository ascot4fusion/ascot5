import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.engine.interpolate import evaluate


def test_neutralradial_init():
    ascot = Ascot()

    rhogrid = np.linspace(0.1, 1.0, 10)
    species = ["H", "He4"]
    density = np.random.random((rhogrid.size, 2)) * unyt.m**(-3)
    temperature = np.random.random((rhogrid.size, 2)) * unyt.eV
    obj = ascot.data.create_neutralradial(
        rhogrid=rhogrid,
        species=species,
        density=density,
        temperature=temperature,
    )

    assert np.allclose(obj.rhogrid, rhogrid)
    assert obj.species == species
    assert np.allclose(obj.density, density)
    assert np.allclose(obj.temperature, temperature)


def test_neutralarbitrary_init():
    ascot = Ascot()

    rgrid = np.linspace(0.1, 1.0, 10) * unyt.m
    zgrid = np.linspace(-1.0, 1.0, 11) * unyt.m
    phigrid = np.linspace(0., 360, 12)[:-1] * unyt.deg
    species = ["H", "He4"]
    density = np.random.random((rgrid.size, phigrid.size, zgrid.size, 2)) * unyt.m**(-3)
    temperature = np.random.random((rgrid.size, phigrid.size, zgrid.size, 2)) * unyt.eV
    obj = ascot.data.create_neutralarbitrary(
        rgrid=rgrid,
        phigrid=phigrid,
        zgrid=zgrid,
        species=species,
        density=density,
        temperature=temperature,
    )

    assert np.allclose(obj.rgrid, rgrid)
    assert np.allclose(obj.phigrid, phigrid)
    assert np.allclose(obj.zgrid, zgrid)
    assert obj.species == species
    assert np.allclose(obj.density, density)
    assert np.allclose(obj.temperature, temperature)
