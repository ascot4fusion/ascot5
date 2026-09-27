import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.engine.interpolate import evaluate


def test_plasmalinear1d_init():
    ascot = Ascot()
    species = ["H", "He4"]
    charge = [1, 2]
    rhogrid = np.linspace(0.1, 1.0, 10)
    ni = np.random.random((rhogrid.size, len(species))) * unyt.m**(-3)
    Ti = np.random.random((rhogrid.size,)) * unyt.eV
    ne = np.random.random((rhogrid.size,)) * unyt.m**(-3)
    Te = np.random.random((rhogrid.size,)) * unyt.eV
    rotation = np.random.random((rhogrid.size,)) * unyt.rad / unyt.s
    obj = ascot.data.create_plasmalinear1d(
        species=species,
        rhogrid=rhogrid,
        ni=ni, ne=ne, Ti=Ti, Te=Te, charge=charge, rotation=rotation
    )
    assert np.allclose(obj.rhogrid, rhogrid)
    assert np.allclose(obj.ni, ni)
    assert np.allclose(obj.Ti, Ti)
    assert np.allclose(obj.ne, ne)
    assert np.allclose(obj.Te, Te)
    assert np.allclose(obj.charge, charge)
    assert np.allclose(obj.rotation, rotation)

def test_plasmalinear2d_init():
    ascot = Ascot()
    species = ["H", "He4"]
    charge = [1, 2]
    rgrid = np.linspace(0.1, 1.0, 10)
    zgrid = np.linspace(-1.0, 1.0, 12)
    ni = np.random.random((rgrid.size, zgrid.size, len(species))) * unyt.m**(-3)
    Ti = np.random.random((rgrid.size, zgrid.size)) * unyt.eV
    ne = np.random.random((rgrid.size, zgrid.size)) * unyt.m**(-3)
    Te = np.random.random((rgrid.size, zgrid.size)) * unyt.eV
    rotation = np.random.random((rgrid.size, zgrid.size)) * unyt.rad / unyt.s
    obj = ascot.data.create_plasmalinear2d(
        species=species,
        rgrid=rgrid, zgrid=zgrid,
        ni=ni, ne=ne, Ti=Ti, Te=Te, charge=charge, rotation=rotation
    )
    assert np.allclose(obj.rgrid, rgrid)
    assert np.allclose(obj.zgrid, zgrid)
    assert np.allclose(obj.ni, ni)
    assert np.allclose(obj.Ti, Ti)
    assert np.allclose(obj.ne, ne)
    assert np.allclose(obj.Te, Te)
    assert np.allclose(obj.charge, charge)
    assert np.allclose(obj.rotation, rotation)

def test_plasmadynamic1d_init():
    ascot = Ascot()
    species = ["H", "He4"]
    charge = [1, 2]
    rhogrid = np.linspace(0.1, 1.0, 10)
    timegrid = np.linspace(-1.0, 1.0, 12)
    ni = np.random.random((rhogrid.size, timegrid.size, len(species))) * unyt.m**(-3)
    Ti = np.random.random((rhogrid.size, timegrid.size)) * unyt.eV
    ne = np.random.random((rhogrid.size, timegrid.size)) * unyt.m**(-3)
    Te = np.random.random((rhogrid.size, timegrid.size)) * unyt.eV
    rotation = np.random.random((rhogrid.size, timegrid.size)) * unyt.rad / unyt.s
    obj = ascot.data.create_plasmadynamic1d(
        species=species,
        rhogrid=rhogrid, timegrid=timegrid,
        ni=ni, ne=ne, Ti=Ti, Te=Te, charge=charge, rotation=rotation
    )
    assert np.allclose(obj.rhogrid, rhogrid)
    assert np.allclose(obj.timegrid, timegrid)
    assert np.allclose(obj.ni, ni)
    assert np.allclose(obj.Ti, Ti)
    assert np.allclose(obj.ne, ne)
    assert np.allclose(obj.Te, Te)
    assert np.allclose(obj.charge, charge)
    assert np.allclose(obj.rotation, rotation)
