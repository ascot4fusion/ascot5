import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.engine.interpolate import evaluate


def test_efieldcartesian_init():
    ascot = Ascot()

    exyz = np.array([1., 2., 3.]) * unyt.V/unyt.m
    obj = ascot.data.create_efieldcartesian(
        exyz=exyz,
        )
    assert np.allclose(obj.exyz, exyz)


def test_efieldpotential1d_init():
    ascot = Ascot()

    rhogrid = np.linspace(0, 1, 5)
    dvdrho = np.array([1., 2., 3., 4., 5.]) * unyt.V
    obj = ascot.data.create_efieldpotential1d(
        rhogrid=rhogrid,
        dvdrho=dvdrho,
        )
    assert np.allclose(obj.rhogrid, rhogrid)
    assert np.allclose(obj.dvdrho, dvdrho)

def test_efieldpotential2d_init():
    ascot = Ascot()

    rgrid = np.linspace(0, 1, 5) * unyt.m
    zgrid = np.linspace(-1, 1, 4) * unyt.m
    potential = np.random.random((5, 4)) * unyt.V
    obj = ascot.data.create_efieldpotential2d(
        rgrid=rgrid,
        zgrid=zgrid,
        potential=potential,
        )
    assert np.allclose(obj.rgrid, rgrid)
    assert np.allclose(obj.zgrid, zgrid)
    assert np.allclose(obj.potential, potential)
