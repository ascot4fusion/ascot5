import unyt
import pytest
import numpy as np

from a5py import Ascot
from a5py.engine.interpolate import evaluate


def test_boozer_init():
    ascot = Ascot()

    psigrid = np.linspace(0, 1, 10)
    boozerpoloidal = np.random.random((10, 8)) * unyt.rad
    boozertoroidal = np.random.random((10, 12)) * unyt.rad
    separatrix = np.array([[0., 0., 1., 1.], [0., 1., 1., 0.]]) * unyt.m
    obj = ascot.data.create_boozer(
        psigrid=psigrid,
        boozerpoloidal=boozerpoloidal,
        boozertoroidal=boozertoroidal,
        separatrix=separatrix,
    )
    assert np.allclose(obj.psigrid, psigrid)
    assert np.allclose(obj.separatrix, separatrix)
    assert np.allclose(obj.boozerpoloidal, boozerpoloidal)
    assert np.allclose(obj.boozertoroidal, boozertoroidal)
