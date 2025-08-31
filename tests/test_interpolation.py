import ctypes

import pytest
import unyt
import numpy as np

from a5py import Ascot
from a5py.libascot import LIBASCOT
from a5py.templates import PremadeMagneticField
from a5py.data.bfield import Bfield

import a5py.engine.functions

def test_interpolation():

    a5 = Ascot()
    template = PremadeMagneticField(a5, field="iter-baseline")
    template.create_input()
    a5.data.bfield.active.stage()

    bfield = Bfield()
    bfield.use(a5.data.bfield.active)

    r = np.array([6.2])
    phi = np.array([0.0])
    z = np.array([0.0])
    t = np.array([0.0])
    b = (np.zeros((3,1), dtype="f8") + np.nan) * unyt.T
    bjac = (np.zeros((9,1), dtype="f8") + np.nan) * unyt.T
    LIBASCOT.libascot_interpolate(
        ctypes.pointer(bfield), None, None, None, None, None, 1, -1,
        r, phi, z, t, b, bjac, *([None]*39)
        )
    print(b[:,0])
    print(bjac[:,0])
    assert 1 == 2
