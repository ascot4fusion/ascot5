import pytest

import unyt
import numpy as np

from a5py import Ascot
from a5py.templates import PremadeMagneticField, FlatPlasma

from a5py.engine import auxiliary


def test_rhotheta2rz():
    """Test rhotheta2rz."""
    a5 = Ascot()
    template = PremadeMagneticField(a5, field="iter-baseline")
    template.create_input()
    r, z = auxiliary.rhotheta2rz(
        a5.data.bfield.active,
        0.8,
        0.0,
        0.0,
        0.0,
        )
    print(r,z)
    assert False
