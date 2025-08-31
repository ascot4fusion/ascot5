"""Test physlib.species."""
import pytest

import a5py.physlib as physlib
from a5py.exceptions import AscotUnitWarning


@pytest.mark.parametrize("charge", [None, 1, 1.5])
def test_species2properties(charge):
    """Test the species2properties function."""
    if charge == 1.5:
        with pytest.raises(ValueError):
            physlib.species2properties("H", charge=charge)
        return
    properties = physlib.species2properties("H", charge=charge)

    assert properties.anum == physlib.KNOWN_SPECIES["H1"].anum
    assert properties.znum == physlib.KNOWN_SPECIES["H1"].znum
    assert properties.mass == physlib.KNOWN_SPECIES["H1"].mass

    if charge is None:
        assert properties.charge == physlib.KNOWN_SPECIES["H1"].charge
    else:
        assert properties.charge == charge


@pytest.mark.parametrize("anumznum", [(1, 1), (-1, 1)])
def test_findmass(anumznum):
    """Test the findmass function."""
    if anumznum == (-1, 1):
        with pytest.raises(ValueError):
            physlib.findmass(*anumznum)
        return
    assert physlib.findmass(*anumznum) == physlib.KNOWN_SPECIES["H1"].mass


@pytest.mark.parametrize("anumznum", [(1, 1), (-1, 1)])
def test_properties2species(anumznum):
    """Test the properties2species function."""
    if anumznum == (-1, 1):
        with pytest.raises(ValueError):
            physlib.findmass(*anumznum)
        return
    assert physlib.properties2species(*anumznum) == "H1"
