"""Test physlib.species."""
import unyt
import pytest

from a5py.physlib import Species

@pytest.mark.parametrize(
        "name, valid",
        [("D", True), ("Deuterium", True), ("Peruna", False)]
        )
def test_from_string(name, valid):
    """Test obtaining species by name."""
    if not valid:
        with pytest.raises(ValueError):
            Species.from_string(name)
        return

    sp = Species.from_string(name)
    assert sp.name == "D"
    assert sp.znum == 1
    assert sp.anum == 2
    assert sp.mass == 2.014*unyt.amu


@pytest.mark.parametrize("z, a, valid", [(1, 2, True), (999, 999, False)])
def test_from_znumanum(z, a, valid):
    """Test obtaining species by Z and A."""
    if not valid:
        with pytest.raises(ValueError):
            Species.from_znumanum(z, a)
        return

    sp = Species.from_znumanum(z, a)
    assert sp.name == "D"
    assert sp.znum == 1
    assert sp.anum == 2
    assert sp.mass == 2.014*unyt.amu

