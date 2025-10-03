"""Tests for a5py.utils.imstorage module."""
import pytest

from a5py.utils import ImmutableStorage

@pytest.fixture()
def icecube():
    """Create an empty ImmutableStorage."""
    return ImmutableStorage()


def test_initialization(icecube):
    """Test that the initialized storage is unfrozen."""
    assert not icecube._frozen


def test_freeze_and_unfreeze(icecube):
    """Test that the storage is immutable when frozen and mutable when unfrozen.
    """
    VALUE = object()
    icecube._freeze()
    assert icecube._frozen
    with pytest.raises(AttributeError):
        icecube.new_attr = VALUE
    with pytest.raises(AttributeError):
        icecube.new_attr

    icecube._unfreeze()
    assert not icecube._frozen
    icecube.new_attr = VALUE
    assert icecube.new_attr is VALUE


def test_modify_attributes_context(icecube):
    """Test that the storage can can be modified inside modify_attributes
    context.
    """
    VALUE = object()
    icecube._freeze()
    with icecube._modify_attributes():
        icecube.new_attr = VALUE
    assert icecube.new_attr is VALUE


def test_setgetitem(icecube):
    """Test that the storage attributes can be assigned and referenced both with
    getattr and getitem methods, and set with either setattr or setitem
    (unless the storage is frozen in which case neither should work).
    """
    VALUE, ANOTHER_VALUE, YET_ANOTHER_VALUE = object(), object(), object()
    icecube["attr"] = VALUE
    assert icecube.attr is VALUE
    assert icecube["attr"] is VALUE

    icecube.attr = ANOTHER_VALUE
    assert icecube.attr is ANOTHER_VALUE
    assert icecube["attr"] is ANOTHER_VALUE

    icecube._freeze()
    with pytest.raises(AttributeError):
        icecube.attr = YET_ANOTHER_VALUE
    with pytest.raises(AttributeError):
        icecube["attr"] = YET_ANOTHER_VALUE
    assert icecube.attr is ANOTHER_VALUE
    assert icecube["attr"] is ANOTHER_VALUE
