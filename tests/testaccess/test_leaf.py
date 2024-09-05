"""Tests for `Leaf` class."""
import pytest

from a5py.data.access.metadata import DEFAULT_TAG
from a5py.exceptions import AscotIOException

from .conftest import (
    QID1, QID2, DATE, NOTE, INPUTVAR, create_leaf,
    )

def test_initialization():
    """Test that attributes were set properly in initialization."""
    leaf = create_leaf(QID1)
    for actual, expected in zip(
        [leaf.qid, leaf.date, leaf.note, leaf.variant],
        [QID1, DATE, NOTE, INPUTVAR],
        ):
        assert actual == expected
    assert leaf.qqid == f"q{QID1}"
    assert leaf.name == f"{INPUTVAR}_{QID1}"


def test_modify_note():
    """Test that description can be changed."""
    leaf = create_leaf(QID1)
    leaf.note = "New note"
    assert leaf.note == "New note"


def test_cannot_change_fixed_attributes():
    """Test that immutable attributes are immutable."""
    leaf = create_leaf(QID1)
    for attr in ["qid", "date", "variant", "qqid", "name"]:
        with pytest.raises(AscotIOException):
            setattr(leaf, attr, QID2)


def test_extract_tag():
    """Test extracting a tag from description."""
    leaf = create_leaf(QID1)
    leaf.note = "get to the choppa"
    assert leaf._extract_tag() == "GET"
    leaf.note = "12345note with numbers"
    assert leaf._extract_tag() == DEFAULT_TAG
    leaf.note = "!no@#$special%^&*()characters"
    assert leaf._extract_tag() == "NOSPECIALCHARACTERS"


def test_methods_trigger_treemanager():
    """Test that the public methods trigger tree manager."""
    leaf = create_leaf(QID1)
    treemanager = leaf._treemanager
    leaf.activate()
    treemanager.activate_leaf.assert_called_with(leaf)
    leaf.destroy(repack=True)
    treemanager.destroy_leaf.assert_called_with(leaf, repack=True)
    leaf.note = "Because I'm going to say... please"
    treemanager.note_changed.assert_called_with(leaf)
