"""Tests for a5py.data.access.Leaf class."""
import pytest
from unittest.mock import MagicMock

from a5py.exceptions import AscotIOException
from a5py.data.access import Leaf
from a5py.data.access.leaf import DEFAULT_TAG

from .conftest import QIDS, NOTE, DATE, INPUTS


@pytest.fixture
def leaf():
    """Create a Leaf instance with mock tree manager."""
    leaf = Leaf(qid=QIDS[0], date=DATE, variant=INPUTS[0], note=NOTE)
    leaf._treemanager = MagicMock()
    return leaf


def test_initialization(leaf):
    """Test that attributes were set properly in initialization."""
    for actual, expected in zip(
        [leaf.qid, leaf.date, leaf.note, leaf.variant],
        [QIDS[0], DATE, NOTE, INPUTS[0]],
        ):
        assert actual == expected
    assert leaf.qqid == f"q{QIDS[0]}"
    assert leaf.name == f"{INPUTS[0]}_{QIDS[0]}"


def test_modify_note(leaf):
    """Test that the note can be changed."""
    leaf.note = "New note"
    assert leaf.note == "New note"


def test_cannot_change_fixed_attributes(leaf):
    """Test that immutable attributes are immutable."""
    for attr in ["qid", "date", "variant", "qqid", "name"]:
        with pytest.raises(AscotIOException):
            setattr(leaf, attr, QIDS[0])


@pytest.mark.parametrize(
    "note, tag",
    [
        ("get to the choppa", "GET"),
        ("12345note with numbers", DEFAULT_TAG),
        ("!no@#$special%^&*()characters", "NOSPECIALCHARACTERS"),
    ],
)
def test_extract_tag(leaf, note, tag):
    """Test extracting a tag from the note."""
    leaf.note = note
    assert leaf._extract_tag() == tag


def test_methods_trigger_treemanager(leaf):
    """Test that the public methods trigger tree manager."""
    treemanager = leaf._treemanager
    leaf.activate()
    treemanager.activate_leaf.assert_called_with(leaf)
    leaf.destroy(repack=True)
    treemanager.destroy_leaf.assert_called_with(leaf, repack=True)
    leaf.note = "Because I'm going to say... please"
    treemanager.note_changed.assert_called_with(leaf)
