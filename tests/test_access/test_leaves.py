"""Tests for a5py.data.access.leaves module."""

import pytest
import difflib
import textwrap
from unittest.mock import MagicMock

from a5py.exceptions import AscotDataException
from a5py.data.access import Leaf, InputVariant, OutputVariant, Status

from .conftest import NOTE, DATE, CATEGORIES


@pytest.fixture
def leaf():
    """Create a Leaf instance with mock tree manager."""
    leaf = Leaf(date=DATE, note=NOTE)
    leaf._treemanager = MagicMock()
    return leaf


@pytest.mark.parametrize(
    "attribute, expected",
    [
        ("name", Leaf.__name__),
        ("date", DATE),
        ("note", NOTE),
        ("variant", Leaf.__name__),
    ],
)
def test_leaf_initialization(leaf, attribute, expected):
    """Test that attributes were set properly in initialization."""
    assert getattr(leaf, attribute) == expected


def test_leaf_modify_note(leaf):
    """Test that the note can be changed."""
    leaf.note = "New note"
    assert leaf.note == "New note"


def test_leaf_cannot_change_fixed_attributes(leaf):
    """Test that immutable attributes are immutable."""
    for attr in ["date", "variant", "name"]:
        with pytest.raises(AttributeError):
            setattr(leaf, attr, "value")


@pytest.mark.parametrize(
    "note, tag",
    [
        ("get to the <choppa>", "CHOPPA"),
        ("12345note with numbers", None),
        ("<!no@#$5pec1al%^&*()character5>", "NO5PEC1ALCHARACTER5"),
        ("<1cannot start with number>", ValueError),
        ("<no> <multiple> tags", ValueError),
    ],
)
def test_leaf_extract_tag(leaf, note, tag):
    """Test extracting a tag from the note."""
    if tag == ValueError:
        with pytest.raises(ValueError):
            leaf.note = note
        return
    leaf.note = note
    assert Leaf.extract_tag(leaf.note)[0] == tag


def test_leaf_methods_trigger_treemanager(leaf):
    """Test that the public methods trigger tree manager."""
    treemanager = leaf._treemanager
    leaf.activate()
    treemanager.activate_leaf.assert_called_with(leaf)
    leaf.destroy(repack=True)
    treemanager.destroy_leaf.assert_called_with(leaf, repack=True)
    leaf.note = "Because I'm going to say... please"
    treemanager.note_changed.assert_called_with(leaf)
    leaf.save()
    treemanager.save_leaf.assert_called_with(leaf)


def test_leaf_status():
    """Test that the leaf status is reported correctly."""
    leaf = Leaf()
    assert leaf.status is Status.STAGED
    leaf._file = object()
    assert leaf.status == Status.SAVED


def test_input_variant_status():
    """Test that the input variant status is reported correctly."""
    leaf = InputVariant()
    with pytest.raises(AscotDataException):
        leaf.status
    leaf._file = object()
    assert leaf.status is Status.SAVED
    leaf._cdata = object()
    assert leaf.status & Status.STAGED
    assert leaf.status & Status.SAVED
    leaf._file = None
    assert leaf.status is Status.STAGED
    leaf._cdata = None
    with pytest.raises(AscotDataException):
        leaf.status


def test_input_variant_stage_unstage():
    """Test staging and unstaging of input variant."""
    leaf = InputVariant()
    leaf._file, leaf._cdata = None, None
    with pytest.raises(AscotDataException):
        leaf.stage()
    with pytest.raises(AscotDataException):
        leaf.unstage()

    leaf._file, leaf._cdata = object(), None
    leaf.stage()
    with pytest.raises(AscotDataException):
        leaf.unstage()

    leaf._file, leaf._cdata = object(), object()
    with pytest.raises(AscotDataException):
        leaf.stage()
    leaf.unstage()

    leaf._file, leaf._cdata = None, object()
    with pytest.raises(AscotDataException):
        leaf.unstage()


def test_outputleaf_inputs():
    """Test that the outputleaf contains references to inputs, but trying to
    query an unused input raises exception.
    """
    leafin = Leaf()
    leafout = OutputVariant(inputs={CATEGORIES[0]: leafin})
    leafout._treemanager = MagicMock()
    leafout._treemanager.nodes = {CATEGORIES[0]: object(), CATEGORIES[1]: object()}
    assert leafout[CATEGORIES[0]] is leafin

    with pytest.raises(AscotDataException):
        getattr(leafout, CATEGORIES[1])


def test_outputleaf_contents():
    """Test that the contents of an simulation output are displayed
    correctly.

    Note the lowercase "tag" as there's no actual tag yet (the leaves are not
    part of a tree).
    """
    leafin = Leaf(date=DATE, note=NOTE)
    leafout = OutputVariant(
        date=DATE,
        note=NOTE,
        inputs={CATEGORIES[0]: leafin},
    )
    expected = textwrap.dedent(
        """
        OutputVariant   1997-08-29 02:14:00
        Let off some steam <Bennett>

        Inputs:
        catX    Leaf            1997-08-29 02:14:00
                Let off some steam <Bennett>
        """
    )
    diff = "\n".join(
        difflib.unified_diff(
            leafout.contents.splitlines(),
            expected.splitlines()[1:],
            fromfile="contents",
            tofile="expected",
            lineterm="",
        )
    )
    assert not diff, f"Strings differ:\n{diff}"
