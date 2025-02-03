"""Tests for `ImmutableNode` class.

These tests also cover the `ImmutableStorage` class implicitly.
"""
import pytest

from unittest.mock import MagicMock

from a5py.exceptions import AscotIOException
from a5py.data.access.treeparts import ImmutableNode

from .conftest import (
    QID1, QID2, QID3, DATE_FRI, DATE_SAT, DATE_SUN, create_leaf,
    )


@pytest.fixture(name="treemanager")
def fixture_treemanager():
    return MagicMock()


@pytest.fixture(name="node")
def fixture_node(treemanager):
    """Set a mock tree manager that tests may use."""
    node = ImmutableNode()
    node._treemanager = treemanager
    return node


def test_initialization(node):
    """Test that the initialized object is unfrozen and is empty."""
    assert not node._frozen
    assert node._qids == []
    assert node._tags == []
    assert node._active is None


def test_freeze_and_unfreeze(node):
    """Test immutability."""
    node._freeze()
    assert node._frozen
    with pytest.raises(AscotIOException):
        node.new_attr = "value"

    node._unfreeze()
    assert not node._frozen
    node.new_attr = "value"
    assert node.new_attr == "value"


def test_modify_attributes_context(node):
    """Test the modify attributes context."""
    node._freeze()
    with node._modify_attributes():
        node.new_attr = "value"
    assert node.new_attr == "value"


def test_setgetitem(node):
    """Test the dictionary-like attribute access."""
    node["attr"] = "value"
    assert node.attr == "value"
    assert node["attr"] == "value"
    node._freeze()
    with pytest.raises(AscotIOException):
        node["attr"] = "another_value"


def test_add_leaf(node):
    """Test adding a single leaf."""
    leaf = create_leaf(QID1)
    node._add_leaf(leaf)
    assert leaf == node[leaf.qqid]
    assert leaf == node[leaf.name]
    assert leaf == node[leaf._extract_tag()]
    assert node._qids == [QID1]


def test_contains(node):
    """Test that the contained leaf can be queried in various ways."""
    leaf = create_leaf(QID1)
    node._add_leaf(leaf)
    assert leaf.qid in node
    assert leaf.qqid in node
    assert leaf.name in node
    assert leaf._extract_tag() in node
    assert leaf in node


def test_remove_leaf(node):
    """Test removing a single leaf."""
    leaf = create_leaf(QID1)
    node._add_leaf(leaf)
    node._remove_leaf(leaf)
    assert leaf.qid not in node
    assert leaf.name not in node
    assert leaf._extract_tag() not in node
    assert leaf not in node
    assert node._qids == []


def test_active_property(node):
    """Test switching active leaf."""
    leaf1 = create_leaf(QID1)
    leaf2 = create_leaf(QID2)
    node._add_leaf(leaf1)
    node._add_leaf(leaf2)
    assert node.active == leaf1
    node._remove_leaf(leaf1)
    assert node.active == leaf2
    node._remove_leaf(leaf2)
    with pytest.raises(AscotIOException):
        _ = node.active


def test_organize_by_date(node):
    """Test that leafs remain organized by date."""
    leaf1 = create_leaf(QID1, date=DATE_SAT, note="Saturday")
    leaf2 = create_leaf(QID2, date=DATE_FRI, note="Friday")
    leaf3 = create_leaf(QID3, date=DATE_SUN, note="Sunday")
    node._add_leaf(leaf1)
    assert node._qids == [QID1]
    node._add_leaf(leaf2)
    assert node._qids == [QID1, QID2]
    node._add_leaf(leaf3)
    assert node._qids == [QID3, QID1, QID2]
    assert node._tags == ["FRIDAY", "SATURDAY", "SUNDAY"]
    node._remove_leaf(leaf3)
    assert node._qids == [QID1, QID2]
    assert node._tags == ["FRIDAY", "SATURDAY"]
    node._remove_leaf(leaf2)
    assert node._qids == [QID1]


def _test_organize_same_tags(node):
    """Test that tags are properly renamed when there are multiple identical
    tags"""
    leaf1 = create_leaf(QID1, date=DATE_SAT, note="TAG")
    leaf2 = create_leaf(QID2, date=DATE_FRI, note="TAG2")
    leaf3 = create_leaf(QID3, date=DATE_SUN, note="TAG")
    node._add_leaf(leaf1)
    node._add_leaf(leaf2)
    assert "TAG" in node
    assert "TAG2" in node
    node._add_leaf(leaf3)
    assert "TAG" not in node
    assert "TAG_0" in node
    assert "TAG_1" in node
    assert "TAG2" in node

    assert node.TAG_0 == leaf3
    assert node.TAG_1 == leaf1
    assert node.TAG2 == leaf2

    node._remove_leaf(leaf3)
    assert "TAG" in node
    assert "TAG_0" not in node
    assert "TAG_1" not in node
    assert node.TAG == leaf1


def _test_update_tag(node):
    """Test that the tag is updated when the note is changed."""
    leaf = create_leaf(QID1, note="tag")
    node._add_leaf(leaf)
    assert "TAG" in node

    leaf.note = "Newtag"
    assert "NEWTAG" in node