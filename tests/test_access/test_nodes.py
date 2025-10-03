"""Tests for ~a5py.data.access.nodes classes.

TreeManager in these tests is a mock object.
"""

import pytest

import difflib
import textwrap

from a5py.exceptions import AscotDataException, AscotMeltdownError
from a5py.data.access import Leaf
from a5py.data.access.nodes import ImmutableNode, InputCategory

from .conftest import DATES, NOTE


@pytest.fixture()
def node():
    """Create an empty ImmutableNode."""
    return ImmutableNode()


@pytest.fixture
def leaf():
    """Create single Leaf instance."""
    leaf = Leaf(date=DATES[0], note=NOTE)
    return leaf


def test_add_leaf(node, leaf):
    """Test that a leaf added to the node can be queried using leaf name or tag."""
    node._add_leaf(leaf)
    assert leaf is node[leaf.name]
    assert leaf is node[Leaf.extract_tag(leaf.note)[0]]
    assert node._names == [leaf.variant]
    with pytest.raises(AscotMeltdownError):
        node._add_leaf(leaf)


def test_contains(node, leaf):
    """Test that if the node contains a leaf can be checked with leaf's QID,
    'q+QID', name, tag, and with the leaf itself, and that the node can be
    iterated over.
    """
    node._add_leaf(leaf)
    assert leaf.name in node
    assert Leaf.extract_tag(leaf.note)[0] in node
    assert leaf in node
    for leaf_in_node in node:
        assert leaf_in_node is leaf


def test_remove_leaf(node, leaf):
    """Test that removing a leaf removes all references to it."""
    node._add_leaf(leaf)
    node._remove_leaf(leaf)
    assert not leaf.name in node
    assert not Leaf.extract_tag(leaf.note)[0] in node
    assert not leaf in node
    assert not len(node._names)


def test_active_property(node):
    """Test that active leaf is by default the first leaf, but removing that
    leaf sets the nextly added leaf active until all leaves are removed.
    """
    leaf1 = Leaf(date=DATES[0])
    leaf2 = Leaf(date=DATES[1])
    leaf3 = Leaf(date=DATES[2])
    leaf1._name += "_1"
    leaf2._name += "_2"
    leaf3._name += "_3"
    node._add_leaf(leaf1)
    assert node.active is leaf1

    node._add_leaf(leaf2)
    node._add_leaf(leaf3)
    assert node.active is leaf1

    node._remove_leaf(leaf1)
    assert node.active is leaf3

    node._remove_leaf(leaf2)
    node._remove_leaf(leaf3)
    with pytest.raises(AscotDataException):
        _ = node.active


def test_organize_by_date(node):
    """Test that leafs remain organized by date."""
    leaf1 = Leaf(date=DATES[1], note="<SATURDAY>")
    leaf2 = Leaf(date=DATES[0], note="<FRIDAY>")
    leaf3 = Leaf(date=DATES[2], note="<SUNDAY>")
    leaf1._name += "_1"
    leaf2._name += "_2"
    leaf3._name += "_3"
    node._add_leaf(leaf1)
    assert node._names == [f"{Leaf.__name__}_1"]
    node._add_leaf(leaf2)
    assert node._names == [f"{Leaf.__name__}_1", f"{Leaf.__name__}_2"]
    node._add_leaf(leaf3)
    assert node._names == [
        f"{Leaf.__name__}_3",
        f"{Leaf.__name__}_1",
        f"{Leaf.__name__}_2",
    ]
    assert node._tags == ["FRIDAY", "SATURDAY", "SUNDAY"]
    node._remove_leaf(leaf3)
    assert node._names == [f"{Leaf.__name__}_1", f"{Leaf.__name__}_2"]
    assert node._tags == ["FRIDAY", "SATURDAY"]
    node._remove_leaf(leaf2)
    assert node._names == [f"{Leaf.__name__}_1"]


def test_organize_same_tags(node):
    """Test that tags are properly renamed when there are multiple identical
    tags.
    """
    leaf1 = Leaf(date=DATES[1], note="<TAG>")
    leaf2 = Leaf(date=DATES[0], note="<TAG2>")
    leaf3 = Leaf(date=DATES[2], note="<TAG>")
    leaf1._name += "_1"
    leaf2._name += "_2"
    leaf3._name += "_3"

    node._add_leaf(leaf1)
    node._add_leaf(leaf2)
    assert "TAG" in node
    assert "TAG2" in node

    node._add_leaf(leaf3)
    assert not "TAG" in node
    assert "TAG_0" in node
    assert "TAG_1" in node
    assert "TAG2" in node
    assert node.TAG_0 is leaf3
    assert node.TAG_1 is leaf1
    assert node.TAG2 is leaf2

    node._remove_leaf(leaf3)
    assert "TAG" in node
    assert not "TAG_0" in node
    assert not "TAG_1" in node
    assert node.TAG is leaf1


def test_input_category_contents():
    """Test that the contents of an input category are displayed
    correctly.
    """
    inputnode = InputCategory()
    assert inputnode.contents == "No data in this category.\n"

    leaf1 = Leaf(date=DATES[1], note=NOTE)
    leaf2 = Leaf(date=DATES[0], note=NOTE)
    leaf1._name += "_1"
    leaf2._name += "_2"
    inputnode._add_leaf(leaf1)
    inputnode._add_leaf(leaf2)

    expected = textwrap.dedent(
        """
        Leaf_1          1997-08-30 02:14:00 [active]
        BENNETT_0
        Let off some steam <Bennett>

        Leaf_2          1997-08-29 02:14:00
        BENNETT_1
        Let off some steam <Bennett>

        """
    )
    diff = "\n".join(
        difflib.unified_diff(
            inputnode.contents.splitlines(),
            expected.splitlines()[1:],
            fromfile="contents",
            tofile="expected",
            lineterm="",
        )
    )
    assert not diff, f"Strings differ:\n{diff}"
