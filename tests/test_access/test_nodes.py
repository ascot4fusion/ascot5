"""Tests for ~a5py.data.access.nodes classes.

TreeManager in these tests is a mock object.
"""
import pytest

import difflib
import textwrap

from a5py.exceptions import AscotIOException
from a5py.data.access import Leaf
from a5py.data.access.nodes import ImmutableNode, InputCategory, OutputLeaf

from .conftest import QIDS, DATES, INPUTS, OUTPUTS, CATEGORIES, NOTE


@pytest.fixture()
def node():
    """Create an empty ImmutableNode."""
    return ImmutableNode()


@pytest.fixture
def leaf():
    """Create single Leaf instance."""
    leaf = Leaf(qid=QIDS[0], date=DATES[0], variant=INPUTS[0], note=NOTE)
    return leaf


def test_initialization(node):
    """Test that the initialized node is unfrozen and empty."""
    assert not node._frozen
    assert not len(node._qids)
    assert not len(node._tags)
    assert node._active is None


def test_freeze_and_unfreeze(node):
    """Test that the node is immutable when frozen and mutable when unfrozen."""
    VALUE = object()
    node._freeze()
    assert node._frozen
    with pytest.raises(AscotIOException):
        node.new_attr = VALUE
    with pytest.raises(AttributeError):
        node.new_attr

    node._unfreeze()
    assert not node._frozen
    node.new_attr = VALUE
    assert node.new_attr is VALUE


def test_modify_attributes_context(node):
    """Test that the node can can be modified inside modify_attributes context.
    """
    VALUE = object()
    node._freeze()
    with node._modify_attributes():
        node.new_attr = VALUE
    assert node.new_attr is VALUE


def test_setgetitem(node):
    """Test that the node attributes can be assigned and referenced both with
    getattr and getitem methods, and set with either setattr or setitem
    (unless the node is frozen in which case neither should work).
    """
    VALUE, ANOTHER_VALUE, YET_ANOTHER_VALUE = object(), object(), object()
    node["attr"] = VALUE
    assert node.attr is VALUE
    assert node["attr"] is VALUE

    node.attr = ANOTHER_VALUE
    assert node.attr is ANOTHER_VALUE
    assert node["attr"] is ANOTHER_VALUE

    node._freeze()
    with pytest.raises(AscotIOException):
        node.attr = YET_ANOTHER_VALUE
    with pytest.raises(AscotIOException):
        node["attr"] = YET_ANOTHER_VALUE
    assert node.attr is ANOTHER_VALUE
    assert node["attr"] is ANOTHER_VALUE


def test_add_leaf(node, leaf):
    """Test that a leaf added to the node can be queried using string 'q+QID',
    leaf name, or leaf tag.
    """
    node._add_leaf(leaf)
    assert leaf is node[leaf.qqid]
    assert leaf is node[leaf.name]
    assert leaf is node[leaf._extract_tag()]
    assert node._qids == [QIDS[0]]


def test_contains(node, leaf):
    """Test that if the node contains a leaf can be checked with leaf's QID,
    'q+QID', name, tag, and with the leaf itself, and that the node can be
    iterated over.
    """
    node._add_leaf(leaf)
    assert leaf.qid in node
    assert leaf.qqid in node
    assert leaf.name in node
    assert leaf._extract_tag() in node
    assert leaf in node
    for leaf_in_node in node:
        assert leaf_in_node is leaf


def test_remove_leaf(node, leaf):
    """Test that removing a leaf removes all references to it."""
    node._add_leaf(leaf)
    node._remove_leaf(leaf)
    assert not leaf.qid in node
    assert not leaf.name in node
    assert not leaf._extract_tag() in node
    assert not leaf in node
    assert not len(node._qids)


def test_active_property(node):
    """Test that active leaf is by default the first leaf, but removing that
    leaf sets the nextly added leaf active until all leaves are removed.
    """
    leaf1 = Leaf(qid=QIDS[0], date=DATES[0], variant=INPUTS[0], note=NOTE)
    leaf2 = Leaf(qid=QIDS[1], date=DATES[0], variant=INPUTS[0], note=NOTE)
    leaf3 = Leaf(qid=QIDS[2], date=DATES[0], variant=INPUTS[0], note=NOTE)
    node._add_leaf(leaf1)
    assert node.active is leaf1

    node._add_leaf(leaf2)
    node._add_leaf(leaf3)
    assert node.active is leaf1

    node._remove_leaf(leaf1)
    assert node.active is leaf2

    node._remove_leaf(leaf2)
    node._remove_leaf(leaf3)
    with pytest.raises(AscotIOException):
        _ = node.active


def test_organize_by_date(node):
    """Test that leafs remain organized by date."""
    leaf1 = Leaf(qid=QIDS[0], date=DATES[1], variant=INPUTS[0], note="<SATURDAY>")
    leaf2 = Leaf(qid=QIDS[1], date=DATES[0], variant=INPUTS[0], note="<FRIDAY>")
    leaf3 = Leaf(qid=QIDS[2], date=DATES[2], variant=INPUTS[0], note="<SUNDAY>")
    node._add_leaf(leaf1)
    assert node._qids == [QIDS[0]]
    node._add_leaf(leaf2)
    assert node._qids == [QIDS[0], QIDS[1]]
    node._add_leaf(leaf3)
    assert node._qids == [QIDS[2], QIDS[0], QIDS[1]]
    assert node._tags == ["FRIDAY", "SATURDAY", "SUNDAY"]
    node._remove_leaf(leaf3)
    assert node._qids == [QIDS[0], QIDS[1]]
    assert node._tags == ["FRIDAY", "SATURDAY"]
    node._remove_leaf(leaf2)
    assert node._qids == [QIDS[0]]


def test_organize_same_tags(node):
    """Test that tags are properly renamed when there are multiple identical
    tags.
    """
    leaf1 = Leaf(qid=QIDS[0], date=DATES[1], variant=INPUTS[0], note="<TAG>")
    leaf2 = Leaf(qid=QIDS[1], date=DATES[0], variant=INPUTS[0], note="<TAG2>")
    leaf3 = Leaf(qid=QIDS[2], date=DATES[2], variant=INPUTS[0], note="<TAG>")

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

    inputnode._add_leaf(
        Leaf(qid=QIDS[0], date=DATES[1], note=NOTE, variant=INPUTS[0])
        )
    inputnode._add_leaf(
        Leaf(qid=QIDS[1], date=DATES[0], note=NOTE, variant=INPUTS[1])
        )

    expected = textwrap.dedent(
        """
        inputA     2991786571 1997-08-30 02:14:00 [active]
        BENNETT_0
        Let off some steam <BENNETT>

        inputB     9753987342 1997-08-29 02:14:00
        BENNETT_1
        Let off some steam <BENNETT>

        """)
    diff = '\n'.join(difflib.unified_diff(
        inputnode.contents.splitlines(),
        expected.splitlines()[1:],
        fromfile="contents", tofile="expected", lineterm="",
    ))
    assert not diff, f"Strings differ:\n{diff}"


def test_outputleaf_inputs():
    """Test that the outputleaf contains references to inputs, but trying to
    query an unused input raises exception.
    """
    leafin = Leaf(qid=QIDS[0], date=DATES[0], note=NOTE, variant=INPUTS[0])
    leafout = OutputLeaf(
        qid=QIDS[1], date=DATES[0], note=NOTE, variant=OUTPUTS[0],
        inputs={CATEGORIES[0]: leafin},
        )
    assert leafout[CATEGORIES[0]] is leafin

    with pytest.raises(AscotIOException):
        getattr(leafout, CATEGORIES[1])

def test_outputleaf_contents():
    """Test that the contents of an simulation output are displayed
    correctly.

    Note the lowercase "tag" as there's no actual tag yet (the leaves are not
    part of a tree).
    """
    leafin = Leaf(qid=QIDS[0], date=DATES[0], note=NOTE, variant=INPUTS[0])
    leafout = OutputLeaf(
        qid=QIDS[1], date=DATES[0], note=NOTE, variant=OUTPUTS[0],
        inputs={CATEGORIES[0]: leafin},
        )
    expected = textwrap.dedent(
        """
        output     9753987342 1997-08-29 02:14:00
        Let off some steam <Bennett>

        Inputs:
        catX    inputA     2991786571 1997-08-29 02:14:00
                Let off some steam <Bennett>
        """)
    diff = '\n'.join(difflib.unified_diff(
        leafout.contents.splitlines(),
        expected.splitlines()[1:],
        fromfile="contents", tofile="expected", lineterm="",
    ))
    assert not diff, f"Strings differ:\n{diff}"
