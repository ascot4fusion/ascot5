# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests for the the tree and it's input category and run variant components."""
from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from a5py import AscotIOException

from a5py.data.access.metadata import input_categories, MetaData
from a5py.data.access.tree import InputCategory, RunVariant, Tree
from a5py.data.access.treeparts import Leaf

from .conftest import (
    QID1, QID2, QID3, QID4, QID5, DATE, NOTE, DATE_FRI, DATE_SAT, INPUTVAR,
    INPUTVAR2, OUTPUTVAR, CATEGORY, DIAGNOSTIC, create_leaf,
    )


@pytest.fixture(name="tree")
def fixture_tree():
    """Create empty tree."""
    return Tree()


@pytest.fixture(name="diagnostic")
def fixture_diagnostic():
    """Create mock diagnostic."""
    return MagicMock()


@pytest.fixture(name="inputA")
def fixture_inputA():
    """Create input of category A."""
    return Leaf(
        qid=QID1,
        date=DATE,
        note=NOTE,
        variant=INPUTVAR,
    )


@pytest.fixture(name="inputB")
def fixture_inputB():
    """Create input of category B."""
    return Leaf(
        qid=QID2,
        date=DATE,
        note=NOTE,
        variant=INPUTVAR2,
    )


@pytest.fixture(name="runvariant")
def fixture_runvariant(diagnostic, inputA, inputB):
    """Create a run variant."""
    return RunVariant(
        inputs={"wall":inputA, "bfield":inputB},
        diagnostics={DIAGNOSTIC:diagnostic},
        qid=QID3,
        date=DATE,
        note=NOTE,
        variant="run",
        )


def test_input_category_contents():
    """Test that the contents of an input category are displayed
    correctly.
    """
    inputnode = InputCategory()
    leaf1 = create_leaf(QID1, date=DATE_SAT)
    leaf2 = create_leaf(QID2, date=DATE_FRI)
    assert inputnode.contents == "No data in this category.\n"

    inputnode._add_leaf(leaf1)
    inputnode._add_leaf(leaf2)
    lines = inputnode.contents.splitlines()
    assert leaf1.variant in lines[0]
    assert leaf1.qid in lines[0]
    assert leaf1.date in lines[0]
    assert "active" in lines[0]
    assert inputnode._tags[0] in lines[1]
    assert leaf1.note in lines[2]
    assert leaf2.variant in lines[4]
    assert leaf2.qid in lines[4]
    assert leaf2.date in lines[4]
    assert "active" not in lines[4]
    assert inputnode._tags[1] in lines[5]
    assert leaf2.note in lines[6]


def test_run_variant_reference_inputs(runvariant, inputA):
    """Test that inputs are referenced correctly."""
    assert runvariant[CATEGORY] == inputA
    with pytest.raises(AscotIOException):
        _ = runvariant.efield


def test_run_variant_reference_diagnostics(runvariant, diagnostic):
    """Test that diagnostics are referenced correctly."""
    assert getattr(runvariant, f"_{DIAGNOSTIC}") == diagnostic
    with pytest.raises(AscotIOException):
        _ = runvariant._endstate


def test_run_variant_contents(runvariant):
    """Test that the contents of an simulation output are displayed
    correctly.
    """
    lines = runvariant.contents.splitlines()
    assert runvariant.variant in lines[0]
    assert runvariant.qid in lines[0]
    assert runvariant.date in lines[0]
    assert runvariant.note in lines[1]
    assert "Diagnostics" in lines[3]
    assert DIAGNOSTIC in lines[4]
    assert "Inputs" in lines[6]
    assert CATEGORY in lines[9]
    assert runvariant[CATEGORY].variant in lines[9]
    assert runvariant[CATEGORY].qid in lines[9]
    assert runvariant[CATEGORY].date in lines[9]
    assert runvariant[CATEGORY].note in lines[10]


def test_tree_initialization(tree):
    """Test initializing (empty) tree."""
    assert not len(tree._qids)
    with pytest.raises(AscotIOException):
        setattr(tree, "attribute", "value")

    for category in input_categories:
        assert category in tree
        with pytest.raises(AscotIOException):
            setattr(tree[category], "attribute", "value")


def test_tree_add_input(tree):
    """Test adding a new input."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    leaf = tree._treemanager.enter_input(meta, dryrun=True)
    assert leaf
    assert leaf not in tree[CATEGORY]
    leaf = tree._treemanager.enter_input(meta)
    assert leaf
    assert leaf in tree[CATEGORY]


def test_tree_add_simulation_output(tree):
    """Test adding a new run."""
    metaout = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_run(metaout, [], [QID1])

    metain = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    tree._treemanager.enter_input(metain)
    leaf = tree._treemanager.enter_run(metaout, [], [QID1], dryrun=True)
    assert leaf
    assert leaf not in tree
    leaf = tree._treemanager.enter_run(metaout, [], [QID1])
    assert leaf
    assert leaf in tree


def test_tree_add_identical_qid(tree):
    """Test adding data when one with identical QID exists."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    tree._treemanager.enter_input(meta)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_input(meta)

    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=OUTPUTVAR)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_run(meta, [], [QID1])

    meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
    tree._treemanager.enter_run(meta, [], [QID1])
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_run(meta, [], [QID2])


def test_tree_data(tree: Tree):
    """Test removing data."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    data = tree._treemanager.enter_input(meta)
    meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
    output = tree._treemanager.enter_run(meta, [], [QID1])

    with pytest.raises(AscotIOException):
        tree.destroy(data=data)

    tree.destroy(data=output)
    assert output not in tree
    tree.destroy(data=data)
    assert data not in tree[CATEGORY]


def test_tree_activate_data(tree):
    """Test setting data active."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    tree._treemanager.enter_input(meta)
    meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=INPUTVAR)
    data = tree._treemanager.enter_input(meta)
    tree.activate(data)
    assert tree[CATEGORY].active == data


def test_tree_change_note(tree):
    """Test changing data's note."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    leaf = tree._treemanager.enter_input(meta)
    leaf.note = "Get to the choppa"
    assert tree[CATEGORY]["GET"] == leaf

    meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
    leaf = tree._treemanager.enter_run(meta, [], [])
    leaf.note = "Get to the choppa"
    assert tree["GET"] == leaf


def test_tree_contents(tree):
    """Test getting the contents in a string."""
    meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
    data = tree._treemanager.enter_input(meta)
    meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=INPUTVAR)
    tree._treemanager.enter_input(meta)
    meta = MetaData(qid=QID4, date=DATE, note=NOTE, variant=INPUTVAR2)
    tree._treemanager.enter_input(meta)
    meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
    output = tree._treemanager.enter_run(meta, [], [QID1])
    meta = MetaData(qid=QID5, date=DATE, note=NOTE, variant=OUTPUTVAR)
    tree._treemanager.enter_run(meta, [], [QID1])

    lines = tree.contents.splitlines()
    assert "Inputs: [only active shown]" in lines[0]
    assert "*no inputs*" in lines[1]
    assert "(no other inputs)" in lines[5]
    assert CATEGORY in lines[13]
    assert data.variant in lines[13]
    assert data.qid in lines[13]
    assert data.date in lines[13]
    assert "+ 1 other(s)" in lines[13]
    assert data.note in lines[14]

    assert "Simulations:" in lines[24]
    assert output.variant in lines[25]
    assert output.qid in lines[25]
    assert output.date in lines[25]
    assert "[active]" in lines[25]
    assert output.note in lines[26]
    assert lines[27]
    assert "[active]" not in lines[27]
