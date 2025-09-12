# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests for the the tree and it's input category and run variant components."""
from __future__ import annotations

import os
import difflib
import textwrap

import pytest

from a5py import AscotIOException

from a5py.data.access import Leaf
from a5py.data.access import Tree
from a5py.data.access.hdf5 import RESULTGROUP, HDF5Interface
from a5py.data.access.nodes import OutputLeaf

from .conftest import (
    QIDS, INPUTS, OUTPUTS, CATEGORIES, DATE, NOTE, FNTEST, FNEMPTY, DATES
    )


def reinit_from_file(tree_to_reinit):
    """Reinitialize tree from file."""
    if tree_to_reinit._treemanager.hdf5manager is None:
        return tree_to_reinit
    file = tree_to_reinit._treemanager.hdf5manager.filename
    return Tree(CATEGORIES, (file, True))


@pytest.fixture
def tree(request):
    """Provide a Tree, optionally backed by an HDF5 file."""
    file = getattr(request, "param", None)
    if file:
        t = Tree(CATEGORIES, (FNEMPTY, False))
    else:
        t = Tree(CATEGORIES)

    yield t

    if file:
        try:
            os.unlink(FNEMPTY)
        except FileNotFoundError:
            pass


def test_init_nofile():
    """Test initializing a tree when no input file is provided.

    Tree-node itself should not contain any QIDs, but all category nodes should
    be present. Both the tree node and category nodes should be immutable.
    """
    tree = Tree(CATEGORIES)
    assert not len(tree._qids)
    with pytest.raises(AscotIOException):
        setattr(tree, "attribute", "value")

    for category in CATEGORIES:
        assert category in tree
        assert not len(tree[category]._qids)
        with pytest.raises(AscotIOException):
            setattr(tree[category], "attribute", "value")

    assert tree._treemanager.hdf5manager is None


def test_init_emptyfile():
    """Test initializing a tree with empty file.

    Tree-node itself should not contain any QIDs, but all category nodes should
    be present both in the instance and in the file.
    """
    try:
        tree = Tree(CATEGORIES, (FNEMPTY, False))
        with HDF5Interface(FNEMPTY, "r") as h5:
            for group in CATEGORIES + ["results"]:
                assert group in h5
            for group in h5.keys():
                assert group in CATEGORIES + ["results"]

        assert not len(tree._qids)
        for category in CATEGORIES:
            assert category in tree
            assert not len(tree[category]._qids)
    finally:
        os.unlink(FNEMPTY)


def test_init_fromfile():
    """Test initializing a tree from a file containing data.

    The file contains single input used by single output. Test that the metadata
    is read correctly and the nodes have correct active fields.
    """
    try:
        with HDF5Interface(FNTEST, "a") as h5:
            h5.set_node(CATEGORIES[0], QIDS[0])
            h5.set_node(RESULTGROUP, QIDS[1])
            h5.set_datagroup(
                CATEGORIES[0], f"{INPUTS[0]}_{QIDS[0]}", date=DATE, note=NOTE,
                )
            h5.set_datagroup(
                RESULTGROUP, f"{OUTPUTS[0]}_{QIDS[1]}", date=DATE, note=NOTE,
                **{CATEGORIES[0]:QIDS[0]},
            )
            h5[f"{CATEGORIES[0]}/{INPUTS[0]}_{QIDS[0]}"]

        tree = Tree(CATEGORIES, (FNTEST, True))
        assert tree[CATEGORIES[0]]._qids == [QIDS[0]]
        leafin = tree[CATEGORIES[0]][f"q{QIDS[0]}"]
        assert tree[CATEGORIES[0]].active is leafin
        assert leafin.date == DATE
        assert leafin.note == NOTE
        assert leafin.variant == INPUTS[0]

        assert tree._qids == [QIDS[1]]
        leafout = tree[f"q{QIDS[1]}"]
        assert leafout[CATEGORIES[0]] is leafin
        assert tree.active is leafout
    finally:
        os.unlink(FNTEST)


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_add_input(tree):
    """Test adding a new input."""
    leaf1 = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])

    tree._treemanager.enter_input(leaf1, CATEGORIES[0], dryrun=True)
    tree = reinit_from_file(tree)
    assert not leaf1.qid in tree[CATEGORIES[0]], "Dryrun was not dry."

    tree._treemanager.enter_input(leaf1, CATEGORIES[0])
    tree = reinit_from_file(tree)
    assert leaf1.qid in tree[CATEGORIES[0]], "Input was not added."

    leaf2 = Leaf(qid=QIDS[1], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf2, CATEGORIES[0])
    tree = reinit_from_file(tree)
    assert leaf1.qid == tree[CATEGORIES[0]].active.qid, (
        "Active input should not change unless explicitly activated.")

    leaf3 = Leaf(qid=QIDS[2], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf3, CATEGORIES[0], activate=True)
    tree = reinit_from_file(tree)
    assert leaf3.qid == tree[CATEGORIES[0]].active.qid, (
        "Active input was not changed even when explicitly activated.")


@pytest.mark.parametrize("tree", [True], indirect=True)
def test_add_input_and_save(tree):
    """Test adding new input without saving it, and then saving retroactively.
    """
    leaf1 = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf1, CATEGORIES[0], save=False)
    tree = reinit_from_file(tree)
    assert not leaf1.qid in tree[CATEGORIES[0]], (
        "Input was saved even though save=False")

    tree._treemanager.enter_input(leaf1, CATEGORIES[0], save=False)
    leaf1.save()
    tree = reinit_from_file(tree)
    assert leaf1.qid in tree[CATEGORIES[0]], (
        "Input was not saved by its save method.")


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_add_output(tree):
    """Test adding a new output."""
    leaf = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    leaf1 = OutputLeaf(
        qid=QIDS[3], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leaf})
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_output(leaf1, [QIDS[0]])

    tree._treemanager.enter_input(leaf, CATEGORIES[0])

    tree._treemanager.enter_output(leaf1, [QIDS[0]], dryrun=True)
    tree = reinit_from_file(tree)
    assert not leaf1.qid in tree, "Dryrun was not dry."

    tree._treemanager.enter_output(leaf1, [QIDS[0]])
    tree = reinit_from_file(tree)
    assert leaf1.qid in tree, "Output was not added."

    leaf2 = OutputLeaf(qid=QIDS[1], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leaf})
    tree._treemanager.enter_output(leaf2, [QIDS[0]])
    tree = reinit_from_file(tree)
    assert leaf1.qid == tree.active.qid, (
        "Active output should not change unless explicitly activated.")

    leaf3 = OutputLeaf(qid=QIDS[2], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leaf})
    tree._treemanager.enter_output(leaf3, [QIDS[0]], activate=True)
    tree = reinit_from_file(tree)
    assert leaf3.qid == tree.active.qid, (
        "Active input was not changed even when explicitly activated.")


@pytest.mark.parametrize("tree", [True], indirect=True)
def test_add_output_and_save(tree):
    """Test adding new output without saving it, and then saving retroactively.

    This test also checks that output cannot be saved unless input is saved
    first.
    """
    leafin = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leafin, CATEGORIES[0], save=False)

    leafout1 = OutputLeaf(qid=QIDS[1], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leafin})
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_output(leafout1, [QIDS[0]])

    assert leafout1.qid in tree, (
        "Output should be in tree even if save failed")

    leafin.save()
    leafout2 = OutputLeaf(qid=QIDS[2], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leafin})
    tree._treemanager.enter_output(leafout2, [QIDS[0]], save=False)
    tree = reinit_from_file(tree)
    assert not leafout2.qid in tree, "Output was saved even though save=False"

    tree._treemanager.enter_output(leafout2, [QIDS[0]], save=False)
    leafout2.save()
    tree = reinit_from_file(tree)
    assert leafout2.qid in tree, "Output was not saved by its save method."


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_add_identical_qid(tree):
    """Test adding data when one with identical QID exists."""
    leafin = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leafin, CATEGORIES[0])
    tree = reinit_from_file(tree)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_input(leafin, CATEGORIES[0])

    leaf = OutputLeaf(qid=QIDS[0], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leafin})
    tree = reinit_from_file(tree)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_output(leaf, [QIDS[0]])

    leaf = OutputLeaf(qid=QIDS[1], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leafin})
    tree._treemanager.enter_output(leaf, [QIDS[0]])
    tree = reinit_from_file(tree)
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_output(leaf, [QIDS[1]])


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_destroy(tree):
    """Test removing data."""
    leafin = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leafin, CATEGORIES[0])
    leafout = OutputLeaf(qid=QIDS[1], date=DATE, note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leafin})
    tree._treemanager.enter_output(leafout, [QIDS[0]])

    with pytest.raises(AscotIOException):
        tree.destroy(data=leafin)

    tree = reinit_from_file(tree)
    assert leafin.qid in tree[CATEGORIES[0]], "Used input was removed."

    tree.destroy(data=leafout)
    tree = reinit_from_file(tree)
    assert leafout not in tree

    tree.destroy(data=leafin)
    tree = reinit_from_file(tree)
    assert leafin not in tree[CATEGORIES[0]]


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_tree_activate_data(tree):
    """Test setting data active."""
    leaf1 = Leaf(qid=QIDS[0], date=DATES[0], note=NOTE, variant=INPUTS[0])
    leaf2 = Leaf(qid=QIDS[1], date=DATES[0], note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf1, CATEGORIES[0])
    tree._treemanager.enter_input(leaf2, CATEGORIES[0])
    assert tree[CATEGORIES[0]].active == leaf1
    tree.activate(leaf2)
    tree = reinit_from_file(tree)
    assert tree[CATEGORIES[0]].active.qid == leaf2.qid


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_tree_change_note(tree):
    """Test changing data's note."""
    leaf1 = Leaf(qid=QIDS[0], date=DATES[0], note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf1, CATEGORIES[0])
    leaf1.note = "Get to the choppa"
    tree = reinit_from_file(tree)
    assert tree[CATEGORIES[0]]["GET"].qid == leaf1.qid

    leaf2 = OutputLeaf(qid=QIDS[1], date=DATES[0], note=NOTE, variant=OUTPUTS[0], inputs={CATEGORIES[0]: leaf1})
    tree._treemanager.enter_output(leaf2, [QIDS[0]])
    leaf2.note = "Get to the choppa"
    tree = reinit_from_file(tree)
    assert tree["GET"].qid == leaf2.qid


@pytest.mark.parametrize("tree", [False], indirect=True)
def test_tree_contents(tree):
    """Test getting the contents in a string."""
    tree._treemanager.enter_input(
        Leaf(qid=QIDS[0], date=DATES[0], note=NOTE, variant=INPUTS[0]),
        CATEGORIES[0]
        )
    tree._treemanager.enter_input(
        Leaf(qid=QIDS[2], date=DATES[1], note=NOTE, variant=INPUTS[0]),
        CATEGORIES[0]
        )
    tree._treemanager.enter_input(
        Leaf(qid=QIDS[3], date=DATES[2], note=NOTE, variant=INPUTS[0]),
        CATEGORIES[1]
        )
    tree._treemanager.enter_output(
        Leaf(qid=QIDS[1], date=DATES[2], note=NOTE, variant=OUTPUTS[0]),
        [QIDS[0]]
        )
    tree._treemanager.enter_output(
        Leaf(qid=QIDS[4], date=DATES[2], note=NOTE, variant=OUTPUTS[0]),
        [QIDS[0]]
        )

    expected = textwrap.dedent(
        """
        Inputs: [only active shown]
        catX      inputA    2991786571 1997-08-29 02:14:00 + 1 other(s)
                  "Let off some steam Bennett"
        catY      inputA    0963810214 1997-08-31 02:14:00 (no other inputs)
                  "Let off some steam Bennett"
        catZ      *no inputs*


        Simulations:
        output    9753987342 1997-08-31 02:14:00 [active]
                  "Let off some steam Bennett"
        output    5960585966 1997-08-31 02:14:00
                  "Let off some steam Bennett"
        """)
    diff = '\n'.join(difflib.unified_diff(
        tree.contents.splitlines(),
        expected.splitlines()[1:],
        fromfile="contents", tofile="expected", lineterm="",
    ))
    assert not diff, f"Strings differ:\n{diff}"


@pytest.mark.parametrize("tree", [False], indirect=True)
def test_save_data_no_file(tree):
    """Test saving input or output when there is no file."""
    leaf = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    with pytest.raises(AscotIOException) as e:
        tree._treemanager.enter_input(leaf, CATEGORIES[0], save=True)
    assert "No HDF5 file was provided." in str(e.value)
    assert tree._treemanager.get_leaf(QIDS[0]).qid == QIDS[0]

    leaf = Leaf(qid=QIDS[2], date=DATE, note=NOTE, variant=OUTPUTS[0])
    with pytest.raises(AscotIOException):
        tree._treemanager.enter_output(leaf, [QIDS[0]], save=True)
    assert "No HDF5 file was provided." in str(e.value)
    assert tree._treemanager.get_leaf(QIDS[2]).qid == QIDS[2]


@pytest.mark.parametrize("tree", [True], indirect=True)
def test_set_active_not_on_file(tree):
    """Test that we don't store active QID of data which is not on disk."""
    leaf1 = Leaf(qid=QIDS[0], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf1, CATEGORIES[0])
    leaf2 = Leaf(qid=QIDS[1], date=DATE, note=NOTE, variant=INPUTS[0])
    tree._treemanager.enter_input(leaf2, CATEGORIES[0], save=False)
    leaf2.activate()

    assert tree[CATEGORIES[0]].active.qid == leaf2.qid

    tree = reinit_from_file(tree)
    assert tree[CATEGORIES[0]].active.qid == leaf1.qid
