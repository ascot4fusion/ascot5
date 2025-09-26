# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests for the the tree and it's input category and run variant components."""
from __future__ import annotations

import os
import difflib
import textwrap

import pytest

from a5py.exceptions import AscotMeltdownError, AscotDataException

from a5py.data.access import InputVariant, OutputVariant
from a5py.data.access import Tree
from a5py.data.access.hdf5io import RESULTGROUP, TreeFile

from .conftest import (
    INPUTS, OUTPUTS, CATEGORIES, DATE, NOTE, FNTEST, FNEMPTY, DATES
    )

# Disable not implemented errors (here we are not storing any actual data).
InputVariant._save_data = lambda self: None


def reinit_from_file(tree_to_reinit):
    """Reinitialize tree from file."""
    if tree_to_reinit._treemanager.filename == "":
        return tree_to_reinit
    file = tree_to_reinit._treemanager.filename
    for leaf in (  tree_to_reinit._treemanager.inputs
                 + tree_to_reinit._treemanager.outputs):
        leaf._treemanager = None
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
    assert not len(tree._names)
    with pytest.raises(AttributeError):
        setattr(tree, "attribute", "value")

    for category in CATEGORIES:
        assert category in tree
        assert not len(tree[category]._names)
        with pytest.raises(AttributeError):
            setattr(tree[category], "attribute", "value")

    assert not tree._treemanager.filename


def test_init_emptyfile():
    """Test initializing a tree with empty file.

    Tree-node itself should not contain any QIDs, but all category nodes should
    be present both in the instance and in the file.
    """
    try:
        tree = Tree(CATEGORIES, (FNEMPTY, False))
        with TreeFile(FNEMPTY, "r") as h5:
            for group in CATEGORIES + ["results"]:
                assert group in h5
            for group in h5.keys():
                assert group in CATEGORIES + ["results"]

        assert not len(tree._names)
        for category in CATEGORIES:
            assert category in tree
            assert not len(tree[category]._names)
    finally:
        os.unlink(FNEMPTY)


def test_init_fromfile():
    """Test initializing a tree from a file containing data.

    The file contains single input used by single output. Test that the metadata
    is read correctly and the nodes have correct active fields.
    """
    try:
        with TreeFile(FNTEST, "a") as h5:
            h5.set_node(CATEGORIES[0], f"{INPUTS[0]}_1")
            h5.set_node(RESULTGROUP, f"{OUTPUTS[0]}_1")
            h5.set_datagroup(
                CATEGORIES[0], f"{INPUTS[0]}_1", date=DATE, note=NOTE,
                )
            h5.set_datagroup(
                RESULTGROUP, f"{OUTPUTS[0]}_1", date=DATE, note=NOTE,
                **{CATEGORIES[0]:f"{INPUTS[0]}_1"},
            )
            h5[f"{CATEGORIES[0]}/{INPUTS[0]}_1"]

        tree = Tree(CATEGORIES, (FNTEST, True))
        assert tree[CATEGORIES[0]]._names == [f"{INPUTS[0]}_1"]
        leafin = tree[CATEGORIES[0]][f"{INPUTS[0]}_1"]
        assert tree[CATEGORIES[0]].active is leafin
        assert leafin.date == DATE
        assert leafin.note == NOTE
        assert leafin.variant == INPUTS[0]

        assert tree._names == [f"{OUTPUTS[0]}_1"]
        leafout = tree[f"{OUTPUTS[0]}_1"]
        assert leafout[CATEGORIES[0]] is leafin
        assert tree.active is leafout
    finally:
        os.unlink(FNTEST)


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_add_input(tree):
    """Test adding a new input."""
    leaf1 = InputVariant()
    leaf1._cdata = object()

    tree._treemanager.enter_leaf(leaf1, CATEGORIES[0])
    tree = reinit_from_file(tree)
    assert leaf1.name in tree[CATEGORIES[0]], "Input was not added."

    leaf2 = InputVariant()
    leaf2._cdata = object()
    tree._treemanager.enter_leaf(leaf2, CATEGORIES[0])
    tree = reinit_from_file(tree)
    assert leaf1.name == tree[CATEGORIES[0]].active.name, (
        "Active input should not change unless explicitly activated.")

    leaf3 = InputVariant()
    leaf3._cdata = object()
    tree._treemanager.enter_leaf(leaf3, CATEGORIES[0], activate=True)
    tree = reinit_from_file(tree)
    assert leaf3.name == tree[CATEGORIES[0]].active.name, (
        "Active input was not changed even when explicitly activated.")


@pytest.mark.parametrize("tree", [True], indirect=True)
def test_add_input_and_save(tree):
    """Test adding new input without saving it, and then saving retroactively.
    """
    leaf1 = InputVariant()
    leaf1._cdata = object()
    tree._treemanager.enter_leaf(leaf1, CATEGORIES[0], save=False)
    tree = reinit_from_file(tree)
    assert not leaf1.name in tree[CATEGORIES[0]], (
        "Input was saved even though save=False")

    tree._treemanager.enter_leaf(leaf1, CATEGORIES[0], save=False)
    leaf1.save()
    tree = reinit_from_file(tree)
    assert leaf1.name in tree[CATEGORIES[0]], (
        "Input was not saved by its save method.")


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_add_output(tree):
    """Test adding a new output."""
    leaf = InputVariant()
    leaf._cdata = object()
    leaf1 = OutputVariant(inputs={CATEGORIES[0]: leaf})
    with pytest.raises(AscotMeltdownError):
        tree._treemanager.enter_leaf(leaf1)

    tree._treemanager.enter_leaf(leaf, CATEGORIES[0])

    tree._treemanager.enter_leaf(leaf1)
    tree = reinit_from_file(tree)
    assert leaf1.name in tree, "Output was not added."

    leaf2 = OutputVariant(inputs={CATEGORIES[0]: tree[CATEGORIES[0]].active})
    tree._treemanager.enter_leaf(leaf2)
    tree = reinit_from_file(tree)
    assert leaf1.name == tree.active.name, (
        "Active output should not change unless explicitly activated.")

    leaf3 = OutputVariant(inputs={CATEGORIES[0]: tree[CATEGORIES[0]].active})
    tree._treemanager.enter_leaf(leaf3, activate=True)
    tree = reinit_from_file(tree)
    assert leaf3.name == tree.active.name, (
        "Active input was not changed even when explicitly activated.")


@pytest.mark.parametrize("tree", [True], indirect=True)
def test_add_output_and_save(tree):
    """Test adding new output without saving it, and then saving retroactively.

    This test also checks that output cannot be saved unless input is saved
    first.
    """
    leafin = InputVariant()
    leafin._cdata = object()
    tree._treemanager.enter_leaf(leafin, CATEGORIES[0], save=False)

    leafout1 = OutputVariant(inputs={CATEGORIES[0]: leafin})
    with pytest.raises(AscotDataException):
        tree._treemanager.enter_leaf(leafout1)

    assert leafout1.name in tree, (
        "Output should be in tree even if save failed")

    leafin.save()
    leafout2 = OutputVariant(inputs={CATEGORIES[0]: tree[CATEGORIES[0]].active})
    tree._treemanager.enter_leaf(leafout2, save=False)
    tree = reinit_from_file(tree)
    assert not leafout2.name in tree, "Output was saved even though save=False"

    leafout2 = OutputVariant(inputs={CATEGORIES[0]: tree[CATEGORIES[0]].active})
    tree._treemanager.enter_leaf(leafout2, save=False)
    leafout2.save()
    tree = reinit_from_file(tree)
    assert leafout2.name in tree, "Output was not saved by its save method."


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_destroy(tree):
    """Test removing data."""
    leafin = InputVariant()
    leafin._cdata = object()
    tree._treemanager.enter_leaf(leafin, CATEGORIES[0])
    leafout = OutputVariant(inputs={CATEGORIES[0]: leafin})
    tree._treemanager.enter_leaf(leafout)

    with pytest.raises(AscotDataException):
        tree.destroy(data=leafin)

    tree = reinit_from_file(tree)
    assert leafin.name in tree[CATEGORIES[0]], "Used input was removed."

    tree.destroy(data=leafout.name)
    tree = reinit_from_file(tree)
    assert leafout not in tree

    tree.destroy(data=leafin.name)
    tree = reinit_from_file(tree)
    assert leafin not in tree[CATEGORIES[0]]


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_tree_activate_data(tree):
    """Test setting data active."""
    leaf1 = InputVariant()
    leaf2 = InputVariant()
    leaf1._cdata = object()
    leaf2._cdata = object()
    tree._treemanager.enter_leaf(leaf1, CATEGORIES[0])
    tree._treemanager.enter_leaf(leaf2, CATEGORIES[0])
    assert tree[CATEGORIES[0]].active == leaf1
    tree.activate(leaf2)
    tree = reinit_from_file(tree)
    assert tree[CATEGORIES[0]].active.name == leaf2.name


@pytest.mark.parametrize("tree", [False, True], indirect=True)
def test_tree_change_note(tree):
    """Test changing data's note."""
    leaf1 = InputVariant()
    leaf1._cdata = object()
    tree._treemanager.enter_leaf(leaf1, CATEGORIES[0])
    leaf1.note = "Get to the <choppa>"
    tree = reinit_from_file(tree)
    assert tree[CATEGORIES[0]]["CHOPPA"].name == leaf1.name

    leaf2 = OutputVariant(inputs={CATEGORIES[0]: tree[CATEGORIES[0]].active})
    tree._treemanager.enter_leaf(leaf2)
    leaf2.note = "Get to the <choppa>"
    tree = reinit_from_file(tree)
    assert tree["CHOPPA"].name == leaf2.name


@pytest.mark.parametrize("tree", [False], indirect=True)
def test_tree_contents(tree):
    """Test getting the contents in a string."""
    tree._treemanager.enter_leaf(
        InputVariant(date=DATES[0], note=NOTE),
        CATEGORIES[0]
        )
    tree._treemanager.enter_leaf(
        InputVariant(date=DATES[1], note=NOTE),
        CATEGORIES[0]
        )
    tree._treemanager.enter_leaf(
        InputVariant(date=DATES[2], note=NOTE),
        CATEGORIES[1]
        )
    leaf = tree._treemanager.get_leaf(f"{INPUTS[0]}_1")
    tree._treemanager.enter_leaf(
        OutputVariant(date=DATES[2], note=NOTE, inputs={CATEGORIES[1]: leaf})
        )
    leaf = tree._treemanager.get_leaf(f"{INPUTS[0]}_1")
    tree._treemanager.enter_leaf(
        OutputVariant(date=DATES[2], note=NOTE, inputs={CATEGORIES[1]: leaf})
        )

    expected = textwrap.dedent(
        """
        Inputs: [only active shown]
        catX      InputVariant_1  1997-08-29 02:14:00 + 1 other(s)
                       "Let off some steam <Bennett>"
        catY      InputVariant_3  1997-08-31 02:14:00 (no other inputs)
                       "Let off some steam <Bennett>"
        catZ      *no inputs*


        Simulations:
        OutputVariant_2 1997-08-31 02:14:00
                       "Let off some steam <Bennett>"
        OutputVariant_1 1997-08-31 02:14:00 [active]
                       "Let off some steam <Bennett>"
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
    leafin = InputVariant()
    leafin._cdata = object()
    with pytest.raises(AscotDataException) as e:
        tree._treemanager.enter_leaf(leafin, CATEGORIES[0], save=True)
    assert "no file has been specified" in str(e.value)
    assert tree._treemanager.get_leaf(f"{INPUTS[0]}_1").name == f"{INPUTS[0]}_1"

    leafout = OutputVariant(inputs={CATEGORIES[0]: leafin})
    with pytest.raises(AscotDataException):
        tree._treemanager.enter_leaf(leafout, save=True)
    assert "no file has been specified" in str(e.value)
    assert tree._treemanager.get_leaf(f"{OUTPUTS[0]}_1").name == f"{OUTPUTS[0]}_1"


@pytest.mark.parametrize("tree", [True], indirect=True)
def test_set_active_not_on_file(tree):
    """Test that we don't store active QID of data which is not on disk."""
    leaf1 = InputVariant()
    leaf1._cdata = object()
    tree._treemanager.enter_leaf(leaf1, CATEGORIES[0])
    leaf2 = InputVariant()
    leaf2._cdata = object()
    tree._treemanager.enter_leaf(leaf2, CATEGORIES[0], save=False)
    leaf2.activate()

    assert tree[CATEGORIES[0]].active.name == leaf2.name

    tree = reinit_from_file(tree)
    assert tree[CATEGORIES[0]].active.name == leaf1.name
