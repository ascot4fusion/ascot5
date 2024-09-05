# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests"""
from __future__ import annotations

import unittest
from unittest.mock import MagicMock

from a5py import AscotIOException

from a5py.data.access.metadata import input_categories, MetaData
from a5py.data.access.tree import InputCategory, RunVariant, Tree
from a5py.data.access.treeparts import Leaf

from .conftest import (
    QID1, QID2, QID3, DATE, NOTE, DATE_FRI, DATE_SAT, INPUTVAR, OUTPUTVAR,
    CATEGORY, create_leaf,
    )


class TestInputCategory(unittest.TestCase):
    """Tests for `InputCategory` class."""

    def test_contents(self):
        """Test that the contents of an input category are displayed
        correctly.
        """
        node = InputCategory()
        leaf1 = create_leaf(QID1, date=DATE_SAT)
        leaf2 = create_leaf(QID2, date=DATE_FRI)
        self.assertEqual(node.contents, "No data in this category.\n")

        node._add_leaf(leaf1)
        node._add_leaf(leaf2)
        lines = node.contents.splitlines()
        self.assertIn(leaf1.variant, lines[0])
        self.assertIn(leaf1.qid, lines[0])
        self.assertIn(leaf1.date, lines[0])
        self.assertIn("active", lines[0])
        self.assertIn(node._tags[0], lines[1])
        self.assertIn(leaf1.note, lines[2])
        self.assertIn(leaf2.variant, lines[4])
        self.assertIn(leaf2.qid, lines[4])
        self.assertIn(leaf2.date, lines[4])
        self.assertNotIn("active", lines[4])
        self.assertIn(node._tags[1], lines[5])
        self.assertIn(leaf2.note, lines[6])


class TestRunVariant(unittest.TestCase):
    """Tests for `RunVariant` class."""

    def setUp(self):
        """Initialize a single simulation output."""
        self.bfield = Leaf(
            qid=QID1,
            date=DATE,
            note=NOTE,
            variant=INPUTVAR,
        )
        self.wall = Leaf(
            qid=QID2,
            date=DATE,
            note=NOTE,
            variant=INPUTVAR,
        )
        self.inistate = MagicMock()
        self.node = RunVariant(
            inputs={"bfield":self.bfield, "wall":self.wall},
            diagnostics={"inistate":self.inistate},
            qid=QID3,
            date=DATE,
            note=NOTE,
            variant="run",
            )

    def test_reference_inputs(self):
        """Test that inputs are referenced correctly."""
        self.assertEqual(self.node[CATEGORY], self.wall)
        with self.assertRaises(AscotIOException):
            _ = self.node.efield

    def test_reference_diagnostics(self):
        """Test that diagnostics are referenced correctly."""
        self.assertEqual(self.node._inistate, self.inistate)
        with self.assertRaises(AscotIOException):
            _ = self.node._endstate

    def test_contents(self):
        """Test that the contents of an simulation output are displayed
        correctly.
        """
        lines = self.node.contents.splitlines()
        self.assertIn(self.node.variant, lines[0])
        self.assertIn(self.node.qid, lines[0])
        self.assertIn(self.node.date, lines[0])
        self.assertIn(self.node.note, lines[1])
        self.assertIn("Diagnostics", lines[3])
        self.assertIn("inistate", lines[4])
        self.assertIn("Inputs", lines[6])
        self.assertIn(CATEGORY, lines[9])
        self.assertIn(self.node[CATEGORY].variant, lines[9])
        self.assertIn(self.node[CATEGORY].qid, lines[9])
        self.assertIn(self.node[CATEGORY].date, lines[9])
        self.assertIn(self.node[CATEGORY].note, lines[10])


class TestTree(unittest.TestCase):
    """Tests for the `Tree` class (not including IO operations)."""

    def test_initialization(self):
        """Test initializing (empty) tree."""
        tree = Tree()
        self.assertFalse(len(tree._qids))
        with self.assertRaises(AscotIOException):
            setattr(tree, "attribute", "value")

        for category in input_categories:
            self.assertIn(category, tree)
            with self.assertRaises(AscotIOException):
                setattr(tree[category], "attribute", "value")

    def test_add_input(self):
        """Test adding a new input."""
        tree = Tree()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        leaf = tree._treemanager.enter_input(meta, dryrun=True)
        self.assertTrue(leaf)
        self.assertNotIn(leaf, tree[CATEGORY])
        leaf = tree._treemanager.enter_input(meta)
        self.assertTrue(leaf)
        self.assertIn(leaf, tree[CATEGORY])

    def test_add_simulation_output(self):
        """Test adding a new run."""
        tree = Tree()
        metaout = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        with self.assertRaises(
            AscotIOException,
            msg="Was able to add run whose input does not exist.",
            ):
            tree._treemanager.enter_run(metaout, [], [QID1])

        metain = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        tree._treemanager.enter_input(metain)
        leaf = tree._treemanager.enter_run(metaout, [], [QID1], dryrun=True)
        self.assertTrue(leaf)
        self.assertNotIn(leaf, tree)
        leaf = tree._treemanager.enter_run(metaout, [], [QID1])
        self.assertTrue(leaf)
        self.assertIn(leaf, tree)

    def test_add_identical_qid(self):
        """Test adding data when one with identical QID exists."""
        tree = Tree()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        tree._treemanager.enter_input(meta)
        with self.assertRaises(AscotIOException):
            tree._treemanager.enter_input(meta)

        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=OUTPUTVAR)
        with self.assertRaises(AscotIOException):
            tree._treemanager.enter_run(meta, [], [QID1])

        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        tree._treemanager.enter_run(meta, [], [QID1])
        with self.assertRaises(AscotIOException):
            tree._treemanager.enter_run(meta, [], [QID2])

    def data(self):
        """Test removing data."""
        tree = Tree()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        data = tree._treemanager.enter_input(meta)
        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        output = tree._treemanager.enter_run(meta, [], [QID1])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to destroy input that is being used in a simulation.",
            ):
            tree.destroy(data)

        tree.destroy(output)
        self.assertNotIn(output, tree)
        tree.destroy(data.qid)
        self.assertNotIn(data, tree[CATEGORY])

    def test_activate_data(self):
        """Test setting data active."""
        tree = Tree()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        tree._treemanager.enter_input(meta)
        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=INPUTVAR)
        data = tree._treemanager.enter_input(meta)
        tree.activate(data)
        self.assertEqual(tree[CATEGORY].active, data)

    def _test_contents(self):
        """Test getting the contents in a string."""
        tree = Tree()
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
        tree.show_contents()
        self.assertIn("Inputs: [only active shown]", lines[0])
        self.assertIn("*no inputs*", lines[1])
        self.assertIn("(no other inputs)", lines[5])
        self.assertIn(CATEGORY, lines[13])
        self.assertIn(data.variant, lines[13])
        self.assertIn(data.qid, lines[13])
        self.assertIn(data.date, lines[13])
        self.assertIn("+ 1 other(s)", lines[13])
        self.assertIn(data.note, lines[14])

        self.assertIn("Simulations:", lines[24])
        self.assertIn(output.variant, lines[25])
        self.assertIn(output.qid, lines[25])
        self.assertIn(output.date, lines[25])
        self.assertIn("[active]", lines[25])
        self.assertIn(output.note, lines[26])
        self.assertTrue(lines[27])
        self.assertNotIn("[active]", lines[27])


if __name__ == "__main__":
    unittest.main()
