# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests that cover the functionality of the data tree formed mainly by `Root`,
`ImmutableNode` and `MetaDataHolder` classes.

Note that this test handles the tree in its abstract form. All IO operations
are tested in separate modules.

If this test passes, there should not be any issues wit the `AscotIO` class when
the data is stored in the memory.
"""
from __future__ import annotations

import unittest
from unittest.mock import MagicMock

from a5py import AscotIOException

from a5py.ascot5io.coreio import metadata
from a5py.ascot5io.coreio.metadata import input_categories, MetaData
from a5py.ascot5io.coreio.treestructure import (
    Leaf, ImmutableNode, InputCategory, SimulationOutput, Root,
    )

QID1 = "2991786571"
QID2 = "9753987342"
QID3 = "4404229430"
QID4 = "0963810214"
QID5 = "5960585966"
DATE_FRI = "1997-08-29 02:14:00"
DATE_SAT = "1997-08-30 02:14:00"
DATE_SUN = "1997-08-31 02:14:00"
INPUTVAR = "input"
INPUTVAR2 = "E_TC"
OUTPUTVAR = "output"
CATEGORY = "wall"
DATE = DATE_FRI
NOTE = "Stick around"

# Make Ascot recognize our dummy variants
metadata.data_variants[CATEGORY] = (
    metadata.data_variants[CATEGORY] + (INPUTVAR,)
)
metadata.run_variants = metadata.run_variants + (OUTPUTVAR,)

def create_leaf(
        qid: str, date: str = DATE, variant: str = INPUTVAR, note: str = NOTE,
        ) -> Leaf:
    """Create a `Leaf` instance."""
    return Leaf(qid=qid, date=date, variant=variant, note=note)

class TestLeaf(unittest.TestCase):
    """Tests for `Leaf` class."""

    def test_initialization(self):
        """Test that attributes were set properly in initialization."""
        leaf = create_leaf(QID1)
        for actual, expected in zip(
            [leaf.qid, leaf.date, leaf.note, leaf.variant],
            [QID1, DATE, NOTE, INPUTVAR],
            ):
            self.assertEqual(actual, expected)
        self.assertEqual(
            leaf.qqid, f"q{QID1}", "No preceding 'q' in qqid.",
            )
        self.assertEqual(
            leaf.name, f"{INPUTVAR}_{QID1}",
            "Property name does not match <variety>_<QID>.",
            )

    def test_can_change_note(self):
        """Test that description can be changed."""
        leaf = create_leaf(QID1)
        leaf.note = "New note"
        self.assertEqual(leaf.note, "New note")

    def test_cannot_change_fixed_attributes(self):
        """Test that immutable attributes are immutable."""
        leaf = create_leaf(QID1)
        for attr in ["qid", "date", "variant", "qqid", "name"]:
            with self.assertRaises(
                AscotIOException,
                msg=f"User is able to change {attr}.",
                ):
                setattr(leaf, attr, QID2)

    def test_adopt_method(self):
        """Test that the fields are set correctly when the leaf is adopted and
        that the leaf cannot be adopted twice."""
        parent = MagicMock()
        leaf = create_leaf(QID1)
        leaf._adopt(parent)
        self.assertEqual(
            leaf._parent, parent,
            "Parent was not set when the leaf was adopted.",
            )
        with self.assertRaises(AscotIOException,
                msg="Was able to adopt same leaf twice.",
            ):
            leaf._adopt(parent)

    def test_extract_tag(self):
        """Test extracting a tag from description."""
        leaf = create_leaf(QID1)
        leaf.note = "tag me"
        self.assertEqual(
            leaf._extract_tag(), "TAG",
            "Extracted tag is not the first word in note capitalized.",
            )
        leaf.note = "12345note with numbers"
        self.assertEqual(
            leaf._extract_tag(), "TAG",
            "Starting note with numbers should result in default tag.",
            )
        leaf.note = "!@#$%^&*()special_characters"
        self.assertEqual(
            leaf._extract_tag(), "SPECIALCHARACTERS",
            "Tag creation should ignore special characters.",
            )

    def test_activate_and_remove_methods(self):
        """Test that activate and remove methods convey the call to the parent
        node.
        """
        leaf = create_leaf(QID1)
        parent = MagicMock()
        with parent._modify_attributes():
            parent._activate_leaf = MagicMock()
            parent._remove_leaf = MagicMock()

        leaf._adopt(parent)
        leaf.activate()
        parent._activate_leaf.assert_called_with(leaf)
        leaf.destroy()
        parent._remove_leaf.assert_called_with(leaf)


class TestImmutableNode(unittest.TestCase):
    """Tests for `ImmutableNode` class.

    There's no separate test for the `ImmutableStorage` class as it is also
    tested here implicitly.
    """

    def setUp(self):
        """Set up a node and mock root."""
        self.root = MagicMock()
        self.node = ImmutableNode(root=self.root)

    def test_initialization(self):
        """Test that the initialized object us unfrozen and is empty."""
        self.assertFalse(self.node._frozen)
        self.assertEqual(self.node._qids, [])
        self.assertEqual(self.node._tags, [])
        self.assertEqual(self.node._root, self.root)
        self.assertIsNone(self.node._active)

    def test_freeze_and_unfreeze(self):
        """Test that the attributes can only be modified when the instance is
        unfrozen.
        """
        self.node._freeze()
        self.assertTrue(self.node._frozen)
        with self.assertRaises(
            AscotIOException,
            msg="Was able to modify node while frozen.",
            ):
            self.node.new_attr = "value"

        self.node._unfreeze()
        self.assertFalse(self.node._frozen,)
        self.node.new_attr = "value"
        self.assertEqual(self.node.new_attr, "value")

    def test_modify_attributes_context(self):
        """Test the modify attributes context."""
        self.node._freeze()
        with self.node._modify_attributes():
            self.node.new_attr = "value"
        self.assertEqual(
            self.node.new_attr, "value",
            "Failed to modify node within _modify_attributes context."
            )

    def test_setitem_and_getitem(self):
        """Test the dictionary-like attribute access."""
        self.node["attr"] = "value"
        self.assertEqual(
            self.node.attr, "value",
            "Failed to set attribute in dictionary-like syntax.",
            )
        self.assertEqual(
            self.node["attr"], "value",
            "Failed to set retrieve attribute using dictionary-like syntax.",
            )
        self.node._freeze()
        with self.assertRaises(
            AscotIOException,
            msg="Able to modify frozen node using dictionary-like syntax.",
            ):
            self.node["attr"] = "another_value"

    def test_add_leaf(self):
        """Test adding a single leaf."""
        leaf = create_leaf(QID1)
        self.node._add_leaf(leaf)
        self.assertEqual(leaf, self.node[leaf.qqid])
        self.assertEqual(leaf, self.node[leaf.name])
        self.assertEqual(leaf, self.node[leaf._extract_tag()])
        self.assertEqual(self.node._qids, [QID1])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to add same leaf twice.",
            ):
            self.node._add_leaf(leaf)

    def test_contains(self):
        """Test that the contained leaf can be queried with QID, name, and with
        the leaf itself"""
        leaf = create_leaf(QID1)
        self.node._add_leaf(leaf)
        self.assertIn(leaf.qid, self.node,)
        self.assertIn(leaf.qqid, self.node,)
        self.assertIn(leaf.name, self.node)
        self.assertIn(leaf, self.node)

    def test_remove_leaf(self):
        """Test removing a single leaf."""
        leaf = create_leaf(QID1)
        self.node._add_leaf(leaf)
        self.node._remove_leaf(leaf)
        self.assertNotIn(leaf.qid, self.node)
        self.assertNotIn(leaf.name, self.node)
        self.assertNotIn(leaf._extract_tag(), self.node)
        self.assertNotIn(leaf, self.node)
        self.assertEqual(self.node._qids, [])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to remove same leaf twice.",
            ):
            self.node._remove_leaf(leaf)

    def test_active_property(self):
        """Test switching active leaf."""
        leaf1 = create_leaf(QID1)
        leaf2 = create_leaf(QID2)
        self.node._add_leaf(leaf1)
        self.node._add_leaf(leaf2)
        self.assertEqual(self.node.active, leaf1)
        self.node._remove_leaf(leaf1)
        self.assertEqual(self.node.active, leaf2)
        self.node._remove_leaf(leaf2)
        with self.assertRaises(AscotIOException):
            _ = self.node.active

    def test_organize_by_date(self):
        """Test that leafs remain organized by date when new ones are added or
        removed.
        """
        leaf1 = create_leaf(QID1, date=DATE_SAT, note="Saturday")
        leaf2 = create_leaf(QID2, date=DATE_FRI, note="Friday")
        leaf3 = create_leaf(QID3, date=DATE_SUN, note="Sunday")
        self.node._add_leaf(leaf1)
        self.assertEqual(self.node._qids, [QID1])
        self.node._add_leaf(leaf2)
        self.assertEqual(self.node._qids, [QID1, QID2])
        self.node._add_leaf(leaf3)
        self.assertEqual(self.node._qids, [QID3, QID1, QID2])
        self.assertEqual(self.node._tags, ["FRIDAY", "SATURDAY", "SUNDAY"])
        self.node._remove_leaf(leaf3)
        self.assertEqual(self.node._qids, [QID1, QID2])
        self.assertEqual(self.node._tags, ["FRIDAY", "SATURDAY"])
        self.node._remove_leaf(leaf2)
        self.assertEqual(self.node._qids, [QID1])

    def test_organize_same_tags(self):
        """Test that tags are properly renamed when there are multiple identical
        tags"""
        leaf1 = create_leaf(QID1, date=DATE_SAT, note="TAG")
        leaf2 = create_leaf(QID2, date=DATE_FRI, note="TAG2")
        leaf3 = create_leaf(QID3, date=DATE_SUN, note="TAG")
        self.node._add_leaf(leaf1)
        self.node._add_leaf(leaf2)
        self.assertIn("TAG", self.node)
        self.assertIn("TAG2", self.node)
        self.node._add_leaf(leaf3)
        self.assertNotIn("TAG", self.node)
        self.assertIn("TAG_0", self.node)
        self.assertIn("TAG_1", self.node)
        self.assertIn("TAG2", self.node)

        self.assertEqual(self.node.TAG_0, leaf3)
        self.assertEqual(self.node.TAG_1, leaf1)
        self.assertEqual(self.node.TAG2, leaf2)

        self.node._remove_leaf(leaf3)
        self.assertIn("TAG", self.node)
        self.assertNotIn("TAG_0", self.node)
        self.assertNotIn("TAG_1", self.node)
        self.assertEqual(self.node.TAG, leaf1)

    def test_update_tag(self):
        """Test that the tag is updated when the note is changed."""
        leaf = create_leaf(QID1, note="tag")
        self.node._add_leaf(leaf)
        self.assertIn("TAG", self.node)

        leaf.note = "Newtag"
        self.assertIn("NEWTAG", self.node)


class TestInputCategory(unittest.TestCase):
    """Tests for `InputCategory` class."""

    def test_contents(self):
        """Test that the contents of an input category are displayed
        correctly.
        """
        root = MagicMock()
        node = InputCategory(root)
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


class TestSimulationOutput(unittest.TestCase):
    """Tests for `SimulationOutput` class."""

    def setUp(self):
        """Initialize a single simulation output."""
        self.bfield = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            note="TAG",
            variant="test",
        )
        self.wall = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            note="TAG",
            variant="wall",
        )
        self.inistate = MagicMock()
        self.node = SimulationOutput(
            inputs={"bfield":self.bfield, "wall":self.wall},
            diagnostics={"inistate":self.inistate},
            qid="0000000002",
            date="1953-12-08 00:00:01",
            note="TAG",
            variant="run",
            )

    def test_reference_inputs(self):
        """Test that inputs are referenced correctly."""
        self.assertEqual(self.node.bfield, self.bfield)
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
        self.assertIn("bfield", lines[7])
        self.assertIn(self.node.bfield.variant, lines[7])
        self.assertIn(self.node.bfield.qid, lines[7])
        self.assertIn(self.node.bfield.date, lines[7])
        self.assertIn(self.node.bfield.note, lines[8])


class TestRoot(unittest.TestCase):
    """Tests for the `Root` class (not including IO operations)."""

    def test_initialization(self):
        """Test initializing (empty) tree."""
        root = Root()
        self.assertFalse(len(root._qids))
        self.assertEqual(root._name, "root")
        with self.assertRaises(AscotIOException):
            setattr(root, "attribute", "value")

        for category in input_categories:
            self.assertIn(category, root)
            self.assertEqual(root[category]._name, category)
            with self.assertRaises(AscotIOException):
                setattr(root[category], "attribute", "value")

    def test_add_input(self):
        """Test adding a new input."""
        root = Root()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        leaf = root._add_input_dataset(meta)
        self.assertIn(leaf, root[CATEGORY])
        self.assertEqual(leaf._parent, root[CATEGORY])

    def test_add_simulation_output(self):
        """Test adding a new run."""
        root = Root()
        metaout = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        with self.assertRaises(
            AscotIOException,
            msg="Was able to add run whose input does not exist.",
            ):
            root._add_simulation_output(metaout, [], [QID1])

        metain = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        root._add_input_dataset(metain)
        output = root._add_simulation_output(metaout, [], [QID1])
        self.assertIn(output, root)

    def test_add_identical_qid(self):
        """Test adding dataset when one with identical QID exists."""
        root = Root()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        root._add_input_dataset(meta)
        with self.assertRaises(AscotIOException):
            root._add_input_dataset(meta)

        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=OUTPUTVAR)
        with self.assertRaises(AscotIOException):
            root._add_simulation_output(meta, [], [QID1])

        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        root._add_simulation_output(meta, [], [QID1])
        with self.assertRaises(AscotIOException):
            root._add_simulation_output(meta, [], [QID2])

    def test_remove_dataset(self):
        """Test removing dataset."""
        root = Root()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        dataset = root._add_input_dataset(meta)
        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        output = root._add_simulation_output(meta, [], [QID1])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to destroy input that is being used in a simulation.",
            ):
            root.destroy_dataset(dataset.qid)

        root.destroy_dataset(output)
        self.assertNotIn(output, root)
        root.destroy_dataset(dataset.qid)
        self.assertNotIn(dataset, root[CATEGORY])

    def test_activate_dataset(self):
        """Test setting dataset active."""
        root = Root()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        root._add_input_dataset(meta)
        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=INPUTVAR)
        dataset = root._add_input_dataset(meta)
        root.activate_dataset(dataset)
        self.assertEqual(root[CATEGORY].active, dataset)

    def test_contents(self):
        """Test getting the contents in a string."""
        root = Root()
        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        dataset = root._add_input_dataset(meta)
        meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=INPUTVAR)
        root._add_input_dataset(meta)
        meta = MetaData(qid=QID4, date=DATE, note=NOTE, variant=INPUTVAR2)
        root._add_input_dataset(meta)
        meta = MetaData(qid=QID2, date=DATE, note=NOTE, variant=OUTPUTVAR)
        output = root._add_simulation_output(meta, [], [QID1])
        meta = MetaData(qid=QID5, date=DATE, note=NOTE, variant=OUTPUTVAR)
        root._add_simulation_output(meta, [], [QID1])

        lines = root.contents.splitlines()
        root.show_contents()
        self.assertIn("Inputs: [only active shown]", lines[0])
        self.assertIn("*no inputs*", lines[1])
        self.assertIn("(no other inputs)", lines[5])
        self.assertIn(CATEGORY, lines[13])
        self.assertIn(dataset.variant, lines[13])
        self.assertIn(dataset.qid, lines[13])
        self.assertIn(dataset.date, lines[13])
        self.assertIn("+ 1 other(s)", lines[13])
        self.assertIn(dataset.note, lines[14])

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
