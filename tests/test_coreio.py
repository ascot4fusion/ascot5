# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests that cover the functionality of the data tree formed mainly by `Root`,
`ImmutableNode` and `MetaDataHolder` classes.

Note that this test handles the tree in its abstract form. All IO operations
are tested in separate modules.

If this test passes, there should not be any issues wit the `AscotIO` class when
the data is stored in the memory.
"""
import unittest
from unittest.mock import MagicMock

from a5py import AscotIOException
from a5py.ascot5io.coreio.metadata import input_categories, MetaData
from a5py.ascot5io.coreio.treestructure import (
    Leaf, ImmutableNode, InputCategory, SimulationOutput, Root,
    )


class TreeTester(unittest.TestCase):
    """Helper class that implements methods for testing the tree structure."""

    QID1 = "0000000001"
    QID2 = "0000000002"
    QID3 = "0000000003"
    QID4 = "0000000004"
    QID5 = "0000000005"
    DATETODAY = "1953-12-08 00:00:00"
    DATETOMORROW = "1953-12-09 00:00:00"
    DATEYESTERDAY = "1953-12-07 00:00:00"
    DEFAULTINPUT = "input"
    DEFAULTOUTPUT = "output"
    DEFAULTNOTE = "TAG"

    def create_leaf(
            self, qid: str, date: str = DATETODAY,
            variant: str = DEFAULTINPUT, note: str = DEFAULTNOTE,
            ) -> Leaf:
        """Create a `Leaf` instance.
        """
        return Leaf(qid=qid, date=date, variant=variant, note=note)

    def assertLeaf(self, obj, **expected_attrs):
        """Assert that the `Leaf` instance has expected values for the
        attributes.

        (This method does not appear in the error stack.)
        """
        for attr, expected in expected_attrs.items():
            actual = getattr(obj, attr, None)
            try:
                self.assertEqual(
                    actual, expected,
                    f"Expected {attr} to be {expected}, but got {actual}."
                    )
            except AssertionError as e:
                err = AssertionError(e).with_traceback(e.__traceback__.tb_next)
                raise err from None

class TestLeaf(TreeTester):
    """Tests for `Leaf` class."""

    def test_initialization(self):
        """Test that attributes were set properly in initialization."""
        leaf = self.create_leaf(TreeTester.QID1)
        self.assertLeaf(
            leaf, qid=TreeTester.QID1, date=TreeTester.DATETODAY,
            note=TreeTester.DEFAULTNOTE, variant=TreeTester.DEFAULTINPUT,
        )
        self.assertEqual(
            leaf.qqid, f"q{TreeTester.QID1}", "No preceding 'q' in qqid.",
            )
        self.assertEqual(
            leaf.name, f"{TreeTester.DEFAULTINPUT}_{TreeTester.QID1}",
            "Property name does not match <variety>_<QID>.",
            )

    def test_can_change_note(self):
        """Test that description can be changed."""
        leaf = self.create_leaf(TreeTester.QID1)
        leaf.note = "New note"
        self.assertEqual(leaf.note, "New note")

    def test_cannot_change_fixed_attributes(self):
        """Test that immutable attributes are immutable."""
        leaf = self.create_leaf(TreeTester.QID1)
        for attr in ["qid", "date", "variant", "qqid", "name"]:
            with self.assertRaises(
                AscotIOException,
                msg=f"User is able to change {attr}.",
                ):
                setattr(leaf, attr, TreeTester.QID2)

    def test_adopt_method(self):
        """Test that the fields are set correctly when the leaf is adopted and
        that the leaf cannot be adopted twice."""
        parent = MagicMock()
        leaf = self.create_leaf(TreeTester.QID1)
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
        leaf = self.create_leaf(TreeTester.QID1)
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
        leaf = self.create_leaf(TreeTester.QID1)
        parent = MagicMock()
        with parent._modify_attributes():
            parent._activate_leaf = MagicMock()
            parent._remove_leaf = MagicMock()

        leaf._adopt(parent)
        leaf.activate()
        parent._activate_leaf.assert_called_with(leaf)
        leaf.destroy()
        parent._remove_leaf.assert_called_with(leaf)


class TestImmutableNode(TreeTester):
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
        leaf = self.create_leaf(TreeTester.QID1)
        self.node._add_leaf(leaf)
        self.assertEqual(leaf, self.node[leaf.qqid])
        self.assertEqual(leaf, self.node[leaf.name])
        self.assertEqual(leaf, self.node[leaf._extract_tag()])
        self.assertEqual(self.node._qids, [TreeTester.QID1])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to add same leaf twice.",
            ):
            self.node._add_leaf(leaf)

    def test_contains(self):
        """Test that the contained leaf can be queried with QID, name, and with
        the leaf itself"""
        leaf = self.create_leaf(TreeTester.QID1)
        self.node._add_leaf(leaf)
        self.assertIn(leaf.qid, self.node,)
        self.assertIn(leaf.qqid, self.node,)
        self.assertIn(leaf.name, self.node)
        self.assertIn(leaf, self.node)

    def test_remove_leaf(self):
        """Test removing a single leaf."""
        leaf = self.create_leaf(TreeTester.QID1)
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
        leaf1 = self.create_leaf(TreeTester.QID1)
        leaf2 = self.create_leaf(TreeTester.QID2)
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
        leaf1 = self.create_leaf(TreeTester.QID1, date=TreeTester.DATETODAY,
                                 note="TODAY")
        leaf2 = self.create_leaf(TreeTester.QID2, date=TreeTester.DATEYESTERDAY,
                                 note="YESTERDAY")
        leaf3 = self.create_leaf(TreeTester.QID3, date=TreeTester.DATETOMORROW,
                                 note="TOMORROW")
        self.node._add_leaf(leaf1)
        self.assertEqual(self.node._qids, [TreeTester.QID1])
        self.node._add_leaf(leaf2)
        self.assertEqual(self.node._qids, [TreeTester.QID1, TreeTester.QID2])
        self.node._add_leaf(leaf3)
        self.assertEqual(self.node._qids,
                         [TreeTester.QID3, TreeTester.QID1, TreeTester.QID2])
        self.assertEqual(self.node._tags, ["TODAY", "TOMORROW", "YESTERDAY"])
        self.node._remove_leaf(leaf3)
        self.assertEqual(self.node._qids, [TreeTester.QID1, TreeTester.QID2])
        self.assertEqual(self.node._tags, ["TODAY", "YESTERDAY"])
        self.node._remove_leaf(leaf2)
        self.assertEqual(self.node._qids, [TreeTester.QID1])

    def test_organize_same_tags(self):
        """Test that tags are properly renamed when there are multiple identical
        tags"""
        leaf1 = self.create_leaf(TreeTester.QID1, date=TreeTester.DATETODAY,
                                 note="TAG")
        leaf2 = self.create_leaf(TreeTester.QID2, date=TreeTester.DATEYESTERDAY,
                                 note="TAG2")
        leaf3 = self.create_leaf(TreeTester.QID3, date=TreeTester.DATETOMORROW,
                                 note="TAG")
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
        leaf = self.create_leaf(TreeTester.QID1)
        self.node._add_leaf(leaf)
        self.assertIn("TAG", self.node)

        leaf.note = "Newtag"
        self.assertIn("NEWTAG", self.node)


class TestInputCategory(TreeTester):
    """Tests for `InputCategory` class."""

    def test_contents(self):
        """Test that the contents of an input category are displayed
        correctly.
        """
        root = MagicMock()
        node = InputCategory(root)
        leaf1 = self.create_leaf(TreeTester.QID1, date=TreeTester.DATETODAY)
        leaf2 = self.create_leaf(TreeTester.QID2, date=TreeTester.DATEYESTERDAY)
        self.assertEqual(
            node.contents,
            "No data in this category.\n"
            )

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
    """Tests for the `Root` class.

    These tests don't cover any IO operations.
    """

    def setUp(self):
        """Initialize an empty tree."""
        self.root = Root()
        self.inistate = MagicMock()
        self.metadata = {}
        self.metadata["efield"] = (
            MetaData("0000000001", "1953-12-08 00:00:00", "note", "E_TC")
        )
        self.metadata["bfield"] = (
            MetaData("0000000002", "1953-12-08 00:00:00", "note", "B_TC")
        )
        self.metadata["bfield2"] = (
            MetaData("0000000003", "1953-12-08 00:00:00", "note", "B_TC")
        )
        self.metadata["bfield_identical_qid"] = (
            MetaData("0000000002", "1953-12-08 00:00:00", "note", "B_TC")
        )
        self.metadata["output"] = (
            MetaData("0000000004", "1953-12-08 00:00:01", "note", "run")
        )
        self.metadata["output_identical_qid"] = (
            MetaData("0000000002", "1953-12-08 00:00:01", "note", "run")
        )

    def test_initialization(self):
        """Test that the tree is initially empty with only input categories
        present.
        """
        self.assertFalse(len(self.root._qids))
        with self.assertRaises(AscotIOException):
            setattr(self.root, "attribute", 0)

        for category in input_categories:
            self.assertIn(category, self.root)
            with self.assertRaises(AscotIOException):
                setattr(self.root[category], "attribute", 0)

    def test_add_input(self):
        """Test adding a new input."""
        leaf = self.root._add_input_dataset(
            self.metadata["bfield"], "note",
            )
        self.assertIn(leaf, self.root.bfield)
        self.assertEqual(leaf._parent, self.root.bfield)

    def test_add_simulation_output(self):
        """Test adding a new run."""
        with self.assertRaises(AscotIOException):
            self.root._add_simulation_output(
                self.metadata["output"], [self.inistate],
                [self.metadata["bfield"].qid], "note",
                )

        self.root._add_input_dataset(self.metadata["bfield"])
        output = self.root._add_simulation_output(
            self.metadata["output"], [self.inistate],
            [self.metadata["bfield"].qid], "note",
            )
        self.assertIn(output, self.root)

    def test_add_identical_qid(self):
        """Test adding an input or output when there is data with identical
        QID.
        """
        self.root._add_input_dataset(self.metadata["bfield"])
        with self.assertRaises(AscotIOException):
            self.root._add_input_dataset(self.metadata["bfield_identical_qid"])
        with self.assertRaises(AscotIOException):
            self.root._add_simulation_output(
                self.metadata["output_identical_qid"], [self.inistate],
                [self.metadata["bfield"].qid], "note",
            )

    def test_remove_dataset(self):
        """Test removing dataset."""
        bfield = self.root._add_input_dataset(self.metadata["bfield"])
        output = self.root._add_simulation_output(
                self.metadata["output"], [self.inistate],
                [self.metadata["bfield"].qid], "note",
            )
        self.root.destroy_dataset(output)
        self.assertNotIn(output, self.root)
        self.root.destroy_dataset(bfield.qid)
        self.assertNotIn(bfield, self.root.bfield)

    def test_remove_dataset_dependent(self):
        """Test removing input which is being used in a run and contents of
        whole nodes at once."""
        bfield = self.root._add_input_dataset(self.metadata["bfield"])
        output = self.root._add_simulation_output(
                self.metadata["output"], [self.inistate],
                [self.metadata["bfield"].qid], "note",
            )
        self.root.destroy_dataset(bfield.qid)
        self.root.destroy_dataset(output)

    def test_activate_dataset(self):
        """Test setting dataset active."""
        self.root._add_input_dataset(self.metadata["bfield"])
        bfield2 = self.root._add_input_dataset(self.metadata["bfield2"])
        self.root.activate_dataset(bfield2)
        self.assertEqual(self.root.bfield.active, bfield2)

    def test_contents(self):
        """Test getting the contents in a string."""
        bfield = self.root._add_input_dataset(self.metadata["bfield"])
        self.root._add_input_dataset(self.metadata["bfield2"])
        self.root._add_input_dataset(self.metadata["efield"])
        output = self.root._add_simulation_output(
                self.metadata["output"], [self.inistate],
                [self.metadata["bfield"].qid], "note",
            )
        lines = self.root.contents.splitlines()
        self.assertIn("Inputs: [only active shown]", lines[0])
        self.assertIn("bfield", lines[3])
        self.assertIn(bfield.variant, lines[3])
        self.assertIn(bfield.qid, lines[3])
        self.assertIn(bfield.date, lines[3])
        self.assertIn("+ 1 other(s)", lines[3])
        self.assertIn(bfield.note, lines[4])
        self.assertIn("(no other inputs)", lines[5])
        self.assertIn("*no inputs*", lines[7])


if __name__ == "__main__":
    unittest.main()
