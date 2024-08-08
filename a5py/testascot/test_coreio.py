# pylint: disable=protected-access, no-member, too-many-instance-attributes
"""Tests that cover the functionality of the data tree formed mainly by `Root`,
`ImmutableNode` and `MetaDataHolder` classes.

Note that this test handles the tree in its abstract form. All IO operations
are tested in separate modules.
"""
import unittest
from unittest.mock import MagicMock

from a5py import AscotIOException
from a5py.ascot5io.coreio.treestructure import (
    Leaf, ImmutableNode, InputCategory, SimulationOutput, Root,
    input_categories,
    )

class TestLeaf(unittest.TestCase):
    """Tests for `Leaf` class."""

    def setUp(self):
        """Initialize single leaf for testing."""
        self.leaf = Leaf(
            qid="0123456789",
            date="1953-12-08 00:00:00",
            description="Test description",
            variant="test",
            )

    def test_initialization(self):
        """Test that attributes were set properly in initialization."""
        self.assertEqual(self.leaf.qid, "0123456789")
        self.assertEqual(self.leaf.date, "1953-12-08 00:00:00")
        self.assertEqual(self.leaf.description, "Test description")
        self.assertEqual(self.leaf.variant, "test")

    def test_description_mutability(self):
        """Test that description can be changed."""
        self.leaf.description = "New description"
        self.assertEqual(
            self.leaf.description, "New description",
            "Failed to set new description.",
            )

    def test_attribute_immutability(self):
        """Test that immutable attributes are immutable."""
        with self.assertRaises(
            AscotIOException,
            msg="User is able to change QID.",
            ):
            self.leaf.qid = "1234567890"
        with self.assertRaises(
            AscotIOException,
            msg="User is able to change date.",
            ):
            self.leaf.date = "1953-12-08 00:00:01"
        with self.assertRaises(
            AscotIOException,
            msg="User is able to change variant.",
            ):
            self.leaf.variant = "newvariant"
        with self.assertRaises(
            AscotIOException,
            msg="User is able to change qQID.",
            ):
            self.leaf.qqid = "q1234567890"
        with self.assertRaises(
            AscotIOException,
            msg="User is able to change name.",
            ):
            self.leaf.name = "1234567890_name"

    def test_adopt_method(self):
        """Test that the fields are set correctly when the leaf is adopted and
        that the leaf cannot be adopted twice."""
        parent = MagicMock()
        self.leaf._adopt(parent)
        self.assertEqual(
            self.leaf._parent, parent,
            "Parent was not set when the leaf was adopted.",
            )
        with self.assertRaises(
            AscotIOException,
            msg="Was able to adopt same leaf twice.",
            ):
            self.leaf._adopt(parent)

    def test_extract_tag(self):
        """Test extracting a tag from description."""
        self.assertEqual(
            self.leaf._extract_tag(), "TEST",
            "Extracted tag is not the first word in description capitalized.",
            )
        self.leaf.description = "12345description with numbers"
        self.assertEqual(
            self.leaf._extract_tag(), "TAG",
            "Starting description with numbers should result in default tag.",
            )
        self.leaf.description = "!@#$%^&*()special_characters"
        self.assertEqual(
            self.leaf._extract_tag(), "SPECIALCHARACTERS",
            "Tag creation should ignore special characters.",
            )

    def test_properties(self):
        """Test derived properties."""
        self.assertEqual(
            self.leaf.qqid, "q0123456789",
            "Property qqid did not add the q prefix."
            )
        self.assertEqual(
            self.leaf.name, "test_0123456789",
            "Property name does not match <variety>_<QID>."
            )

    def test_activate_and_remove_methods(self):
        """Test that activate and remove methods convey the call to the parent
        node.
        """
        parent = Root()
        with parent._modify_attributes():
            parent.activate_dataset = MagicMock()
            parent._remove_leaf = MagicMock()

        self.leaf._adopt(parent)
        self.leaf.activate()
        parent.activate_dataset.assert_called_with("0123456789")
        self.leaf.destroy()
        parent._remove_leaf.assert_called_with(self.leaf)


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

    def test_contains(self):
        """Test that the contained leaf can be queried with QID, name, and with
        the leaf itself"""
        leaf = Leaf(
            qid="0123456789",
            date="1953-12-08 00:00:00",
            description="Test description",
            variant="test",
            )
        self.node._add_leaf(leaf)
        self.assertIn("0123456789", self.node,)
        self.assertIn("q0123456789", self.node,)
        self.assertIn("test_0123456789", self.node)
        self.assertIn(leaf, self.node)

    def test_add_leaf(self):
        """Test adding a single leaf."""
        leaf = Leaf(
            qid="0123456789",
            date="1953-12-08 00:00:00",
            description="Test description",
            variant="test",
            )
        self.node._add_leaf(leaf)
        self.assertEqual(leaf, self.node.q0123456789)
        self.assertEqual(leaf, self.node.test_0123456789)
        self.assertEqual(leaf, self.node.TEST)
        self.assertEqual(self.node._qids, ["0123456789"])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to add same leaf twice.",
            ):
            self.node._add_leaf(leaf)

    def test_remove_leaf(self):
        """Test removing a single leaf."""
        leaf = Leaf(
            qid="0123456789",
            date="1953-12-08 00:00:00",
            description="Test description",
            variant="test",
            )
        self.node._add_leaf(leaf)
        self.node._remove_leaf(leaf)
        self.assertNotIn("q0123456789", self.node)
        self.assertNotIn("leaf_0123456789", self.node)
        self.assertNotIn("TEST", self.node)
        self.assertNotIn(leaf, self.node)
        self.assertEqual(self.node._qids, [])

        with self.assertRaises(
            AscotIOException,
            msg="Was able to remove same leaf twice.",
            ):
            self.node._remove_leaf(leaf)

    def test_active_property(self):
        """Test switching active leaf."""
        leaf1 = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="Test description",
            variant="test",
            )
        leaf2 = Leaf(
            qid="0000000002",
            date="1953-12-08 00:00:01",
            description="Test description",
            variant="test",
            )
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
        leaf1 = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="Today",
            variant="test",
            )
        leaf2 = Leaf(
            qid="0000000002",
            date="1953-12-07 00:00:00",
            description="Yesterday",
            variant="test",
            )
        leaf3 = Leaf(
            qid="0000000003",
            date="1953-12-09 00:00:00",
            description="Tomorrow",
            variant="test",
            )
        self.node._add_leaf(leaf1)
        self.assertEqual(self.node._qids, ["0000000001"])
        self.node._add_leaf(leaf2)
        self.assertEqual(self.node._qids, ["0000000001", "0000000002"])
        self.node._add_leaf(leaf3)
        self.assertEqual(self.node._qids,
                         ["0000000003", "0000000001", "0000000002"])
        self.assertEqual(self.node._tags, ["TODAY", "TOMORROW", "YESTERDAY"])
        self.node._remove_leaf(leaf3)
        self.assertEqual(self.node._qids, ["0000000001", "0000000002"])
        self.node._remove_leaf(leaf2)
        self.assertEqual(self.node._qids, ["0000000001"])

    def test_organize_same_tags(self):
        """Test that tags are properly renamed when there are multiple identical
        tags"""
        leaf1 = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="TAG Today",
            variant="test",
            )
        leaf2 = Leaf(
            qid="0000000002",
            date="1953-12-07 00:00:00",
            description="TAG2 Yesterday",
            variant="test",
            )
        leaf3 = Leaf(
            qid="0000000003",
            date="1953-12-09 00:00:00",
            description="TAG Tomorrow",
            variant="test",
            )
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


class TestInputCategory(unittest.TestCase):
    """Tests for `InputCategory` class."""

    def test_contents(self):
        """Test that the contents of an input category are displayed
        correctly.
        """
        root = MagicMock()
        node = InputCategory(root)
        leaf1 = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="Today",
            variant="test",
            )
        leaf2 = Leaf(
            qid="0000000002",
            date="1953-12-07 00:00:00",
            description="Yesterday",
            variant="test",
            )
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
        self.assertIn(leaf1.description, lines[2])
        self.assertIn(leaf2.variant, lines[4])
        self.assertIn(leaf2.qid, lines[4])
        self.assertIn(leaf2.date, lines[4])
        self.assertNotIn("active", lines[4])
        self.assertIn(node._tags[1], lines[5])
        self.assertIn(leaf2.description, lines[6])


class TestSimulationOutput(unittest.TestCase):
    """Tests for `SimulationOutput` class."""

    def setUp(self):
        """Initialize a single simulation output."""
        self.bfield = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="TAG",
            variant="test",
        )
        self.wall = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="TAG",
            variant="wall",
        )
        self.inistate = MagicMock()
        self.node = SimulationOutput(
            inputs={"bfield":self.bfield, "wall":self.wall},
            diagnostics={"inistate":self.inistate},
            qid="0000000002",
            date="1953-12-08 00:00:01",
            description="TAG",
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
        self.assertIn(self.node.description, lines[1])
        self.assertIn("Diagnostics", lines[3])
        self.assertIn("inistate", lines[4])
        self.assertIn("Inputs", lines[6])
        self.assertIn("bfield", lines[7])
        self.assertIn(self.node.bfield.variant, lines[7])
        self.assertIn(self.node.bfield.qid, lines[7])
        self.assertIn(self.node.bfield.date, lines[7])
        self.assertIn(self.node.bfield.description, lines[8])

class TestRoot(unittest.TestCase):
    """Tests for the `Root` class.

    These tests don't cover any IO operations.
    """

    def setUp(self):
        """Initialize an empty tree."""
        self.root = Root()
        self.efield = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="Not used in a simulation",
            variant="E_TC",
        )
        self.bfield = Leaf(
            qid="0000000002",
            date="1953-12-08 00:00:00",
            description="Not used in a simulation",
            variant="B_TC",
        )
        self.bfield2 = Leaf(
            qid="0000000003",
            date="1953-12-08 00:00:01",
            description="Not used in a simulation",
            variant="B_TC",
        )
        self.bfield_identical_qid = Leaf(
            qid="0000000002",
            date="1953-12-08 00:00:01",
            description="Not used in a simulation",
            variant="B_TC",
        )
        self.inistate = MagicMock()
        self.output = SimulationOutput(
            inputs={"bfield":self.bfield},
            diagnostics={"inistate":self.inistate},
            qid="0000000004",
            date="1953-12-08 00:00:01",
            description="TAG",
            variant="run",
            )
        self.output_identical_qid = SimulationOutput(
            inputs={"bfield":self.bfield},
            diagnostics={"inistate":self.inistate},
            qid="0000000002",
            date="1953-12-08 00:00:01",
            description="TAG",
            variant="run",
            )

    def test_initialization(self):
        """Test that the tree is initially empty with only input categories
        present.
        """
        self.assertFalse(len(self.root._qids))
        for category in input_categories:
            self.assertIn(category, self.root)

    def test_add_input(self):
        """Test adding a new input."""
        leaf = Leaf(
            qid="0123456789",
            date="1953-12-08 00:00:00",
            description="Test description",
            variant="B_TC",
            )
        self.root._add_input(
            leaf, "description", dryrun=False,
            )
        self.assertIn(leaf, self.root.bfield)
        self.assertEqual(leaf._parent, self.root.bfield)

    def test_add_run(self):
        """Test adding a new run."""
        with self.assertRaises(AscotIOException):
            self.root._add_run(
                self.output, "description",
                )

        self.root._add_input(self.bfield)
        self.root._add_run(
            self.output, "description",
            )
        self.assertIn(self.output, self.root)

    def test_add_identical_qid(self):
        """Test adding an input or output when there is data with identical
        QID.
        """
        self.root._add_input(self.bfield)
        with self.assertRaises(AscotIOException):
            self.root._add_input(self.bfield_identical_qid)
        with self.assertRaises(AscotIOException):
            self.root._add_run(
                self.output_identical_qid, "description",
            )

    def test_remove_dataset(self):
        """Test removing dataset."""
        self.root._add_input(self.bfield)
        self.root._add_run(
            self.output, "description",
            )
        self.root.destroy_dataset(self.output)
        self.assertNotIn(self.output, self.root)
        self.root.destroy_dataset(self.bfield.qid)
        self.assertNotIn(self.bfield, self.root.bfield)

    def test_remove_dataset_dependent(self):
        """Test removing input which is being used in a run and contents of
        whole nodes at once."""
        self.root._add_input(self.bfield)
        self.root._add_run(
            self.output, "description",
            )
        self.root.destroy_dataset(self.bfield.qid)
        self.root.destroy_dataset(self.output)

    def test_activate_dataset(self):
        """Test setting dataset active."""
        self.root._add_input(self.bfield)
        self.root._add_input(self.bfield2)
        self.root.activate_dataset(self.bfield2)
        self.assertEqual(self.root.bfield.active, self.bfield2)

    def test_contents(self):
        """Test getting the contents in a string."""
        self.root._add_input(self.bfield)
        self.root._add_input(self.bfield2)
        self.root._add_input(self.efield)
        self.root._add_run(
            self.output, "description",
            )
        lines = self.root.contents.splitlines()
        self.assertIn("Inputs: [only active shown]", lines[0])
        self.assertIn("bfield", lines[3])
        self.assertIn(self.bfield.variant, lines[3])
        self.assertIn(self.bfield.qid, lines[3])
        self.assertIn(self.bfield.date, lines[3])
        self.assertIn("+ 1 other(s)", lines[3])
        self.assertIn(self.bfield.description, lines[4])
        self.assertIn("(no other inputs)", lines[5])
        self.assertIn("*no inputs*", lines[7])


if __name__ == "__main__":
    unittest.main()
