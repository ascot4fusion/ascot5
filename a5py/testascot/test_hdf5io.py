# pylint: disable=protected-access, no-member
"""Tests that cover the HDF5 IO operations performed by the data tree structure.

Note that the abstract functionality of the tree (`Root`, `ImmutableNode`, and
`MetaDataHolder`) are covered in a separate test. Therefore here we can assume
that the tree operates correctly, and can focus on making sure that the data is
written and the tree is built properly when based on a HDF5 file.
"""
import os
import hashlib
import unittest
from unittest.mock import MagicMock

import h5py
import numpy as np

from a5py import AscotIOException
from a5py.ascot5io.coreio.hdf5interface import HDF5Interface, HDF5Manager
from a5py.ascot5io.coreio.datamanager import DataManager
from a5py.ascot5io.coreio.treestructure import Leaf, SimulationOutput, Root

import a5py.ascot5io.coreio.metadata as metadata

from a5py.testascot.test_coreio import TreeTester

# Make Ascot recognize our dummy variants
metadata.data_variants["bfield"] = metadata.data_variants["bfield"] + ("input",)
metadata.run_variants = metadata.run_variants + ("output",)

class DummyInput(Leaf):
    """Dummy input input variant for testing."""
    pass

class DummyOutput(Leaf):
    """Dummy output variant for testing."""
    pass

class DummyDataTree(Root):
    """Tree that creates dummy input and output data."""

    @classmethod
    def _leaf_factory(cls, metadata, **kwargs):
        """Create dummy data."""
        return super()._leaf_factory(metadata, **kwargs)

class TestDataManager(unittest.TestCase):
    """Tests for the `DataManager` class."""

    def setUp(self) -> None:
        with h5py.File("test.h5", "w") as h5:
            group = h5.create_group("data")
            group.create_dataset("attribute", data=np.array([1]))

    def tearDown(self) -> None:
        os.remove("test.h5")

    def test_memory_access(self):
        dataholder = DataManager(attribute=np.array([1]))
        self.assertEqual(dataholder.get("attribute"), np.array([1]))

    def test_hdf5_access(self):
        dataholder = DataManager(hdf5_filename="test.h5",
                                 path_within_hdf5="data")
        self.assertEqual(dataholder.get("attribute"), np.array([1]))


class TestHDF5Interface(unittest.TestCase):
    """Test the HDF5Interface class."""

    def setUp(self):
        "Create temporary file for testing and open it."
        self.hdf5_interface = HDF5Interface("test.h5")

    def tearDown(self):
        "Close the file and remove it."
        self.hdf5_interface.close()
        os.unlink("test.h5")

    def test_set_node(self):
        """Test setting a node with an active attribute."""
        self.hdf5_interface.set_node("node")
        self.assertIn("node", self.hdf5_interface)

        node = self.hdf5_interface["node"]
        self.assertEqual(node.attrs["active"], np.bytes_(""))

        self.hdf5_interface.set_node("node", active="value")
        self.assertEqual(node.attrs["active"], np.bytes_("value"))

        self.hdf5_interface.set_node("node")
        self.assertEqual(node.attrs["active"], np.bytes_("value"))

    def test_get_node(self):
        """Test getting the active attribute of a node."""
        self.hdf5_interface.set_node("node", active="value")
        value = self.hdf5_interface.get_node(node="node")
        self.assertEqual(value, "value")

    def test_set_dataset(self):
        """Test setting metadata to a dataset."""
        self.hdf5_interface.set_dataset("/", "dataset")
        self.assertIn("dataset", self.hdf5_interface)

        attrs = self.hdf5_interface["dataset"].attrs
        self.assertEqual(attrs["date"], np.bytes_(""))
        self.assertEqual(attrs["desc"], np.bytes_(""))

        self.hdf5_interface.set_dataset(
            "/", "dataset", date="2023-08-12 10:00:00", note="note"
            )
        self.assertEqual(attrs["date"], np.bytes_("2023-08-12 10:00:00"))
        self.assertEqual(attrs["desc"], np.bytes_("note"))

        self.hdf5_interface.set_dataset("/", "dataset")
        self.assertEqual(attrs["date"], np.bytes_("2023-08-12 10:00:00"))
        self.assertEqual(attrs["desc"], np.bytes_("note"))

        self.hdf5_interface.set_dataset("/", "dataset", note="New note")
        self.assertEqual(attrs["date"], np.bytes_("2023-08-12 10:00:00"))
        self.assertEqual(attrs["desc"], np.bytes_("New note"))

        self.hdf5_interface.set_dataset("/", "dataset", extra_attribute="value")
        self.assertIn("extra_attribute", attrs)
        self.assertEqual(attrs["extra_attribute"], np.bytes_("value"))

    def test_get_dataset(self):
        """Test getting metadata from a dataset."""
        self.hdf5_interface.set_node("node")
        self.hdf5_interface.set_dataset("node", "dataset")
        self.hdf5_interface.set_dataset(
            "node", "dataset", date="2023-08-12 10:00:00", note="note"
            )
        date, note = self.hdf5_interface.get_dataset("node", "dataset")
        self.assertEqual(date, "2023-08-12 10:00:00")
        self.assertEqual(note, "note")

        self.hdf5_interface.set_dataset(
            "node", "dataset", date="2023-08-12 10:00:00", note="note",
            extra_attribute="value"
            )
        _, _, extra = self.hdf5_interface.get_dataset("node", "dataset")
        self.assertIn("extra_attribute", extra)
        self.assertEqual(extra["extra_attribute"], "value")


class TestHDF5Manager(unittest.TestCase):
    """Tests for the HDF5Manager class."""

    def setUp(self):
        """Specify metadata and write a file with one input and one output."""
        self.testfile = "test.h5"
        self.emptyfile = "empty.h5"
        self.date = "2008-10-20 00:00:00"
        self.note = "Note"

        self.inputqid = "0123456789"
        self.inputvar = "input"
        self.inputnode = "bfield"
        self.outputqid = "1234567890"
        self.outputvar = "output"
        with HDF5Interface(self.testfile) as h5:
            h5.set_node(self.inputnode, self.inputqid)
            h5.set_node("/", self.outputqid)
            h5.set_dataset(
                self.inputnode, f"{self.inputvar}_{self.inputqid}",
                date=self.date, note=self.note,
                )
            h5.set_dataset(
                "/", f"{self.outputvar}_{self.outputqid}", date=self.date,
                note=self.note, **{self.inputnode:self.inputqid},
            )

    def tearDown(self):
        """Remove both the test file and empty file (if one was created by
        a test)."""
        os.unlink(self.testfile)
        try:
            os.unlink(self.emptyfile)
        except FileNotFoundError:
            pass

    def test_init(self):
        """Test initialization from an empty file and file with data."""
        with self.assertRaises(FileNotFoundError):
            HDF5Manager(self.emptyfile, file_exists=True)
        with self.assertRaises(FileExistsError):
            HDF5Manager(self.testfile, file_exists=False)

        HDF5Manager(self.emptyfile, file_exists=False)
        with HDF5Interface(self.emptyfile) as h5:
            for node in list(metadata.input_categories) + ["/"]:
                self.assertIn(node, h5)
                self.assertEqual(h5.get_node(node), "")

        HDF5Manager(self.testfile, file_exists=True)
        with HDF5Interface(self.testfile) as h5:
            for node in list(metadata.input_categories) + ["/"]:
                self.assertIn(node, h5)
                if node == self.inputnode:
                    self.assertEqual(h5.get_node(node), self.inputqid)
                elif node == "/":
                    self.assertEqual(h5.get_node(node), self.outputqid)
                else:
                    self.assertEqual(h5.get_node(node), "")

    def test_read_node(self):
        """Test reading metadata field from a node."""
        manager = HDF5Manager(self.testfile, file_exists=True)
        active, qids, variants = manager.read_node(self.inputnode)
        self.assertEqual(active, self.inputqid)
        self.assertEqual(qids, [self.inputqid])
        self.assertEqual(variants, [self.inputvar])

        active, qids, variants = manager.read_node("root")
        self.assertEqual(active, self.outputqid)
        self.assertEqual(qids, [self.outputqid])
        self.assertEqual(variants, [self.outputvar])

    def test_set_active(self):
        """Test setting active QId to a node."""
        manager = HDF5Manager(self.testfile, file_exists=True)

        manager.set_active("root", "0000000001")
        active, _, _ = manager.read_node("root")
        self.assertEqual(active, "0000000001")

        manager.set_active("bfield", "0000000002")
        active, _, _ = manager.read_node("bfield")
        self.assertEqual(active, "0000000002")

    def test_read_input(self):
        """Test reading input dataset."""
        manager = HDF5Manager(self.testfile, file_exists=True)
        date, note = manager.read_input(self.inputqid, self.inputvar)
        self.assertEqual(date, self.date)
        self.assertEqual(note, self.note)

    def test_read_output(self):
        """Test reading output dataset."""
        manager = HDF5Manager(self.testfile, file_exists=True)
        out = manager.read_output(self.outputqid, self.outputvar)
        date, note, inputqids = out
        self.assertEqual(date, self.date)
        self.assertEqual(note, self.note)
        self.assertEqual(inputqids, [self.inputqid])

    def test_write_input(self):
        """Test writing input dataset."""
        manager = HDF5Manager(self.emptyfile, file_exists=False)
        meta = metadata.MetaData(
            qid=self.inputqid, date=self.date, note=self.note,
            variant=self.inputvar,
            )
        manager.write_input(meta)
        date, note = manager.read_input(meta.qid, meta.variant)
        self.assertEqual(date, self.date)
        self.assertEqual(note, self.note)

    def test_write_output(self):
        """Test writing output dataset."""
        manager = HDF5Manager(self.emptyfile, file_exists=False)
        meta = metadata.MetaData(
            qid=self.outputqid, date=self.date, note=self.note,
            variant=self.outputvar,
            )
        manager.write_output(meta, {self.inputnode:self.inputqid})
        out = manager.read_output(meta.qid, meta.variant)
        date, note, inputqids = out
        self.assertEqual(date, self.date)
        self.assertEqual(note, self.note)
        self.assertEqual(inputqids, [self.inputqid])

    def test_set_note(self):
        """Test modifying the note on the file."""
        newnote = "New note"
        manager = HDF5Manager(self.testfile, file_exists=True)

        manager.set_note(self.inputqid, self.inputvar, newnote)
        _, note = manager.read_input(self.inputqid, self.inputvar)
        self.assertEqual(note, newnote)

        manager.set_note(self.outputqid, self.outputvar, newnote)
        _, note, _ = manager.read_output(self.outputqid, self.outputvar)
        self.assertEqual(note, newnote)

    def test_remove_dataset(self):
        """Test removing a dataset."""
        manager = HDF5Manager(self.testfile, file_exists=True)
        with HDF5Interface(self.testfile) as h5:
            name = f"{self.inputvar}_{self.inputqid}"
            self.assertIn(name, h5[self.inputnode])
            name = f"{self.outputvar}_{self.outputqid}"
            self.assertIn(name, h5)

        manager.remove_dataset(self.inputqid, self.inputvar)
        with HDF5Interface(self.testfile) as h5:
            name = f"{self.inputvar}_{self.inputqid}"
            self.assertNotIn(name, h5[self.inputnode])

        manager.remove_dataset(self.outputqid, self.outputvar)
        with HDF5Interface(self.testfile) as h5:
            name = f"{self.outputvar}_{self.outputqid}"
            self.assertNotIn(name, h5)

    def test_repack(self):
        """Test repacking the file."""
        manager = HDF5Manager(self.testfile, file_exists=True)
        with HDF5Interface(self.testfile) as h5:
            name = f"{self.inputvar}_{self.inputqid}"
            h5[self.inputnode][name].create_dataset(
                "data", data=np.ones(10000,),
                )

        originalsize = os.path.getsize(self.testfile)
        manager.remove_dataset(self.inputqid, self.inputvar)
        unpackedsize = os.path.getsize(self.testfile)
        manager.repack()
        repackedsize = os.path.getsize(self.testfile)

        self.assertLessEqual(unpackedsize, originalsize)
        self.assertLess(repackedsize, unpackedsize)


class TestTree(unittest.TestCase):
    """Test tree operations related to HDF5."""

    def setUp(self):
        """Prepare a mock of HDF5Manager that manages a (immutable) file with
        fixed datasets."""
        self.testfile = "test.h5"
        self.emptyfile = "empty.h5"
        self.date = TreeTester.DATETODAY
        self.note = TreeTester.DEFAULTNOTE

        self.inputqid = TreeTester.QID1
        self.inputvar = "input"
        self.inputnode = "bfield"
        self.outputqid = TreeTester.QID2
        self.outputvar = "output"

        self.hdf5manager = MagicMock()

        def read_node(node):
            if node == "root":
                return (self.outputqid, [self.outputqid], [self.outputvar])
            elif node == self.inputnode:
                return (self.inputqid, [self.inputqid], [self.inputvar])
            return (None, [], [])
        self.hdf5manager.read_node.side_effect = read_node

        def read_input(qid, var):
            return (self.date, self.note)
        self.hdf5manager.read_input.side_effect = read_input

        def read_output(qid, var):
            return (self.date, self.note, [self.inputqid])
        self.hdf5manager.read_output.side_effect = read_output

    def tearDown(self):
        pass

    def test_initialization(self):
        """Test that the tree is initialized properly when data is read from
        the HDF5 file."""
        a5 = Root()
        with a5._modify_attributes():
            a5._set_hdf5manager(self.hdf5manager)
            a5._init_from_hdf5()

        self.assertEqual(a5.bfield.active.qid, TreeTester.QID1)
        self.assertEqual(a5.bfield.active.date, TreeTester.DATETODAY)
        self.assertEqual(a5.bfield.active.variant, "input")
        self.assertEqual(a5.bfield.active.note, TreeTester.DEFAULTNOTE)

        with self.assertRaises(AscotIOException):
            _ = a5.efield.active

        self.assertEqual(a5.active.qid, TreeTester.QID2)
        self.assertEqual(a5.active.date, TreeTester.DATETODAY)
        self.assertEqual(a5.active.variant, "output")
        self.assertEqual(a5.active.note, TreeTester.DEFAULTNOTE)
        self.assertEqual(a5.active.bfield.qid, TreeTester.QID1)

        with self.assertRaises(AscotIOException):
            _ = a5.active.efield

    def test_add_input(self):
        """Test adding a new input."""
        a5 = Root()

        meta = metadata.MetaData(qid=TreeTester.QID1, date=TreeTester.DATETODAY,
                                 note=TreeTester.DEFAULTNOTE, variant="input")
        with self.assertRaises(
            AscotIOException,
            msg="Error was not raised when trying to write to nonexistent file."
            ):
            a5._add_input_dataset(meta, store_hdf5=True)

        with a5._modify_attributes():
            a5._set_hdf5manager(self.hdf5manager)
            a5._init_from_hdf5()

        a5._add_input_dataset(meta)
        meta = metadata.MetaData(qid=TreeTester.QID1, date=TreeTester.DATETODAY,
                                 note=TreeTester.DEFAULTNOTE, variant="input")
        try:
            a5._add_input_dataset(meta, store_hdf5=True)
        except:
            # Can't write dataset with same QID
            pass
        file_contents_changed = self.file_hash("empty.h5") != hash
        self.assertFalse(file_contents_changed,
                         "File was modified when writing input failed.")

    # def test_add_simulation_output(self):
    #     """Test adding a new run."""
    #     with self.assertRaises(AscotIOException):
    #         self.root._add_simulation_output(
    #             self.metadata["output"], [self.inistate],
    #             [self.metadata["bfield"].qid], "note",
    #             )

    #     self.root._add_input_dataset(self.metadata["bfield"])
    #     output = self.root._add_simulation_output(
    #         self.metadata["output"], [self.inistate],
    #         [self.metadata["bfield"].qid], "note",
    #         )
    #     self.assertIn(output, self.root)

    # def test_add_identical_qid(self):
    #     """Test adding an input or output when there is data with identical
    #     QID.
    #     """
    #     self.root._add_input_dataset(self.metadata["bfield"])
    #     with self.assertRaises(AscotIOException):
    #         self.root._add_input_dataset(self.metadata["bfield_identical_qid"])
    #     with self.assertRaises(AscotIOException):
    #         self.root._add_simulation_output(
    #             self.metadata["output_identical_qid"], [self.inistate],
    #             [self.metadata["bfield"].qid], "note",
    #         )

    # def test_remove_dataset(self):
    #     """Test removing dataset."""
    #     bfield = self.root._add_input_dataset(self.metadata["bfield"])
    #     output = self.root._add_simulation_output(
    #             self.metadata["output"], [self.inistate],
    #             [self.metadata["bfield"].qid], "note",
    #         )
    #     self.root.destroy_dataset(output)
    #     self.assertNotIn(output, self.root)
    #     self.root.destroy_dataset(bfield.qid)
    #     self.assertNotIn(bfield, self.root.bfield)

    # def test_remove_dataset_dependent(self):
    #     """Test removing input which is being used in a run and contents of
    #     whole nodes at once."""
    #     bfield = self.root._add_input_dataset(self.metadata["bfield"])
    #     output = self.root._add_simulation_output(
    #             self.metadata["output"], [self.inistate],
    #             [self.metadata["bfield"].qid], "note",
    #         )
    #     self.root.destroy_dataset(bfield.qid)
    #     self.root.destroy_dataset(output)

    # def test_activate_dataset(self):
    #     """Test setting dataset active."""
    #     self.root._add_input_dataset(self.metadata["bfield"])
    #     bfield2 = self.root._add_input_dataset(self.metadata["bfield2"])
    #     self.root.activate_dataset(bfield2)
    #     self.assertEqual(self.root.bfield.active, bfield2)

    # def test_set_active(self):
    #     """Test that changing the active input or output is reflected on disk.
    #     """
    #     pass

    # def test_set_note(self):
    #     """Test that changing the note of an input or output is reflected
    #     on disk.
    #     """
    #     pass

    # def test_destroy(self):
    #     """Test that the datasets are removed from the HDF5 file."""
    #     pass

    # def test_hybrid(self):
    #     """Test a case where one dataset is on disk and the other is on memory.

    #     Changing an active group to the one that is not stored on disk should
    #     not be reflected on HDF5 (the file cannot contain a QID that does not
    #     exists on file). Likewise, outputs which use memory-only inputs cannot
    #     be stored before the input has been stored.
    #     """
    #     pass


class TestInputDataset(unittest.TestCase):

    def test_read_write_hdf5(self):
        """
        """
        pass

    def test_init(self):
        """
        """
        pass

    def test_eval(self):
        """
        """
        pass


if __name__ == "__main__":
    unittest.main()