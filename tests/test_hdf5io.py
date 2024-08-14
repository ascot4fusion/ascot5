# pylint: disable=protected-access, no-member
"""Tests that cover the HDF5 IO operations performed by the data tree structure.

We assume that the abstract functionality of the tree operates properly (since
it is covered by a separate test) and only focus making sure that HDF5 IO
operations are performed correctly.

Tested classes: HDF5Interface, HDF5Manager, and Root (that it calls HDF5Manager
properly).
"""
import os
import unittest
from unittest.mock import MagicMock

import h5py
import numpy as np

from a5py import AscotIOException

import a5py.ascot5io.coreio.metadata as metadata
from a5py.ascot5io.coreio.metadata import MetaData
from a5py.ascot5io.coreio.hdf5interface import HDF5Interface, HDF5Manager
from a5py.ascot5io.coreio.datamanager import DataManager
from a5py.ascot5io.coreio.treestructure import Leaf, Root

FNTEST = "test.h5"
FNEMPTY = "empty.h5"
QID1 = "2991786571"
QID2 = "9753987342"
QID3 = "4404229430"
QID4 = "0963810214"
QID5 = "5960585966"
DATE_FRI = "1997-08-29 02:14:00"
DATE_SAT = "1997-08-30 02:14:00"
DATE_SUN = "1997-08-31 02:14:00"
INPUTVAR = "input"
OUTPUTVAR = "output"
CATEGORY = "wall"
DATE = DATE_FRI
NOTE = "Let off some steam Bennett"

# Make Ascot recognize our dummy variants
metadata.data_variants[CATEGORY] = (
    metadata.data_variants[CATEGORY] + (INPUTVAR,)
)
metadata.run_variants = metadata.run_variants + (OUTPUTVAR,)

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
        self.hdf5_interface = HDF5Interface(FNTEST)

    def tearDown(self):
        "Close the file and remove it."
        self.hdf5_interface.close()
        os.unlink(FNTEST)

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

        self.hdf5_interface.set_dataset("/", "dataset", date=DATE, note=NOTE)
        self.assertEqual(attrs["date"], np.bytes_(DATE))
        self.assertEqual(attrs["desc"], np.bytes_(NOTE))

        self.hdf5_interface.set_dataset("/", "dataset")
        self.assertEqual(attrs["date"], np.bytes_(DATE))
        self.assertEqual(attrs["desc"], np.bytes_(NOTE))

        self.hdf5_interface.set_dataset("/", "dataset", note="New note")
        self.assertEqual(attrs["date"], np.bytes_(DATE))
        self.assertEqual(attrs["desc"], np.bytes_("New note"))

        self.hdf5_interface.set_dataset("/", "dataset", extra_attribute="value")
        self.assertIn("extra_attribute", attrs)
        self.assertEqual(attrs["extra_attribute"], np.bytes_("value"))

    def test_get_dataset(self):
        """Test getting metadata from a dataset."""
        self.hdf5_interface.set_node("node")
        self.hdf5_interface.set_dataset("node", "dataset")
        self.hdf5_interface.set_dataset("node", "dataset", date=DATE, note=NOTE)
        date, note = self.hdf5_interface.get_dataset("node", "dataset")
        self.assertEqual(date, DATE)
        self.assertEqual(note, NOTE)

        self.hdf5_interface.set_dataset(
            "node", "dataset", date=DATE, note=NOTE, extra_attribute="value",
            )
        _, _, extra = self.hdf5_interface.get_dataset("node", "dataset")
        self.assertIn("extra_attribute", extra)
        self.assertEqual(extra["extra_attribute"], "value")


class TestHDF5Manager(unittest.TestCase):
    """Tests for the HDF5Manager class."""

    def setUp(self):
        """Write a file with one input and one output."""
        with HDF5Interface(FNTEST) as h5:
            h5.set_node(CATEGORY, QID1)
            h5.set_node("/", QID2)
            h5.set_dataset(
                CATEGORY, f"{INPUTVAR}_{QID1}", date=DATE, note=NOTE,
                )
            h5.set_dataset(
                "/", f"{OUTPUTVAR}_{QID2}", date=DATE, note=NOTE,
                **{CATEGORY:QID1},
            )

    def tearDown(self):
        """Remove both the test file and empty file (if one was created by
        a test)."""
        os.unlink(FNTEST)
        try:
            os.unlink(FNEMPTY)
        except FileNotFoundError:
            pass

    def test_init(self):
        """Test initialization from an empty file and file with data."""
        with self.assertRaises(FileNotFoundError):
            HDF5Manager(FNEMPTY, file_exists=True)
        with self.assertRaises(FileExistsError):
            HDF5Manager(FNTEST, file_exists=False)

        HDF5Manager(FNEMPTY, file_exists=False)
        with HDF5Interface(FNEMPTY) as h5:
            for node in list(metadata.input_categories) + ["/"]:
                self.assertIn(node, h5)
                self.assertEqual(h5.get_node(node), "")

        HDF5Manager(FNTEST, file_exists=True)
        with HDF5Interface(FNTEST) as h5:
            for node in list(metadata.input_categories) + ["/"]:
                self.assertIn(node, h5)
                if node == CATEGORY:
                    self.assertEqual(h5.get_node(node), QID1)
                elif node == "/":
                    self.assertEqual(h5.get_node(node), QID2)
                else:
                    self.assertEqual(h5.get_node(node), "")

    def test_read_node(self):
        """Test reading metadata field from a node."""
        manager = HDF5Manager(FNTEST, file_exists=True)
        active, qids, variants = manager.read_node(CATEGORY)
        self.assertEqual(active, QID1)
        self.assertEqual(qids, [QID1])
        self.assertEqual(variants, [INPUTVAR])

        active, qids, variants = manager.read_node("root")
        self.assertEqual(active, QID2)
        self.assertEqual(qids, [QID2])
        self.assertEqual(variants, [OUTPUTVAR])

    def test_set_active(self):
        """Test setting active QID to a node."""
        manager = HDF5Manager(FNTEST, file_exists=True)

        manager.set_active("root", QID1)
        active, _, _ = manager.read_node("root")
        self.assertEqual(active, QID1)

        manager.set_active(CATEGORY, QID2)
        active, _, _ = manager.read_node(CATEGORY)
        self.assertEqual(active, QID2)

    def test_read_input(self):
        """Test reading input dataset."""
        manager = HDF5Manager(FNTEST, file_exists=True)
        date, note = manager.read_input(QID1, INPUTVAR)
        self.assertEqual(date, DATE)
        self.assertEqual(note, NOTE)

    def test_read_output(self):
        """Test reading output dataset."""
        manager = HDF5Manager(FNTEST, file_exists=True)
        out = manager.read_output(QID2, OUTPUTVAR)
        date, note, inputqids = out
        self.assertEqual(date, DATE)
        self.assertEqual(note, NOTE)
        self.assertEqual(inputqids, [QID1])

    def test_write_input(self):
        """Test writing input dataset."""
        manager = HDF5Manager(FNEMPTY, file_exists=False)
        meta = metadata.MetaData(
            qid=QID1, date=DATE, note=NOTE,
            variant=INPUTVAR,
            )
        manager.write_input(meta)
        date, note = manager.read_input(meta.qid, meta.variant)
        self.assertEqual(date, DATE)
        self.assertEqual(note, NOTE)

    def test_write_output(self):
        """Test writing output dataset."""
        manager = HDF5Manager(FNEMPTY, file_exists=False)
        meta = metadata.MetaData(
            qid=QID2, date=DATE, note=NOTE,
            variant=OUTPUTVAR,
            )
        manager.write_output(meta, {CATEGORY:QID1})
        out = manager.read_output(meta.qid, meta.variant)
        date, note, inputqids = out
        self.assertEqual(date, DATE)
        self.assertEqual(note, NOTE)
        self.assertEqual(inputqids, [QID1])

    def test_set_note(self):
        """Test modifying the note on the file."""
        newnote = "New note"
        manager = HDF5Manager(FNTEST, file_exists=True)

        manager.set_note(QID1, INPUTVAR, newnote)
        _, note = manager.read_input(QID1, INPUTVAR)
        self.assertEqual(note, newnote)

        manager.set_note(QID2, OUTPUTVAR, newnote)
        _, note, _ = manager.read_output(QID2, OUTPUTVAR)
        self.assertEqual(note, newnote)

    def test_remove_dataset(self):
        """Test removing a dataset."""
        manager = HDF5Manager(FNTEST, file_exists=True)
        with HDF5Interface(FNTEST) as h5:
            name = f"{INPUTVAR}_{QID1}"
            self.assertIn(name, h5[CATEGORY])
            name = f"{OUTPUTVAR}_{QID2}"
            self.assertIn(name, h5)

        manager.remove_dataset(QID1, INPUTVAR)
        with HDF5Interface(FNTEST) as h5:
            name = f"{INPUTVAR}_{QID1}"
            self.assertNotIn(name, h5[CATEGORY])

        manager.remove_dataset(QID2, OUTPUTVAR)
        with HDF5Interface(FNTEST) as h5:
            name = f"{OUTPUTVAR}_{QID2}"
            self.assertNotIn(name, h5)

    def test_repack(self):
        """Test repacking the file."""
        manager = HDF5Manager(FNTEST, file_exists=True)
        with HDF5Interface(FNTEST) as h5:
            name = f"{INPUTVAR}_{QID1}"
            h5[CATEGORY][name].create_dataset(
                "data", data=np.ones(10000,),
                )

        originalsize = os.path.getsize(FNTEST)
        manager.remove_dataset(QID1, INPUTVAR)
        unpackedsize = os.path.getsize(FNTEST)
        manager.repack()
        repackedsize = os.path.getsize(FNTEST)

        self.assertLessEqual(unpackedsize, originalsize)
        self.assertLess(repackedsize, unpackedsize)


class TestTreeWithHDF5(unittest.TestCase):
    """Test tree operations related to HDF5."""

    def setUp(self):
        """Prepare a mock of HDF5Manager that manages a (immutable) file with
        fixed datasets."""
        self.hdf5manager = MagicMock()

        def read_node(node):
            if node == "root":
                return (QID2, [QID2], [OUTPUTVAR])
            elif node == CATEGORY:
                return (QID1, [QID1], [INPUTVAR])
            return (None, [], [])
        self.hdf5manager.read_node.side_effect = read_node

        def read_input(qid, var):
            return (DATE, NOTE)
        self.hdf5manager.read_input.side_effect = read_input

        def read_output(qid, var):
            return (DATE, NOTE, [QID1])
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

        self.assertEqual(a5[CATEGORY].active.qid, QID1)
        self.assertEqual(a5[CATEGORY].active.date, DATE)
        self.assertEqual(a5[CATEGORY].active.variant, INPUTVAR)
        self.assertEqual(a5[CATEGORY].active.note, NOTE)

        with self.assertRaises(AscotIOException):
            _ = a5.efield.active

        self.assertEqual(a5.active.qid, QID2)
        self.assertEqual(a5.active.date, DATE)
        self.assertEqual(a5.active.variant, OUTPUTVAR)
        self.assertEqual(a5.active.note, NOTE)
        self.assertEqual(a5.active[CATEGORY].qid, QID1)

        with self.assertRaises(AscotIOException):
            _ = a5.active.efield

    def test_add_input(self):
        """Test adding a new input."""
        a5 = Root()
        meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=INPUTVAR)
        with self.assertRaises(
            AscotIOException,
            msg="Error was not raised when trying to write to nonexistent file."
            ):
            a5._add_input_dataset(meta, store_hdf5=True)

        with a5._modify_attributes():
            a5._set_hdf5manager(self.hdf5manager)
            a5._init_from_hdf5()

        a5._add_input_dataset(meta)
        self.hdf5manager.write_input.assert_called_once_with(meta)
        self.hdf5manager.reset_mock()

        meta = MetaData(qid=QID4, date=DATE, note=NOTE, variant=INPUTVAR)
        dataset = a5._add_input_dataset(meta, store_hdf5=False)
        self.hdf5manager.write_input.assert_not_called()

        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=INPUTVAR)
        try:
            a5._add_input_dataset(meta)
        except:
            # Can't write dataset with same QID
            pass
        self.hdf5manager.write_input.assert_not_called()

    def test_add_simulation_output(self):
        """Test adding a new run."""
        a5 = Root()
        meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=OUTPUTVAR)
        with self.assertRaises(
            AscotIOException,
            msg="Error was not raised when trying to write to nonexistent file."
            ):
            a5._add_simulation_output(meta, [], [QID1], store_hdf5=True)

        with a5._modify_attributes():
            a5._set_hdf5manager(self.hdf5manager)
            a5._init_from_hdf5()

        a5._add_simulation_output(meta, [], [QID1])
        self.hdf5manager.write_output.assert_called_once_with(
            meta, {CATEGORY:QID1},
            )
        self.hdf5manager.reset_mock()

        meta = MetaData(qid=QID4, date=DATE, note=NOTE, variant=OUTPUTVAR)
        a5._add_simulation_output(meta, [], [QID1], store_hdf5=False)
        self.hdf5manager.write_output.assert_not_called()

        meta = MetaData(qid=QID1, date=DATE, note=NOTE, variant=OUTPUTVAR)
        try:
            a5._add_simulation_output(meta, [], [QID1])
        except:
            # Can't write dataset with same QID
            pass
        self.hdf5manager.write_output.assert_not_called()

    def test_add_simulation_output_hybrid(self):
        """Test adding simulation output when the used input is not stored on
        disk.
        """
        a5 = Root()
        with a5._modify_attributes():
            a5._set_hdf5manager(self.hdf5manager)
            a5._init_from_hdf5()

        meta = MetaData(qid=QID3, date=DATE, note=NOTE, variant=INPUTVAR)
        a5._add_input_dataset(meta, store_hdf5=False)

        meta = MetaData(qid=QID4, date=DATE, note=NOTE, variant=OUTPUTVAR)
        with self.assertRaises(AscotIOException):
            a5._add_simulation_output(meta, [], [QID3])
        self.hdf5manager.write_output.assert_not_called()

        a5._add_simulation_output(meta, [], [QID3], store_hdf5=False)
        self.hdf5manager.write_output.assert_not_called()

    def test_remove_dataset(self):
        """Test removing dataset."""
        a5 = Root()
        with a5._modify_attributes():
            a5._set_hdf5manager(self.hdf5manager)
            a5._init_from_hdf5()

        try:
            a5.destroy_dataset(QID1)
        except AscotIOException:
            # Can't destroy dataset that is being used
            pass
        self.hdf5manager.remove_dataset.assert_not_called()

        a5.destroy_dataset(QID2)
        self.hdf5manager.remove_dataset.assert_called_once_with(QID2, OUTPUTVAR)
        self.hdf5manager.reset_mock()

        a5[CATEGORY].active.destroy()
        self.hdf5manager.remove_dataset.assert_called_once_with(QID1, INPUTVAR)

    def test_remove_node(self):
        """Test removing all datasets within a node."""
        pass

    def test_reset_active(self):
        """Test that the active """
        pass

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
    #     self.assertEqual(self.root[CATEGORY].active, bfield2)

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