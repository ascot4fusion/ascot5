# pylint: disable=protected-access, no-member
"""Tests that cover the HDF5 IO operations performed by the data tree structure.

Note that the abstract functionality of the tree (`Root`, `ImmutableNode`, and
`MetaDataHolder`) are covered in a separate test. Therefore here we can assume
that the tree operates correctly, and can focus on making sure that the data is
written and the tree is built properly when based on a HDF5 file.
"""
import os
import unittest
from unittest.mock import MagicMock

import h5py
import numpy as np

from a5py import AscotIOException
from a5py.ascot5io.coreio.treedata import DataManager, Leaf, SimulationOutput
from a5py.ascot5io import Ascot5IO


class TestDataManager(unittest.TestCase):

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


class TestTree(unittest.TestCase):

    def setUp(self) -> None:
        with h5py.File("test.h5", "w") as h5:
            category = h5.create_group("bfield")
            category.attrs["active"] = np.bytes_("0000000001")
            group = category.create_group("B_TC_0000000001")
            group.attrs["date"] = np.bytes_("2015-12-23 00:00:00")
            group.attrs["desc"] = np.bytes_("Test")

            h5.attrs["active"] = np.bytes_("0000000002")
            run = h5.create_group("run_0000000002")
            run.attrs["date"] = np.bytes_("2015-12-23 00:00:00")
            run.attrs["desc"] = np.bytes_("Test")

    def tearDown(self) -> None:
        os.remove("test.h5")

    def test_build_tree(self):
        """Test that the tree is initialized properly when data is read from
        the HDF5 file."""
        a5 = Ascot5IO("test.h5")

        self.assertEqual(a5.bfield.active.qid, "0000000001")
        self.assertEqual(a5.bfield.active.date, "2015-12-23 00:00:00")
        self.assertEqual(a5.bfield.active.variant, "B_TC")
        self.assertEqual(a5.bfield.active.description, "Test")

        with self.assertRaises(AscotIOException):
            _ = a5.efield.active

        self.assertEqual(a5.active.qid, "0000000002")
        self.assertEqual(a5.active.date, "2015-12-23 00:00:00")
        self.assertEqual(a5.active.variant, "run")
        self.assertEqual(a5.active.description, "Test")
        self.assertEqual(a5.active.bfield.qid, "0000000001")

        with self.assertRaises(AscotIOException):
            _ = a5.active.efield

    def test_write_tree(self):
        """Test that the tree is written correctly to the HDF5 file.
        """
        a5 = Ascot5IO()

        efield = Leaf(
            qid="0000000001",
            date="1953-12-08 00:00:00",
            description="Not used in a simulation",
            variant="E_TC",
        )
        bfield = Leaf(
            qid="0000000002",
            date="1953-12-08 00:00:00",
            description="Not used in a simulation",
            variant="B_TC",
        )
        bfield2 = Leaf(
            qid="0000000003",
            date="1953-12-08 00:00:01",
            description="Not used in a simulation",
            variant="B_TC",
        )
        inistate = MagicMock()
        output = SimulationOutput(
            inputs={"bfield":bfield},
            diagnostics={"inistate":inistate},
            qid="0000000004",
            date="1953-12-08 00:00:01",
            description="TAG",
            variant="run",
            )
        a5._add_input(efield, store_hdf5=True)
        a5._add_input(bfield, store_hdf5=True)
        a5._add_input(bfield2, store_hdf5=True)
        a5._add_run(output, store_hdf5=True)

    def test_set_active(self):
        """Test that changing the active input or output is reflected on disk.
        """
        pass

    def test_set_description(self):
        """Test that changing the description of an input or output is reflected
        on disk.
        """
        pass

    def test_destroy(self):
        """Test that the datasets are removed from the HDF5 file."""
        pass

    def test_hybrid(self):
        """Test a case where one dataset is on disk and the other is on memory.

        Changing an active group to the one that is not stored on disk should
        not be reflected on HDF5 (the file cannot contain a QID that does not
        exists on file). Likewise, outputs which use memory-only inputs cannot
        be stored before the input has been stored.
        """
        pass


class TestBTC(unittest.TestCase):

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