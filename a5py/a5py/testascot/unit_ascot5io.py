"""Test Ascot5IO.
"""
import numpy as np
import h5py
import time
import unittest
import datetime
import subprocess
import importlib

from a5py import Ascot
from a5py.exceptions import AscotIOException
from a5py.ascot5io._iohelpers import fileapi
from a5py.ascot5io.state import State

class TestAscot5IO(unittest.TestCase):
    """Class to test `ascot5io` module.
    """

    def setUp(self):
        """Initialize name for the file used in testing.
        """
        self.testfilename = "unittest.h5"

    def test_io_operations(self):
        """Test Ascot5IO operations that affect the contents of the HDF5 file.
        """

        # Start with an empty object.
        a5 = Ascot(self.testfilename, create=True)
        ETC = {"exyz" : np.array([0,0,0])}

        # Add one group
        grp = a5.data.create_input("E_TC", ETC)

        # Verify it was added to Ascot5IO
        self.assertTrue("efield" in a5.data,
                        "Parent group not found in Ascot5IO.")
        self.assertTrue(grp in a5.data["efield"],
                        "Data group not found in Ascot5IO.")

        # Verify it was added to HDF5.
        with h5py.File(self.testfilename, "r") as h5:
            self.assertTrue("efield" in h5,
                            "Parent group not found in HDF5.")
            self.assertTrue(grp in h5["efield"],
                            "Data group not found in HDF5.")

        # Verify it became the active group.
        self.assertTrue("active" in a5.data["efield"],
                        "No active field in Ascot5IO.")
        self.assertEqual(a5.data["efield"]["active"], a5.data["efield"][grp],
                         "Active group not properly set.")

        # Remove group
        a5.data["efield"][grp].destroy(repack=False)

        # Since it was the only group, the parent should have been removed.
        self.assertFalse("efield" in a5.data,
                        "Parent group not removed from Ascot5IO.")

        # Verify the parent was also removed from the HDF5.
        with h5py.File(self.testfilename, "r") as h5:
            self.assertFalse("efield" in h5,
                            "Parent group not removed from HDF5.")

        # Add three new groups (sleep between so that dates are not same)
        grp1 = a5.data.create_input("E_TC", ETC)
        time.sleep(1.0)
        grp2 = a5.data.create_input("E_TC", ETC)
        time.sleep(1.0)
        grp3 = a5.data.create_input("E_TC", ETC)

        # Groups should be ordered by date
        dates = a5.data["efield"].get_contents()["date"]
        self.assertEqual(dates, sorted(dates)[::-1],
                         "Groups are not sorted by date.")

        # First group should be set as active
        self.assertEqual(a5.data["efield"]["active"], a5.data["efield"][grp1],
                         "Active group not properly set.")

        # Switch active
        a5.data["efield"][grp2].activate()
        self.assertEqual(a5.data["efield"]["active"], a5.data["efield"][grp2],
                         "Active group not properly set.")

        # Remove active group
        a5.data["efield"][grp2].destroy(repack=False)

        # Parent group and active field should still be present
        with h5py.File(self.testfilename, "r") as h5:
            self.assertTrue("efield" in h5,
                            "Parent group not found in HDF5.")
        self.assertTrue("efield" in a5.data,
                        "Parent group not found in Ascot5IO.")
        self.assertTrue("active" in a5.data["efield"],
                        "No active field in Ascot5IO.")

        # Group itself should have been removed
        self.assertFalse(grp2 in a5.data["efield"],
                        "Data group not removed from Ascot5IO.")
        with h5py.File(self.testfilename, "r") as h5:
            self.assertFalse(grp2 in h5["efield"],
                            "Data group not removed from HDF5.")

        # The most recent group should be active
        self.assertEqual(a5.data["efield"]["active"], a5.data["efield"][grp3],
                         "Active group not properly set.")

        # Remove parent
        a5.data["efield"].destroy(repack=False)

        # Parent group should have been removed now
        self.assertFalse("efield" in a5.data,
                        "Parent group not removed from Ascot5IO.")
        with h5py.File(self.testfilename, "r") as h5:
            self.assertFalse("efield" in h5,
                            "Parent group not removed from HDF5.")

    def test_metadata(self):
        """Test object metadata.
        """
        # Create group and store its type, QID, and date.
        a5 = Ascot(self.testfilename, create=True)
        ETC = {"exyz" : np.array([0,0,0])}
        date1 = datetime.datetime.strptime(str(datetime.datetime.now())[0:19],
                                           "%Y-%m-%d %H:%M:%S")
        grp = a5.data.create_input("E_TC", ETC)
        qid = grp.split("_")[-1]
        gtype = grp[:-11]
        date2 = datetime.datetime.strptime(str(datetime.datetime.now())[0:19],
                                           "%Y-%m-%d %H:%M:%S")

        # Set description
        desc = "Test"
        a5.data.efield.active.set_desc(desc)

        # Check that the values in HDF5 correspond to ones that were set.
        self.assertEqual(a5.data.efield.active.get_qid(), qid,
                         "QID invalid.")
        self.assertEqual(a5.data.efield.active.get_type(), gtype,
                         "Type invalid.")
        date = datetime.datetime.strptime(a5.data.efield.active.get_date(),
                                          "%Y-%m-%d %H:%M:%S")
        self.assertTrue(date >= date1 and date <=date2,
                         "Date invalid.")
        self.assertEqual(a5.data.efield.active.get_desc(), desc,
                         "Desc invalid.")

    def test_reference(self):
        """Test object referencing in treeview.
        """
        # Create object with multiple groups.
        a5 = Ascot(self.testfilename, create=True)
        ETC = {"exyz" : np.array([0,0,0])}
        grp1 = a5.data.create_input("E_TC", ETC)
        time.sleep(1.0)
        grp2 = a5.data.create_input("E_TC", ETC)
        time.sleep(1.0)
        grp3 = a5.data.create_input("E_TC", ETC)

        # Test with parent group that dict and attribute references are equal
        self.assertEqual(a5.data.efield, a5.data["efield"],
                         "Referencing parent group with dict and attr failed.")

        # Reference by name
        qid   = a5.data.efield[grp1].get_qid()
        gtype = a5.data.efield[grp1].get_type()
        self.assertEqual(grp1, gtype + "_" + qid,
                         "Reference by name failed.")

        # Reference by QID
        self.assertEqual(a5.data.efield[grp1], a5.data.efield["q"+qid],
                         "Reference by QID failed.")

        # Reference by tag: same tags should have running index from newest
        self.assertEqual(a5.data.efield[grp3], a5.data.efield["TAG_0"],
                         "Reference by tag failed.")
        self.assertEqual(a5.data.efield[grp2], a5.data.efield["TAG_1"],
                         "Reference by tag failed.")
        self.assertEqual(a5.data.efield[grp1], a5.data.efield["TAG_2"],
                         "Reference by tag failed.")

        # Reference by tag: same tag running index should update dynamically
        a5.data.efield[grp2].set_desc("test test")
        self.assertEqual(a5.data.efield[grp3], a5.data.efield["TAG_0"],
                         "Reference by tag failed.")
        self.assertEqual(a5.data.efield[grp1], a5.data.efield["TAG_1"],
                         "Reference by tag failed.")

        # Reference by tag: tag is the first word and on uppercase
        self.assertEqual(a5.data.efield[grp2], a5.data.efield["TEST"],
                         "Reference by tag failed.")

        # Reference by tag: no more than ten chars and special chars omitted
        a5.data.efield[grp2].set_desc("N1C3T^GBUT_THISISOMITTED")
        self.assertFalse("TEST" in a5.data.efield,
                         "Reference by tag failed.")
        self.assertEqual(a5.data.efield[grp2], a5.data.efield["N1C3TGBUT"],
                         "Reference by tag failed.")

        # Reference by tag: no tag reverts to default and if first char is num.
        a5.data.efield[grp2].set_desc("")
        self.assertEqual(a5.data.efield[grp2], a5.data.efield["TAG_1"],
                         "Reference by tag failed.")
        a5.data.efield[grp2].set_desc("1TAGGY")
        self.assertEqual(a5.data.efield[grp2], a5.data.efield["TAG_1"],
                         "Reference by tag failed.")

    def test_resultnode(self):
        """Test result nodes.
        """
        # Create file and two inputs
        a5 = Ascot(self.testfilename, create=True)
        ETC = {"exyz" : np.array([0,0,0])}
        grp1 = a5.data.create_input("E_TC", ETC)
        time.sleep(1.0)
        grp2 = a5.data.create_input("E_TC", ETC)

        # Create two runs; one for each input
        qid1 = a5.data.efield[grp1].get_qid()
        qid2 = a5.data.efield[grp2].get_qid()
        with h5py.File(self.testfilename, "a") as h5:
            group = fileapi.add_group(h5, "results", "run", desc="RUN1")
            group.attrs["qid_efield"] = np.string_(qid1)
            time.sleep(1.0)
            group = fileapi.add_group(h5, "results", "run", desc="RUN2")
            group.attrs["qid_efield"] = np.string_(qid2)

        a5 = Ascot(self.testfilename)
        run1 = a5.data["RUN1"].get_qid()
        run2 = a5.data["RUN2"].get_qid()

        # Test print functions (that they don't raise errors)
        a5.data.efield.ls()
        a5.data.active.ls()
        a5.data.ls()

        # Run-to-input reference
        self.assertEqual(a5.data["RUN1"].efield.get_qid(), qid1,
                         "Run does not reference input properly.")
        self.assertEqual(a5.data["RUN2"].efield.get_qid(), qid2,
                         "Run does not reference input properly.")

        # Input-to-run reference
        self.assertEqual(a5.data.get_runs(qid1)[0], run1,
                         "Failed to find the run where input was used.")
        self.assertEqual(a5.data.get_runs(qid2)[0], run2,
                         "Failed to find the run where input was used.")

        # One cannot remove input that is used by a run
        with self.assertRaises(AscotIOException,
                               msg="Input was incorrectly removed."):
            a5.data["efield"][grp1].destroy(repack=False)

        # Removing one run should preserve results group
        a5.data["q"+run1].destroy(repack=False)
        self.assertFalse(("q"+run1) in a5.data,
                        "Run group not removed in Ascot5IO.")
        with h5py.File(self.testfilename, "r") as h5:
            self.assertTrue("results" in h5,
                            "Results not found in HDF5.")
            self.assertFalse(run1 in h5["results"],
                            "Run not removed in HDF5.")

        # Removing both runs removes the group as well
        a5.data["q"+run2].destroy(repack=False)
        self.assertFalse("results" in a5.data,
                        "Results not removed in Ascot5IO.")
        with h5py.File(self.testfilename, "r") as h5:
            self.assertFalse("results" in h5,
                            "Results not removed in HDF5.")

        # Create a new run and test the remove all results method
        with h5py.File(self.testfilename, "a") as h5:
            group = fileapi.add_group(h5, "results", "run", desc="RUN1")
            group.attrs["qid_efield"] = np.string_(qid1)

        a5.data.destroy(repack=False)
        self.assertFalse("results" in a5.data,
                        "Results not removed in Ascot5IO.")
        with h5py.File(self.testfilename, "r") as h5:
            self.assertFalse("results" in h5,
                            "Results not removed in HDF5.")

    def test_endcond(self):
        """Test end conditions are parsed properly
        """

        state = State(None, None)
        self.assertTrue(state._endcond_check(0x2, "aborted"),
                        "Endcond did not match its binary repr.")
        self.assertFalse(state._endcond_check(0x2, "none"),
                         "Endcond matched wrong binary repr.")
        self.assertFalse(state._endcond_check(0x2, "aborted none"),
                         "Endcond AND operation failed.")
        self.assertFalse(state._endcond_check(0x2, "not aborted"),
                         "Endcond NOT operation failed.")
        self.assertTrue(state._endcond_check(0x2, "aborted not none"),
                        "Endcond AND NOT operation failed.")
        self.assertTrue(state._endcond_check(0x2 | 0x4, "aborted"),
                        "Endcond failed when multiple end conditions active.")
        with self.assertRaises(
                ValueError,
                msg="Failed to raise exception when endcond unknown"):
            state._endcond_check(0x2 | 0x4, "vanished")

    def test_inputs(self):
        inputs = {
            "bfield"  : ["B_TC", "B_GS", "B_2DS", "B_3DS", "B_3DST", "B_STS",],
            "efield"  : ["E_TC", "E_1DS", "E_3D", "E_3DS", "E_3DST",],
            "marker"  : ["Prt", "GC", "FL",],
            "wall"    : ["wall_2D", "wall_3D",],
            "plasma"  : ["plasma_1D", "plasma_1DS",],
            "neutral" : ["N0_3D",],
            "boozer"  : ["Boozer",],
            "mhd"     : ["MHD_STAT", "MHD_NONSTAT",],
            "options" : ["Opt",],
            #"nbi"     : ["NBI",],
            }
        a5 = Ascot(self.testfilename, create=True)
        for parent, inp in inputs.items():
            mod = importlib.import_module("a5py.ascot5io." + parent)

            for grp in inp:
                fun = getattr(getattr(mod, grp), "write_hdf5_dummy")
                fun(self.testfilename)
                a5 = Ascot(self.testfilename)
                data = a5.data[parent].DUMMY.read()

                fun = getattr(getattr(mod, grp), "write_hdf5")
                fun(fn=self.testfilename, desc="DUMMY2", **data)
                a5 = Ascot(self.testfilename)
                a5.data[parent].DUMMY.destroy()
                a5.data[parent].DUMMY2.destroy()

    def tearDown(self):
        """Remove the file used in testing.
        """
        subprocess.run(["rm", "-f", self.testfilename])

if __name__ == '__main__':
    unittest.main()
