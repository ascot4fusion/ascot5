"""Test functionality of a5py.
"""
import numpy as np
import h5py
import time
import unyt
import unittest
import datetime
import subprocess
import importlib
import matplotlib.pyplot as plt

import a5py.routines.plotting as a5plt
from a5py import Ascot, AscotInitException, AscotIOException
from a5py.ascot5io.coreio import fileapi
from a5py.ascot5io.state import State
from a5py.ascot5io.marker import Marker
from a5py.ascot5io import RF2D_fields

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
        dates = a5.data["efield"].get_contents()[1]
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
        grp1 = a5.data.create_input("E_TC", **ETC)
        time.sleep(1.0)
        grp2 = a5.data.create_input("E_TC", **ETC)
        time.sleep(1.0)
        grp3 = a5.data.create_input("E_TC", **ETC)

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

        # Reference by tag: special chars omitted
        a5.data.efield[grp2].set_desc("N1C3T^GBA_T")
        self.assertFalse("TEST" in a5.data.efield,
                         "Reference by tag failed.")
        self.assertEqual(a5.data.efield[grp2], a5.data.efield["N1C3TGBAT"],
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
            group.attrs["qid_efield"] = np.bytes_(qid1)
            time.sleep(1.0)
            group = fileapi.add_group(h5, "results", "run", desc="RUN2")
            group.attrs["qid_efield"] = np.bytes_(qid2)

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
            group.attrs["qid_efield"] = np.bytes_(qid1)

        a5.data.destroy(repack=False)
        self.assertFalse("results" in a5.data,
                        "Results not removed in Ascot5IO.")
        with h5py.File(self.testfilename, "r") as h5:
            self.assertFalse("results" in h5,
                            "Results not removed in HDF5.")

    def test_endcond(self):
        """Test end conditions are parsed properly
        """
        self.assertTrue(State.endcond_check(0x2, "aborted"),
                        "Endcond did not match its binary repr.")
        self.assertFalse(State.endcond_check(0x2, "none"),
                         "Endcond matched wrong binary repr.")
        self.assertFalse(State.endcond_check(0x2, "aborted none"),
                         "Endcond AND operation failed.")
        self.assertFalse(State.endcond_check(0x2, "not aborted"),
                         "Endcond NOT operation failed.")
        self.assertTrue(State.endcond_check(0x2, "aborted not none"),
                        "Endcond AND NOT operation failed.")
        self.assertTrue(State.endcond_check(0x2 | 0x4, "aborted"),
                        "Endcond failed when multiple end conditions active.")
        with self.assertRaises(
                ValueError,
                msg="Failed to raise exception when endcond unknown"):
            State.endcond_check(0x2 | 0x4, "vanished")

    def test_inputs(self):
        inputs = {
            "bfield"  : ["B_TC", "B_GS", "B_2DS", "B_3DS", "B_STS",],
            "efield"  : ["E_TC", "E_1DS"],
                         #"E_3D", "E_3DS", "E_3DST",],
            "marker"  : ["prt", "gc", "fl",],
            "wall"    : ["wall_2D", "wall_3D",],
            "plasma"  : ["plasma_1D", "plasma_1DS", "plasma_1Dt"],
            "neutral" : ["N0_1D", "N0_3D",],
            "boozer"  : ["Boozer",],
            "mhd"     : ["MHD_STAT", "MHD_NONSTAT",],
            "asigma"  : ["asigma_loc"],
            "options" : ["opt",],
            #"nbi"     : ["NBI",],
            }
        a5 = Ascot(self.testfilename, create=True)
        for parent, inp in inputs.items():
            for grp in inp:
                a5.data.create_input(grp, desc="DUMMY")

                a5.data.create_input(grp, desc="DUMMY2")
                data = a5.data[parent].DUMMY2.read()
                a5.data[parent].DUMMY.destroy()
                a5.data[parent].DUMMY2.destroy()

                if parent not in ["marker", "options"]:
                    a5.input_init(**{parent:data})
                    a5.input_free()

    def tearDown(self):
        """Remove the file used in testing.
        """
        subprocess.run(["rm", "-f", self.testfilename],
                       stdout=subprocess.DEVNULL)

class TestAscot(unittest.TestCase):
    """Class for testing :class:`Ascot` object.
    """

    @classmethod
    def setUpClass(cls):
        """Create and run a test case.

        Assumes ascot5_main is located in the same folder.
        """
        super(TestAscot, cls).setUpClass()
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("options tutorial")
        a5.data.create_input("bfield analytical iter circular", splines=False)
        a5.data.create_input("wall rectangular")
        a5.data.create_input("plasma flat")

        from a5py.ascot5io.marker import Marker
        mrk = Marker.generate("gc", n=100, species="alpha")
        mrk["energy"][:] = 3.5e6
        #mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(100,)
        mrk["pitch"][:]  = -0.99 - 0 * np.random.rand(100,)
        mrk["r"][:]      = np.linspace(6.2, 8.2, 100)
        a5.data.create_input("gc", **mrk)
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]))
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]), desc="UNUSED")

        a5.data.create_input("N0_3D")
        a5.data.create_input("Boozer")
        a5.data.create_input("MHD_STAT")
        a5.data.create_input("asigma_loc")

        name = a5.data.options.active.new(
            ENDCOND_MAX_MILEAGE=0.5e-4, DIST_MIN_CHARGE=1.5,
            DIST_MAX_CHARGE=2.5, DIST_NBIN_PPE=10, DIST_NBIN_PPA=20,
            DIST_NBIN_PHI=3, DIST_NBIN_R=10, DIST_NBIN_Z=10, DIST_NBIN_RHO=10,
            desc="Fast")
        a5.data.options[name].activate()

        subprocess.run(["./../../build/ascot5_main", "--in=unittest.h5"],
                       stdout=subprocess.DEVNULL)

    @classmethod
    def tearDownClass(cls):
        super(TestAscot, cls).tearDownClass()
        subprocess.run(["rm", "-f", "unittest.h5"])

    def test_initandsim(self):
        """Test :class:`Ascotpy` initialization and simulation routines.
        """
        a5 = Ascot("unittest.h5")
        a5.data.active.efield.activate()
        allinps = ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd", "asigma"]

        # Test basic init and free
        a5.input_init()
        inps = a5.input_initialized()
        self.assertCountEqual(inps.keys(), allinps,
                              "Not all inputs were initialized.")
        a5.input_free()
        inps = a5.input_initialized()
        self.assertTrue(len(inps) == 0,
                        "Not all inputs were free'd.")

        # Test explicit init and free
        a5.input_init(efield=True)
        inps = a5.input_initialized()
        self.assertCountEqual(inps.keys(), ["efield"],
                              "Explicit initialization failed.")
        a5.input_free(bfield=True)
        inps = a5.input_initialized()
        self.assertCountEqual(inps.keys(), ["efield"],
                              "Explicit free failed.")
        a5.input_free(efield=True)
        inps = a5.input_initialized()
        self.assertTrue(len(inps) == 0,
                        "Explicit free failed.")

        # Test init and free with multiple inputs of same type
        a5.input_init(efield=True)
        inps = a5.input_initialized()
        self.assertEqual(inps["efield"], a5.data.efield.active.get_qid(),
                         "Initialized input was not the active one.")

        try:
            a5.input_init(efield=True)
        except AscotInitException:
            self.fail("Re-initializing same input failed.")

        with self.assertRaises(
                AscotInitException,
                msg="Initializing occupied input did not raise exception."):
            a5.input_init(efield=a5.data.efield["UNUSED"].get_qid())

        with self.assertRaises(
                AscotIOException,
                msg="Initializing nonexistent input did not raise exception."):
            a5.input_init(efield="0123456789")

        try:
            a5.input_init(efield=a5.data.efield["UNUSED"].get_qid(),
                          switch=True)
        except AscotInitException:
            self.fail("Switch did not work.")

        a5.input_free(efield=True)
        inps = a5.input_initialized()
        self.assertTrue(len(inps) == 0,
                        "Not all inputs were free'd.")

        # Double free should be ok
        a5.input_free(bfield=True)
        a5.input_free()

        # Test init from run
        a5.data.efield["UNUSED"].activate()
        a5.input_init(run=True, efield=True)
        inps = a5.input_initialized()
        self.assertEqual(inps["efield"], a5.data.active.efield.get_qid(),
                         "Initialized input was not the correct one.")
        a5.input_free(efield=True)

        # Give run QID, init all inputs
        a5.input_init(run=a5.data.active.get_qid())
        inps = a5.input_initialized()
        self.assertEqual(inps["efield"], a5.data.active.efield.get_qid(),
                         "Initialized input was not the correct one.")
        a5.input_free()
        a5.data.active.efield.activate()

        # Test packing
        a5.simulation_initinputs()
        with self.assertRaises(
                AscotInitException,
                msg="Trying to init packed data should raise an exception."):
            a5.input_init(efield=a5.data.efield["UNUSED"].get_qid())

        # Test marker and options initialization
        opt = a5.data.options.active.read()
        opt.update({"ENDCOND_MAX_MILEAGE" : 1e-10})
        a5.simulation_initoptions(**opt)

        mrk = a5.data.marker.active.read()
        a5.simulation_initmarkers(**mrk)

        # Run simulation and test output initialization
        a5.simulation_run(printsummary=False)

        # Verify that simulations can't be run when data is occupied
        with self.assertRaises(
                AscotInitException,
                msg="Trying to run with occupied output should cause an error"):
            a5.simulation_run(printsummary=False)

        # This should run
        a5.simulation_free(diagnostics=True)
        a5.simulation_run(printsummary=False)
        a5.simulation_free(diagnostics=True)

        # Can't run when input is missing
        a5.simulation_free(markers=True, diagnostics=True)
        with self.assertRaises(
                AscotInitException,
                msg="Trying to run without markers should cause an error"):
            a5.simulation_run(printsummary=False)

        # Init markers but free inputs
        a5.simulation_initmarkers(**mrk)
        a5.simulation_free(inputs=True)
        with self.assertRaises(
                AscotInitException,
                msg="Trying to run without input should cause an error"):
            a5.simulation_run(printsummary=False)

        # Free everyting
        a5.simulation_free()
        with self.assertRaises(
                AscotInitException,
                msg="Trying to run without input should cause an error"):
            a5.simulation_run(printsummary=False)

    def test_eval(self):
        """Test :class:`Ascotpy` evaluation routines.
        """
        a5 = Ascot("unittest.h5")

        # Input evaluations
        a5.input_init()
        a5.input_init(rffield=RF2D_fields.create_test_fields())
        inputqnts = a5.input_eval_list(show=False)
        r = 6.2 * unyt.m; phi = 0 * unyt.deg; z = 0 * unyt.m; t = 0 * unyt.s
        br, bphi, bz = a5.input_eval(r, phi, z, t, "br", "bphi", "bz")
        for q in inputqnts:
            out = a5.input_eval(r, phi, z, t, q)
        a5.input_free()

        # State evaluations
        a5.input_init(bfield=True)
        e, v, m = a5.data.active.getstate(
            "ekin", "vpar", "mu", endcond="tlim", ids=1)
        self.assertTrue(
            e.units == unyt.eV and v.units == unyt.m/unyt.s and \
            m.units == unyt.eV/unyt.T and e.size==1 and v.size==1 and \
            m.size==1,
            "State evaluation failed")
        outputqnts = a5.data.active.getstate_list()
        for q in outputqnts:
            out = a5.data.active.getstate(q)

        # Orbit evaluations
        e, v, m = a5.data.active.getorbit(
            "ekin", "vpar", "mu", endcond="tlim", ids=1)
        self.assertTrue(
            e.units == unyt.eV and v.units == unyt.m/unyt.s and \
            m.units == unyt.eV/unyt.T,
            "Orbit evaluation failed")
        outputqnts = a5.data.active.getorbit_list()
        for q in outputqnts:
            out = a5.data.active.getorbit(q)
        a5.input_free()

        # Packing should not prevent input evaluation
        a5.simulation_initinputs()
        out = a5.input_eval(r, phi, z, t, "bnorm")

        # Prepare simulation data
        opt = a5.data.options.active.read()
        opt.update({"ENDCOND_MAX_MILEAGE" : 1e-7})
        a5.simulation_initoptions(**opt)
        mrk = a5.data.marker.active.read()
        a5.simulation_initmarkers(**mrk)
        vrun = a5.simulation_run(printsummary=False)

        # Virtual state evaluations
        e, v, m = vrun.getstate(
            "ekin", "vpar", "mu", endcond="tlim", ids=1)
        self.assertTrue(
            e.units == unyt.eV and v.units == unyt.m/unyt.s and \
            m.units == unyt.eV/unyt.T and e.size==1 and v.size==1 and \
            m.size==1,
            "Virtual state evaluation failed")
        outputqnts = a5.data.active.getstate_list()
        for q in outputqnts:
            out = vrun.getstate(q)

        # Virtual orbit evaluations
        out = vrun.getorbit("r", "phi", "z")
        e, v, m = vrun.getorbit(
            "ekin", "vpar", "mu", endcond="tlim", ids=1)
        self.assertTrue(
            e.units == unyt.eV and v.units == unyt.m/unyt.s and \
            m.units == unyt.eV/unyt.T,
            "Virtual orbit evaluation failed")
        outputqnts = vrun.getorbit_list()
        for q in outputqnts:
            out = a5.data.active.getorbit(q)
        a5.simulation_free()

    def test_dist(self):
        """Test distribution postprocessing.
        """
        a5 = Ascot("unittest.h5")

        # Test volume calculation
        R,Z = np.meshgrid(np.linspace(3,7,101), np.linspace(-2,2,101))
        psi = np.sqrt((R-5)**2 + (Z-0)**2)
        bdata = {"rmin":3, "rmax":7, "nr":101, "zmin":-2, "zmax":2, "nz":101,
                 "axisr":5, "axisz":0, "psi":psi, "psi0":0, "psi1":1,
                 "br":psi*0, "bphi":psi*0, "bz":psi*0}
        qid = a5.data.create_input("B_2DS", **bdata).split("_")[-1]

        try:
            a5.input_init(bfield=qid)
            vol1, area1, r1, p1, z1 = a5.input_rhovolume(
                method="prism", nrho=5, ntheta=4, nphi=3, return_area=True,
                return_coords=True)
            vol2, area2, r2, p2, z2 = a5.input_rhovolume(
                method="mc", nrho=5, ntheta=4, nphi=3, return_area=True,
                return_coords=True)
            a5.input_free()
        except:
            a5.input_free()
            raise

        print(r1.shape, r2.shape)
        print(p1.shape, p2.shape)
        print(z1.shape, z2.shape)
        delta = 1
        vol0 = 2*np.pi*5*np.pi*1.0**2
        print(
            """
            Volume calculation:
            %3.3f %3.3f %3.3f
            """ % (np.sum(vol1).v, np.sum(vol2).v, vol0)
        )
        self.assertAlmostEqual(
            np.sum(vol1).v, vol0, None, "Volume calculation with prism failed",
            delta)
        self.assertAlmostEqual(
            np.sum(vol2).v, vol0, None, "Volume calculation with MC failed",
            delta)
        area0 = np.pi*1.0**2
        print(
            """
            Area calculation:
            %3.3f %3.3f %3.3f
            """ % (np.sum(area1[:,:,0]).v, np.sum(area2[:,:,0]).v, area0)
        )
        self.assertAlmostEqual(
            np.sum(area1[:,:,0]).v, area0, None,
            "Area calculation with prism failed", delta)
        self.assertAlmostEqual(
            np.sum(area2[:,:,0]).v, area0, None,
            "Area calculation with MC failed", delta)

        # Test distribution initializing
        a5.data.active.getdist("5d")
        a5.data.active.getdist("rho5d")
        a5.data.active.getdist("6d")
        a5.data.active.getdist("rho6d")
        a5.data.active.getdist("com")

        # Test slicing, interpolating, and integrating
        dist = a5.data.active.getdist("5d")
        dist.integrate(charge=np.s_[:], time=np.s_[:])
        self.assertFalse("charge" in dist.abscissae, "Integration failed")
        self.assertFalse("time"   in dist.abscissae, "Integration failed")

        dist.slice(ppar=np.s_[0], pperp=np.s_[0])
        self.assertTrue(dist.abscissa("ppar").size == 1, "Slicing failed")
        self.assertTrue(dist.abscissa("pperp").size == 1, "Slicing failed")

        out = dist.interpolate(phi=2*unyt.deg, r=6.2*unyt.m)
        self.assertEqual(out.size, dist.abscissa("z").size,
                         "Integration/slicing/interpolation failed")

        # Get 5D distribution in both momentum spaces and verify the number of
        # particles is same
        dist    = a5.data.active.getdist("5d")
        exidist = a5.data.active.getdist("5d", exi=True, ekin_edges=10,
                                         pitch_edges=20, plotexi=False)
        time, weight = a5.data.active.getstate("mileage", "weight", state="end")
        self.assertAlmostEqual(np.sum(dist.histogram().v),
                               np.sum(time*weight).v, None,
                               "Distribution not valid",
                               0.001)
        self.assertAlmostEqual(np.sum(exidist.histogram().v),
                               np.sum(time*weight).v, None,
                               "Converting ppar-pperp to ekin-pitch failed",
                               0.001)

        a5.input_init(bfield=True)
        rhodist    = a5.data.active.getdist("rho5d")
        rhoexidist = a5.data.active.getdist("rho5d", exi=True)
        mom1 = a5.data.active.getdist_moments(dist, "density")
        mom2 = a5.data.active.getdist_moments(exidist, "density")
        mom3 = a5.data.active.getdist_moments(rhodist, "density")
        mom4 = a5.data.active.getdist_moments(rhoexidist, "density")
        a5.input_free()
        print(np.sum(mom1.ordinate("density") * mom1.volume),
              np.sum(mom2.ordinate("density") * mom2.volume),
              np.sum(mom3.ordinate("density") * mom3.volume),
              np.sum(mom4.ordinate("density") * mom4.volume))

        self.assertAlmostEqual(np.sum(mom1.ordinate("density") * mom1.volume).v,
                               np.sum(time*weight).v, None,
                               "Converting ppar-pperp to ekin-pitch failed",
                               0.001)
        self.assertAlmostEqual(np.sum(mom2.ordinate("density") * mom2.volume).v,
                               np.sum(time*weight).v, None,
                               "Converting ppar-pperp to ekin-pitch failed",
                               0.001)
        self.assertAlmostEqual(np.sum(mom3.ordinate("density") * mom3.volume).v,
                               np.sum(time*weight).v, None,
                               "Converting ppar-pperp to ekin-pitch failed",
                               0.001)
        self.assertAlmostEqual(np.sum(mom4.ordinate("density") * mom4.volume).v,
                               np.sum(time*weight).v, None,
                               "Converting ppar-pperp to ekin-pitch failed",
                               0.001)

class TestMoments(unittest.TestCase):
    """Class for testing distribution moments.
    """

    @classmethod
    def setUpClass(cls):
        """Create and run a test case.

        Assumes ascot5_main is located in the same folder.
        """
        super(TestMoments, cls).setUpClass()
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("options tutorial")
        a5.data.create_input("bfield analytical iter circular")
        a5.data.create_input("wall rectangular")
        a5.data.create_input("plasma flat")

        from a5py.ascot5io.marker import Marker
        mrk = Marker.generate("gc", n=1, species="alpha")
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.7
        mrk["r"][:]      = 7.8
        a5.data.create_input("gc", **mrk)
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]))

        a5.data.create_input("N0_3D")
        a5.data.create_input("Boozer")
        a5.data.create_input("MHD_STAT")
        a5.data.create_input("asigma_loc")

        name = a5.data.options.active.new(
            ENDCOND_MAX_MILEAGE=1.5e-3, DIST_NBIN_TIME=1, DIST_MIN_CHARGE=1.5,
            DIST_MAX_CHARGE=2.5, DIST_NBIN_CHARGE=1, DIST_MIN_R=4, DIST_MAX_R=8,
            DIST_MIN_Z=-2, DIST_MAX_Z=2,
            ENABLE_DIST_6D=0, ENABLE_DIST_RHO6D=1, ENABLE_DIST_COM=0,
            DIST_NBIN_PHI=1, DIST_NBIN_R=50, DIST_NBIN_Z=100,
            DIST_NBIN_PPE=50, DIST_NBIN_PPA=50,
            DIST_NBIN_THETA=2,  DIST_NBIN_RHO=10, ENABLE_COULOMB_COLLISIONS=1,
            ORBITWRITE_NPOINT=60000)
        a5.data.options[name].activate()

        subprocess.run(["./../../build/ascot5_main", "--in=unittest.h5"],
                       stdout=subprocess.DEVNULL)

    @classmethod
    def tearDownClass(cls):
        super(TestMoments, cls).tearDownClass()
        subprocess.run(["rm", "-f", "unittest.h5"])

    def test_moments(self):
        a5 = Ascot("unittest.h5")
        ordinates = ["density", "chargedensity", "energydensity", "pressure",
                     "toroidalcurrent", "parallelcurrent", "powerdep",
                     "electronpowerdep", "ionpowerdep"]
        #ordinates = ["jxbtorque"]

        a5.input_init(bfield=True, plasma=True)
        dist    = a5.data.active.getdist("5d")
        mom     = a5.data.active.getdist_moments(dist, *ordinates)
        rhodist = a5.data.active.getdist("rho5d")
        rhomom  = a5.data.active.getdist_moments(rhodist, *ordinates)

        re = dist.abscissa_edges("r")
        ze = dist.abscissa_edges("z")
        v1, v2 = np.meshgrid(re[1:]**2 - re[:-1]**2, ze[1:] - ze[:-1])
        vol = np.pi * v1 * v2

        rho, r, z, phi, psi, weight, time,  \
        charge, energy, vphi, vnorm, vpar,  \
        mass, p, pitch, bnorm, bphi, ppar = \
            a5.data.active.getorbit("rho", "r", "z", "phi", "psi", "weight",
                                    "mileage", "charge", "ekin", "vphi",
                                    "vnorm", "vpar", "mass", "pnorm", "pitch",
                                    "bnorm", "bphi", "ppar")
        ei, psii, Pphii = a5.data.active.getstate("ekin", "psi", "pphi",
                                                  state="ini")
        ef, tf = a5.data.active.getstate("ekin", "mileage", state="end")
        dt = np.diff(time, prepend=0)

        k, nu = a5.input_eval_collcoefs(
            mass[0], charge[0], r, phi, z, time, vnorm, "k", "nu",
            grid=False)
        dEtot_d = p * dt * np.sum(k, axis=0)
        dEele_d = p * dt * k[0,:]
        dEion_d = p * dt * np.sum(k[1:,:], axis=0)
        dpsi    = np.zeros(psi.shape) * psi.units
        dpsi[0] = psi[0] - psii
        dpsi[1:] = psi[1:] - psi[:-1]
        nu      = np.sum(nu, axis=0)
        k       = np.sum(k, axis=0)
        dppar   = mass * k * pitch - p * pitch * nu
        Pphi    = ppar * r * (bphi/bnorm) + charge * psi
        dPphi   = np.diff(Pphi, prepend=Pphi[0])

        a5.input_free()

        mrkdist = {}
        mrkdist["density"]          = weight * dt
        mrkdist["chargedensity"]    = weight * charge * dt
        mrkdist["energydensity"]    = weight * (energy * dt).to("J*s")
        mrkdist["pressure"]         = weight * (mass * vnorm**2*dt).to("J*s")/3
        mrkdist["toroidalcurrent"]  = weight * ( charge * vphi * dt ).to("A*s*m")
        mrkdist["parallelcurrent"]  = weight * ( charge * vpar * dt ).to("A*s*m")
        mrkdist["powerdep"]         = weight * dEtot_d.to("J")
        mrkdist["electronpowerdep"] = weight * dEele_d.to("J")
        mrkdist["ionpowerdep"]      = weight * dEion_d.to("J")
        mrkdist["jxbtorque"]        = weight * (-charge * dpsi/unyt.s).to("N*m")
        mrkdist["colltorque"]       = weight * (r*dppar*(bphi/bnorm)*dt).to("J*s")
        mrkdist["canmomtorque"]     = weight * -charge * dPphi

        print(weight[0]*tf)
        print(((ef-ei)*weight[0]).to("W"))
        for o in ordinates:
            a1 = np.sum(rhomom.ordinate(o) * rhomom.volume)
            a2 = np.sum(mom.ordinate(o) * mom.volume)
            a3 = np.sum(mrkdist[o])
            print(o, a1, a2, a3)

class TestPlotting(unittest.TestCase):

    def tearDown(self):
        try:
            a5 = Ascot("unittest.h5")
            a5.simulation_free()
        except:
            pass
        subprocess.run(["rm", "-f", "unittest.h5"])

    def initandrunsim(self, simtime):
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("options tutorial")
        a5.data.create_input("bfield analytical iter circular")
        a5.data.create_input("wall rectangular")
        a5.data.create_input("plasma flat")
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]))
        a5.data.create_input("N0_3D")
        a5.data.create_input("Boozer")
        a5.data.create_input("MHD_STAT")
        a5.data.create_input("asigma_loc")

        a5.simulation_initinputs()

        mrk = Marker.generate("gc", n=100, species="alpha")
        mrk["weight"][:] = 10**(10 +10*np.random.rand(100,))
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(100,)
        mrk["r"][:]      = np.linspace(6.2, 8.5, 100)
        a5.simulation_initmarkers(**mrk)

        opt = a5.data.options.active.read()
        opt.update(ENDCOND_MAX_MILEAGE=simtime)
        a5.simulation_initoptions(**opt)
        return a5.simulation_run(printsummary=False)

    def test_input_plotrz(self):
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("bfield analytical iter circular", splines=True,
                             axisymmetric=True)
        a5.data.create_input("plasma flat")
        a5.input_init(bfield=True, plasma=True)

        fig = plt.figure(figsize=(30,40))

        r = np.linspace( 3.0, 9.0,  50) * unyt.m
        z = np.linspace(-2.0, 2.0, 100) * unyt.m

        ax = fig.add_subplot(2,4,1)
        a5.input_plotrz(r, z, "bnorm", clim=[None, None], cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,5)
        a5.input_plotrz(r, z, "bnorm", clim=[0, 10], cmap="Greys", axes=ax)

        ax = fig.add_subplot(2,4,2)
        a5.input_plotrz(r, z, "br", clim=[None, None], cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,6)
        a5.input_plotrz(r, z, "br", clim=[-1, 1], cmap="PuOr", axes=ax)

        ax = fig.add_subplot(2,4,3)
        a5.input_plotrz(r, z, "log ne", clim=[None, None], cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,7)
        a5.input_plotrz(r, z, "log ne", clim=[1e10, 1e22],
                        cmap="Greys", axes=ax)

        ax = fig.add_subplot(2,4,4)
        a5.input_plotrz(r, z, "log divb", clim=[None, None],
                        cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,8)
        a5.input_plotrz(r, z, "log divb",clim=[-1e-14, 1e-14],
                        cmap="PuOr",axes=ax)

        a5.input_free()
        plt.show(block=False)

    def test_plotstate_scatter(self):
        vr = self.initandrunsim(1e-3)
        fig = plt.figure(figsize=(30,40))

        ax = fig.add_subplot(2,3,1)
        vr.plotstate_scatter("end rho", "log end ekin",
                             axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,2)
        vr.plotstate_scatter("ini rho", "log diff ekin",
                             axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,3)
        vr.plotstate_scatter("log end x", "log end y", c="log diff phi",
                             axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,4, projection="3d")
        vr.plotstate_scatter("end x", "end y", "end z",
                             axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,5, projection="3d")
        vr.plotstate_scatter("end x", "end y", "end z", c="end rho",
                             axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,6, projection="3d")
        vr.plotstate_scatter("log end x", "log end y", "log end z",
                             c="log diff phi", axesequal=True, axes=ax)

        plt.show(block=False)

    def test_plotstate_histogram(self):
        vr = self.initandrunsim(1e-4)
        fig = plt.figure(figsize=(30,40))

        ax = fig.add_subplot(2,3,1)
        vr.plotstate_histogram("end phimod", axes=ax)

        ax = fig.add_subplot(2,3,2)
        vr.plotstate_histogram("log end phi", xbins=10, endcond="WALL", axes=ax)

        ax = fig.add_subplot(2,3,3)
        vr.plotstate_histogram("end ekin", xbins=np.linspace(0,4e6,100),
                               logscale=True, weight=True, axes=ax)

        ax = fig.add_subplot(2,3,4)
        vr.plotstate_histogram("end phimod", "end thetamod", axes=ax)

        ax = fig.add_subplot(2,3,5)
        vr.plotstate_histogram("log end phi", "log end theta", endcond="WALL",
                               axes=ax)

        ax = fig.add_subplot(2,3,6)
        vr.plotstate_histogram("end phimod", "end thetamod", logscale=True,
                               weight=True, axes=ax)

        plt.show(block=False)

    def test_plotorbit_trajectory(self):
        vr = self.initandrunsim(1e-4)
        fig = plt.figure(figsize=(30,40))

        ax = fig.add_subplot(2,3,1)
        vr.plotorbit_trajectory("r", "z",
                                axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,2)
        vr.plotorbit_trajectory("mileage", "reldiff ekin", c="rho",
                                axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,3)
        vr.plotorbit_trajectory("log mileage", "log reldiff ekin",
                                c="log rho",
                                axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,4, projection="3d")
        vr.plotorbit_trajectory("x", "y", "z",
                                axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,5, projection="3d")
        vr.plotorbit_trajectory("x", "y", "z", c="diff ekin",
                                axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,6, projection="3d")
        vr.plotorbit_trajectory("log x", "log y", "log z", c="log diff ekin",
                                axesequal=True, axes=ax)

        plt.show(block=False)

if __name__ == '__main__':
    unittest.main()
