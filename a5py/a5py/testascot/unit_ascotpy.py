"""Test `Ascot` object.
"""
import numpy as np
import unittest
import subprocess

from a5py import Ascot, AscotInitException, AscotIOException

class TestAscot(unittest.TestCase):
    """Class for testing `Ascot` object.
    """

    @classmethod
    def setUpClass(cls):
        super(TestAscot, cls).setUpClass()
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_premade(predef="analytical")
        a5.data.create_input("N0_3D", None)
        a5.data.create_input("Boozer", None)
        a5.data.create_input("MHD_STAT", None)
        a5.data.create_input("E_TC", {"exyz" : np.array([0,0,0]),
                                      "desc" : "UNUSED"})

        subprocess.run(["./ascot5_main", "--in=unittest.h5"])

    @classmethod
    def tearDownClass(cls):
        super(TestAscot, cls).tearDownClass()
        #subprocess.run(["rm", "-f", "unittest.h5"])

    def test_initandsim(self):
        """Test `ascotpy` initialization and simulation routines.
        """
        a5 = Ascot("unittest.h5")
        a5.data.active.efield.activate()
        allinps = ["bfield", "efield", "plasma", "wall", "neutral", "boozer",
                    "mhd"]

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

        a5.input_init(run=a5.data.active.get_qid())
        inps = a5.input_initialized()
        self.assertEqual(inps["efield"], a5.data.active.efield.get_qid(),
                         "Initialized input was not the correct one.")
        a5.input_free()
        a5.data.active.efield.activate()

        # Test packing
        a5.simulation_initinput()

        with self.assertRaises(
                AscotInitException,
                msg="Trying to init packed data should raise an exception."):
            a5.input_init(efield=a5.data.efield["UNUSED"].get_qid())

        # Test marker and options initialization
        opt = a5.data.options.active.read()
        opt.update({"ENDCOND_MAX_MILEAGE" : 1e-10})
        a5.simulation_initoptions(opt)

        mrk = a5.data.marker.active.read()
        a5.simulation_initmarkers(mrk)

        # Run simulation and test output initialization
        a5.simulation_run(printsummary=False)

        # Verify that simulations can't be run when data is occupied or missing
        a5.simulation_free(inputs=True, markers=True, diagnostics=True)

    def test_eval(self):
        """Test `ascotpy` evaluation routines.
        """
        a5 = Ascot("unittest.h5")

        # Input evaluations
        a5.input_init()
        inputqnts = a5.input_eval_list()
        for q in inputqnts:
            out = a5.input_eval(6.2, 0.0, 0.0, 0.0, q)
            #print(out)
        a5.input_free()

        # State evaluations
        a5.input_init(bfield=True)
        out = a5.data.active.getstate("psi")
        print(out)
        a5.input_free()

        # Orbit evaluations

        # Simulation evaluations
        a5.simulation_initinput()
        for q in ["bphi"]:
            out = a5.input_eval(6.2, 0.0, 0.0, 0.0, q)
            #print(out)

        opt = a5.data.options.active.read()
        opt.update({"ENDCOND_MAX_MILEAGE" : 1e-7})
        a5.simulation_initoptions(opt)

        mrk = a5.data.marker.active.read()
        a5.simulation_initmarkers(mrk)

        a5.simulation_run(printsummary=False)

        out = a5.simulation_getstate("mileage")
        #print(out)

        out = a5.simulation_getorbit("r",ids=2)
        #print(out)

    def test_postpro(self):
        """Test all postprocessing routines.
        """
        pass

if __name__ == '__main__':
    # Simple simulation with results is required for the tests.
    unittest.main()
