"""Unit tests for most post-processing routines.
"""
import numpy as np
import unittest
import subprocess
import unyt

from a5py import Ascot, AscotInitException, AscotIOException

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
        a5.data.create_input("bfield analytical iter circular")
        a5.data.create_input("wall rectangular")
        a5.data.create_input("plasma flat")

        from a5py.ascot5io.marker import Marker
        mrk = Marker.generate("gc", n=100, species="alpha")
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(100,)
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
            DIST_MAX_CHARGE=2.5, DIST_NBIN_PPE=50, DIST_NBIN_PPA=100,
            desc="Fast")
        a5.data.options[name].activate()

        subprocess.run(["./ascot5_main", "--in=unittest.h5"])

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

        # This hould run
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
            vol1, area1 = a5.input_rhovolume(
                method="prism", nrho=20, ntheta=20, nphi=2, return_area=True)
            vol2, area2 = a5.input_rhovolume(
                method="mc", nrho=2, ntheta=5, nphi=2, return_area=True)
            a5.input_free()
        except:
            a5.input_free()
            raise

        delta = 1
        vol0 = 2*np.pi*5*np.pi*1.0**2
        self.assertAlmostEqual(
            np.sum(vol1).v, vol0, None, "Volume calculation with prism failed",
            delta)
        self.assertAlmostEqual(
            np.sum(vol2).v, vol0, None, "Volume calculation with MC failed",
            delta)
        area0 = np.pi*1.0**2
        self.assertAlmostEqual(
            np.sum(area1).v, area0, None, "Area calculation with prism failed",
            delta)
        self.assertAlmostEqual(
            np.sum(area2).v, area0, None, "Area calculation with MC failed",
            delta)

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
        self.assertFalse("time" in dist.abscissae, "Integration failed")

        dist.slice(ppar=np.s_[0], pperp=np.s_[0])
        self.assertTrue(dist.abscissa("ppar").size == 1, "Slicing failed")
        self.assertTrue(dist.abscissa("pperp").size == 1, "Slicing failed")

        out = dist.interpolate(phi=2*unyt.deg, r=6.2*unyt.m)
        self.assertEqual(out.size, dist.abscissa("z").size,
                         "Integration/slicing/interpolation failed")

        # Get 5D distribution in both momentum spaces and verify the number of
        # particles is same
        dist    = a5.data.active.getdist("5d")
        exidist = a5.data.active.getdist("5d", exi=True)
        time, weight = a5.data.active.getstate("mileage", "weight", state="end")
        #print(np.sum(dist.v), np.sum(exidist.v), np.sum(time*weight).v)
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
        #print(np.sum(mom1.ordinate("density") * mom1.volume),
        #      np.sum(mom2.ordinate("density") * mom2.volume),
        #      np.sum(mom3.ordinate("density") * mom3.volume),
        #      np.sum(mom4.ordinate("density") * mom4.volume))

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
if __name__ == '__main__':
    unittest.main()
