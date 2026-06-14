"""Module for verifying that ASCOT5 can model the physics it claims.

Implemented tests:

  - elementary: verify that ASCOT5 reproduces correct gyromotion and drifts
    in a simple electromagnetic field.
  - orbitfollowing: verify the conservation properties of orbit-integrators.
  - gctransform: verify the guiding center transformation implementation.
  - ccoll: verify the Coulomb collision operator.
  - classical: verify that ASCOT5 reproduces classical transport correctly.
  - neoclassical: verify that ASCOT5 reproduces the neoclassical transport
    correctly
  - boozer: verify Boozer coordinate transformation.
  - mhd: verify inclusion of MHD modes.
  - atomic: verify implementation of ionization and neutralization reactions.
"""
import copy
import unyt
import subprocess
import warnings
import scipy
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec

from a5py import Ascot, physlib
from a5py.plotting import plotting as a5plt
from a5py.ascot5io.options import Opt
from a5py.ascot5io.marker import Marker
from a5py.ascot5io.bfield import B_2DS
from a5py.ascot5io.nbi import Injector

class PhysTest():

    tag_elementary_gyro    = "TESTELEMENTARYGYRO"
    tag_elementary_exbgo   = "TESTELEMENTARYEXBGO"
    tag_elementary_exbgc   = "TESTELEMENTARYEXBGC"
    tag_elementary_gradbgo = "TESTELEMENTARYGRADBGO"
    tag_elementary_gradbgc = "TESTELEMENTARYGRADBGC"
    tag_orbfol_go          = "TESTORBFOLGO"
    tag_orbfol_gcf         = "TESTORBFOLGCF"
    tag_orbfol_gca         = "TESTORBFOLGCA"
    tag_gctransform_go     = "TESTGCTRANSFORMGO"
    tag_gctransform_gc     = "TESTGCTRANSFORMGC"
    tag_gctransform_go2gc  = "TESTGCTRANSFORMGO2GC"
    tag_gctransform_zeroth = "TESTGCTRANSFORMZEROTH"
    tag_gctransform_first  = "TESTGCTRANSFORMFIRST"
    tag_ccoll_thermalgo    = "TESTCCOLLTHERMALGO"
    tag_ccoll_thermalgcf   = "TESTCCOLLTHERMALGCF"
    tag_ccoll_thermalgca   = "TESTCCOLLTHERMALGCA"
    tag_ccoll_slowinggo    = "TESTCCOLLSLOWINGGO"
    tag_ccoll_slowinggcf   = "TESTCCOLLSLOWINGGCF"
    tag_ccoll_slowinggca   = "TESTCCOLLSLOWINGGCA"
    tag_classical_go       = "TESTCLASSGO"
    tag_classical_gcf      = "TESTCLASSGCF"
    tag_classical_gca      = "TESTCLASSGCA"
    tag_neoclassical_go    = "TESTNEOCLASSGO"
    tag_neoclassical_gcf   = "TESTNEOCLASSGCF"
    tag_neoclassical_gca   = "TESTNEOCLASSGCA"
    tag_boozer             = "TESTBOOZER"
    tag_mhd_go             = "TESTMHDGO"
    tag_mhd_gcf            = "TESTMHDGCF"
    tag_mhd_gca            = "TESTMHDGCA"
    tag_atomic_cx          = "TESTATOMICCX"
    tag_atomic_ionz        = "TESTATOMICIONZ"
    tag_afsi_thermal       = "TESTAFSITHERMAL"
    tag_afsi_beamthermal   = "TESTAFSIBEAMTHERMAL"
    tag_afsi_beambeam      = "TESTAFSIBEAMBEAM"
    tag_bbnbi              = "TESTBBNBI"
    tag_biosaw             = "TESTBIOSAW"

    def __init__(self, fn="testascot.h5"):
        try:
            self.ascot = Ascot(fn)
        except FileNotFoundError:
            self.ascot = Ascot(fn, create=True, mute="no")
            init = self.ascot.data.create_input
            init("opt",        desc="DUMMY")
            init("gc",         desc="DUMMY")
            init("B_TC",       desc="DUMMY")
            init("E_TC",       desc="DUMMY")
            init("wall_2D",    desc="DUMMY")
            init("plasma_1D",  desc="DUMMY")
            init("N0_1D",      desc="DUMMY")
            init("Boozer",     desc="DUMMY")
            init("MHD_STAT",   desc="DUMMY")
            init("asigma_loc", desc="DUMMY")

    def execute(self, init=True, run=True, check=True, tests=None):
        """Execute test(s).

        Parameters
        ----------
        """
        if tests is None:
            tests = ["elementary", "orbitfollowing", "gctransform", "ccoll",
                     "classical", "neoclassical", "boozer", "mhd", "atomic",
                     "afsi", "bbnbi", "biosaw"]
        elif isinstance(tests, str):
            tests = [tests]

        import time

        failed = False
        duration = 0
        for test in tests:
            if init:
                getattr(self, "init_" + test)()
                print("Test %s initialized" % test)
            if run:
                start = time.time()
                getattr(self, "run_" + test)()
                print("Test %s simulation complete" % test)
                duration = time.time() - start
            if check:
                a5plt.setpaperstyle()
                passed = getattr(self, "check_" + test)()
                if passed:
                    print("Test %s check passed in %d s" % (test,duration))
                else:
                    print("Test %s check FAILED in %d s" % (test,duration))
                    failed = True
        return failed


    def _activateinputs(self, tag):
        data = self.ascot.data
        tag0 = tag if hasattr(data.bfield, tag) else "DUMMY"
        data.bfield[tag0].activate()
        tag0 = tag if hasattr(data.efield, tag) else "DUMMY"
        data.efield[tag0].activate()
        tag0 = tag if hasattr(data.wall, tag) else "DUMMY"
        data.wall[tag0].activate()
        tag0 = tag if hasattr(data.plasma, tag) else "DUMMY"
        data.plasma[tag0].activate()
        tag0 = tag if hasattr(data.neutral, tag) else "DUMMY"
        data.neutral[tag0].activate()
        tag0 = tag if hasattr(data.boozer, tag) else "DUMMY"
        data.boozer[tag0].activate()
        tag0 = tag if hasattr(data.mhd, tag) else "DUMMY"
        data.mhd[tag0].activate()
        tag0 = tag if hasattr(data.asigma, tag) else "DUMMY"
        data.asigma[tag0].activate()
        tag0 = tag if hasattr(data.options, tag) else "DUMMY"
        data.options[tag0].activate()
        tag0 = tag if hasattr(data.marker, tag) else "DUMMY"
        data.marker[tag0].activate()

    def _runascot(self, test):
        subprocess.call(
            ["./../../build/ascot5_main", "--in=testascot.h5", "--d="+test],
            stdout=subprocess.DEVNULL)
        self.ascot = Ascot(self.ascot.file_getpath())

def _linfit(x, y):
    """Calculate linear fit and return the slope.
    """
    try:
        return np.polyfit(x, y, 1)[0]
    except TypeError:
        return np.polyfit(x, y.v, 1)[0]

if __name__ == '__main__':
    test = PhysTest()
    # Atomic does not work standalone yet
    failed = test.execute(
        init=True, run=True, check=True,
        tests=["elementary", "orbitfollowing", "gctransform", "ccoll",
               "classical", "neoclassical", "boozer", "mhd", "afsi",
               "bbnbi", "biosaw", "atomic"])
    plt.show(block=False)
    if failed: raise Exception("Verification failed")
