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
from a5py.routines import plotting as a5plt
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
            self.ascot = Ascot(fn, create=True)
            init = self.ascot.data.create_input
            init("opt")
            init("gc")
            init("B_TC")
            init("E_TC")
            init("wall_2D")
            init("plasma_1D")
            init("N0_1D")
            init("Boozer")
            init("MHD_STAT")
            init("asigma_loc")

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

    def init_elementary(self):
        """Initialize data for the elementary test.

        This test covers elementary results such as the gyromotion and drifts by
        verifying that ASCOT5 reproduces:

        1. correct gyroradius and gyrofrequency in an uniform magnetic field.
        2. correct E X B drift in uniform electromagnetic field.
        3. correct gradB drift in a magnetic field with constant gradient.

        Test 1 is done with GO mode and 2 and 3 with both GO and GC modes,
        latter using the fixed step scheme. Tests are done in Cartesian
        electromagnetic field (B_TC and E_TC) and without collisions.
        The test particles are an energetic electron and a positron so
        the tests also verify that ASCOT5 is valid in the relativistic regime.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_elementary_gyro):
            warnings.warn("Inputs already present: Test elementary")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_SIMTIMELIM" : 1,
            "ENDCOND_LIM_SIMTIME" : 2e-9, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
            "ORBITWRITE_INTERVAL" : 1e-11, "ORBITWRITE_NPOINT" : 202
        })
        init("opt", **opt, desc=PhysTest.tag_elementary_gyro)
        opt.update({
            "FIXEDSTEP_USERDEFINED" : 1e-10, "ENDCOND_LIM_SIMTIME" : 1e-7,
            "ORBITWRITE_NPOINT" : 10002
        })
        init("opt", **opt, desc=PhysTest.tag_elementary_exbgo)
        init("opt", **opt, desc=PhysTest.tag_elementary_gradbgo)
        opt.update({
            "SIM_MODE" : 2, "FIXEDSTEP_USERDEFINED" : 1e-9,
            "ORBITWRITE_INTERVAL" : 1e-9, "ORBITWRITE_NPOINT" : 102
        })
        init("opt", **opt, desc=PhysTest.tag_elementary_exbgc)
        init("opt", **opt, desc=PhysTest.tag_elementary_gradbgc)

        # Magnetic field
        d = {"bxyz" : np.array([5, 0, 0]),
             "jacobian" : np.array([0,0,0,0,0,0,0,0,0]),
             "rhoval" : 1.5}
        init("B_TC", **d, desc=PhysTest.tag_elementary_gyro)
        init("B_TC", **d, desc=PhysTest.tag_elementary_exbgo)
        init("B_TC", **d, desc=PhysTest.tag_elementary_exbgc)

        d.update({"jacobian" : np.array([0,0,0.1,0,0,0,0,0,0])})
        init("B_TC", **d, desc=PhysTest.tag_elementary_gradbgo)
        init("B_TC", **d, desc=PhysTest.tag_elementary_gradbgc)

        # Electric field
        d = {"exyz" : np.array([0, 0, 0])}
        init("E_TC", **d, desc=PhysTest.tag_elementary_gyro)
        init("E_TC", **d, desc=PhysTest.tag_elementary_gradbgo)
        init("E_TC", **d, desc=PhysTest.tag_elementary_gradbgc)

        d.update({"exyz" : np.array([0, 1e6, 0])})
        init("E_TC", **d, desc=PhysTest.tag_elementary_exbgo)
        init("E_TC", **d, desc=PhysTest.tag_elementary_exbgc)

        # Marker input is an electron and positron
        mrk = Marker.generate("gc", n=2, species="electron")
        mrk["charge"]    = np.array([1, -1])
        mrk["r"][:]      = 5
        mrk["phi"][:]    = 90
        mrk["z"][:]      = 0
        mrk["zeta"][:]   = 0
        mrk["energy"][:] = 100e6
        mrk["pitch"][:]  = 0.5

        for tag in [PhysTest.tag_elementary_gyro,
                    PhysTest.tag_elementary_exbgo,
                    PhysTest.tag_elementary_exbgc,
                    PhysTest.tag_elementary_gradbgo,
                    PhysTest.tag_elementary_gradbgc]:
            init("gc", **mrk, desc=tag)

    def run_elementary(self):
        """Run elementary test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_elementary_gyro):
            warnings.warn("Results already present: Test elementary")
            return
        for tag in [PhysTest.tag_elementary_gyro,
                    PhysTest.tag_elementary_exbgo,
                    PhysTest.tag_elementary_exbgc,
                    PhysTest.tag_elementary_gradbgo,
                    PhysTest.tag_elementary_gradbgc]:
            self._activateinputs(tag)
            self._runascot(tag)

    def check_elementary(self):
        """Verify and plot elementary test results.

        Returns
        -------
        passed : bool
            True if the test passed.
        """
        run_gyro    = self.ascot.data[PhysTest.tag_elementary_gyro]
        run_exbgo   = self.ascot.data[PhysTest.tag_elementary_exbgo]
        run_exbgc   = self.ascot.data[PhysTest.tag_elementary_exbgc]
        run_gradbgo = self.ascot.data[PhysTest.tag_elementary_gradbgo]
        run_gradbgc = self.ascot.data[PhysTest.tag_elementary_gradbgc]

        # Initialize plots
        fig = a5plt.figuredoublecolumn()
        gs = GridSpec(2, 3, figure=fig)
        h1a = fig.add_subplot(gs[0,0])
        h2a = fig.add_subplot(gs[0,1])
        h2b = fig.add_subplot(gs[1,1])
        h3a = fig.add_subplot(gs[0,2])
        h3b = fig.add_subplot(gs[1,2])

        h1a.set_title("Gyro-motion")
        h2a.set_title("ExB drift")
        h3a.set_title("Grad-B drift")

        h1a.set_ylabel("z [m]")
        h1a.set_xlabel("y [m]")

        h1a.set_aspect("equal", adjustable="box")
        h2a.set_aspect("equal", adjustable="box")
        h2b.set_aspect("equal", adjustable="box")
        h3a.set_aspect("equal", adjustable="box")
        h3b.set_aspect("equal", adjustable="box")

        def plotarrow(ax, ini, end, txt):
            """Plot line with an arrow
            """
            ax.annotate("", xytext=ini, xy=end,
                        arrowprops={"arrowstyle":"->", "color":"black"})
            ax.plot([ini[0], end[0]], [ini[1], end[1]], color="black")
            ax.annotate(txt, ini + .005)

        ## Gyro test ##
        m, q, pitch, ekin, x0, y0 = \
            run_gyro.getstate("mass", "charge", "pitch", "ekin", "y", "z")
        x, y, time = run_gyro.getorbit("y", "z", "time", ids=1)

        # Analytical values
        bnorm = np.sqrt(np.sum(run_gyro.bfield.read()["bxyz"]**2)) * unyt.T
        larmorrad_ana = physlib.gyrolength(m, q, ekin, pitch, bnorm).to("m")
        gyrofreq_ana  = physlib.gyrofrequency(m, q, ekin, bnorm).to("rad/s")

        # Numerical values
        larmorrad_go = np.mean( np.sqrt( (x - x0[0])**2 + (y - y0[0])**2 ) )
        gyrofreq_go  = np.sum( np.sqrt( np.diff(x)**2 + np.diff(y)**2 ) ) \
            * unyt.rad / (larmorrad_go * time[-1])

        # Plot
        orbx = x - x0[0]
        orby = y - y0[0]

        anax   = larmorrad_ana[0] * np.sin(np.linspace(0, 2*np.pi, 360))
        anay   = larmorrad_ana[0] * np.cos(np.linspace(0, 2*np.pi, 360))
        fangle = np.arctan2(orby[0], orbx[0]) * unyt.rad \
            + gyrofreq_ana[0] * time[-1]
        fx = larmorrad_ana[0] * np.cos(-fangle) # Minus sign because this we are
        fy = larmorrad_ana[0] * np.sin(-fangle) # in LHS coordinate system.

        h1a.plot(orbx, orby, linewidth=3)
        h1a.plot(anax, anay, linestyle="--", color="black", alpha=0.7)
        h1a.scatter(orbx[-1], orby[-1], marker="o", color="C0", zorder=3)
        h1a.scatter(fx, fy, marker="x", color="black", zorder=4)

        ## E x B drift
        yi_go            = run_exbgo.getstate("z", state="ini")
        yf_go, deltat_go = run_exbgo.getstate("z", "mileage", state="end")
        x_go1, y_go1     = run_exbgo.getorbit("y", "z", ids=1)
        x_go2, y_go2     = run_exbgo.getorbit("y", "z", ids=2)

        yi_gc            = run_exbgc.getstate("z", state="ini")
        yf_gc, deltat_gc = run_exbgc.getstate("z", "mileage", state="end")
        x_gc1, y_gc1     = run_exbgc.getorbit("y", "z", ids=1)
        x_gc2, y_gc2     = run_exbgc.getorbit("y", "z", ids=2)

        # Analytical values
        bvec  = run_exbgo.bfield.read()["bxyz"].ravel()
        evec  = run_exbgo.efield.read()["exyz"].ravel()
        v_ExB = np.cross(evec, bvec) / np.inner(bvec, bvec) * unyt.m / unyt.s
        time  = run_exbgo.options.read()["ENDCOND_LIM_SIMTIME"] * unyt.s

        # Numerical values
        vgo_ExB = (yf_go - yi_go) / deltat_go
        vgc_ExB = (yf_gc - yi_gc) / deltat_gc

        # Plot
        h2a.plot(x_go1, y_go1)
        h2b.plot(x_go2, y_go2)
        h2a.plot(x_gc1, y_gc1)
        h2b.plot(x_gc2, y_gc2)

        ini = np.array([x_gc1[0].v + .01, y_gc1[0].v])
        end = np.array([x_gc1[0].v + .01, (y_gc1[0] + v_ExB[2]*deltat_gc[1]).v])
        plotarrow(h2a, ini, end, r"$e^+$")
        ini = np.array([x_gc2[0].v + .01, y_gc2[0].v])
        end = np.array([x_gc2[0].v + .01, (y_gc2[0] + v_ExB[2]*deltat_gc[0]).v])
        plotarrow(h2b, ini, end, r"$e^-$")

        ## gradB drift
        pitch, ekin, q, m = run_gradbgo.getstate(
            "pitch", "ekin", "charge", "mass", ids=1)
        xi_go            = run_gradbgo.getstate("y", state="ini")
        xf_go, deltat_go = run_gradbgo.getstate("y", "mileage", state="end")
        x_go1, y_go1     = run_gradbgo.getorbit("y", "z", ids=1)
        x_go2, y_go2     = run_gradbgo.getorbit("y", "z", ids=2)

        xi_gc            = run_gradbgc.getstate("y", state="ini")
        xf_gc, deltat_gc = run_gradbgc.getstate("y", "mileage", state="end")
        x_gc1, y_gc1     = run_gradbgc.getorbit("y", "z", ids=1)
        x_gc2, y_gc2     = run_gradbgc.getorbit("y", "z", ids=2)

        # Analytical values
        gamma   = physlib.gamma_energy(m, ekin)
        vnorm   = physlib.vnorm_gamma(gamma)
        bvec    = run_gradbgo.bfield.read()["bxyz"].ravel()
        gradb   = run_gradbgo.bfield.read()["jacobian"][0,:]
        bnorm   = np.sqrt(np.sum(bvec**2)) * unyt.T
        v_gradb = gamma * m * ( (1.0 - pitch**2) * vnorm**2 / ( q * bnorm ) ) \
            * 0.5 * np.cross(bvec, gradb) / bnorm**2 * unyt.T**2 / unyt.m
        v_gradb.convert_to_units("m/s")
        time    = run_gradbgo.options.read()["ENDCOND_LIM_SIMTIME"] * unyt.s

        # Numerical values
        vgo_gradb = (xf_go - xi_go) / deltat_go
        vgc_gradb = (xf_gc - xi_gc) / deltat_gc

        # Plot
        h3a.plot(x_go1, y_go1)
        h3b.plot(x_go2, y_go2)
        h3a.plot(x_gc1, y_gc1)
        h3b.plot(x_gc2, y_gc2)

        ini = np.array([x_gc1[0].v, y_gc1[0].v + .01])
        end = np.array([x_gc1[0] + v_gradb[1]*deltat_gc[0], y_gc1[0].v + .01])
        plotarrow(h3a, ini, end, r"$e^+$")
        ini = np.array([x_gc2[0].v, y_gc2[0] + .01])
        end = np.array([x_gc2[0] - v_gradb[1]*deltat_gc[0], y_gc2[0].v + .01])
        plotarrow(h3b, ini, end, r"$e^-$")

        print("Test elementary:")
        passed = True

        print("  Gyroradius %e, Gyrofrequency %e (gyro-orbit)" \
              % (larmorrad_go, gyrofreq_go))
        err = np.abs(larmorrad_go - larmorrad_ana[0]) / larmorrad_ana[0]
        if err > 5e-8: print("Error: %e (FAILED)" % err); passed = False

        print("  Gyroradius %e, Gyrofrequency %e (expected)" \
              % (larmorrad_ana[0], gyrofreq_ana[0]))
        err = np.abs(gyrofreq_go - gyrofreq_ana[0]) / gyrofreq_ana[0]
        if err > 5e-4: print("Error: %e (FAILED)" % err); passed = False

        print("  ExB drift for positrons %e electrons %e (gyro-orbit)" \
              % (vgo_ExB[0], vgo_ExB[1]) )
        print("  ExB drift for positrons %e electrons %e (guiding center)" \
              % (vgc_ExB[0], vgc_ExB[1]) )
        print("  ExB drift for positrons %e electrons %e (expected)" \
              % (v_ExB[2], v_ExB[2]) )
        err = np.amax([np.abs( ( vgo_ExB[0] - v_ExB[2] ) / v_ExB[2] ),
                       np.abs( ( vgo_ExB[1] - v_ExB[2] ) / v_ExB[2] ),
                       np.abs( ( vgc_ExB[0] - v_ExB[2] ) / v_ExB[2] ),
                       np.abs( ( vgc_ExB[1] - v_ExB[2] ) / v_ExB[2] )])
        if err > 5e-5: print("Error: %e (FAILED)" % err); passed = False

        print("  Grad-B drift for positrons %e electrons %e (gyro-orbit)" \
              % (vgo_gradb[0], vgo_gradb[1]) )
        print("  Grad-B drift for positrons %e electrons %e (guiding center)" \
              % (vgc_gradb[0], vgc_gradb[1]) )
        print("  Grad-B drift for positrons %e electrons %e (expected)" \
              % (v_gradb[1], -v_gradb[1]) )
        err = np.amax([np.abs( ( vgo_gradb[0] - v_gradb[1] ) / v_gradb[1] ),
                       np.abs( ( vgo_gradb[1] + v_gradb[1] ) / v_gradb[1] ),
                       np.abs( ( vgc_gradb[0] - v_gradb[1] ) / v_gradb[1] ),
                       np.abs( ( vgc_gradb[1] + v_gradb[1] ) / v_gradb[1] )])
        if err > 1e-3: print("Error: %e (FAILED)" % err); passed = False

        return passed

    def init_orbitfollowing(self):
        """Initialize data for the  test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_orbfol_go):
            warnings.warn("Inputs already present: Test orbit-following")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_SIMTIMELIM" : 1,
            "ENDCOND_LIM_SIMTIME" : 5e-6, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
            "ORBITWRITE_INTERVAL" : 1e-10, "ORBITWRITE_NPOINT" : 50002
        })
        init("opt", **opt, desc=PhysTest.tag_orbfol_go)
        opt.update({
            "SIM_MODE" : 2, "FIXEDSTEP_USERDEFINED" : 1e-12,
            "ENABLE_ADAPTIVE" : 0,
            "ORBITWRITE_INTERVAL" : 1e-8, "ORBITWRITE_NPOINT" : 502
        })
        init("opt", **opt, desc=PhysTest.tag_orbfol_gcf)
        opt.update({
            "ENABLE_ADAPTIVE" : 1, "ADAPTIVE_MAX_DRHO" : 0.1,
            "ADAPTIVE_MAX_DPHI" : 10, "ADAPTIVE_TOL_ORBIT" : 1e-11,
            "FIXEDSTEP_USERDEFINED" : 1e-8
        })
        init("opt", **opt, desc=PhysTest.tag_orbfol_gca)

        # Magnetic field is just some tokamak
        init("bfield_analytical_iter_circular", desc=PhysTest.tag_orbfol_go)
        init("bfield_analytical_iter_circular", desc=PhysTest.tag_orbfol_gcf)
        init("bfield_analytical_iter_circular", desc=PhysTest.tag_orbfol_gca)

        # Marker input is a trapped positron and a passing electron
        mrk = Marker.generate("gc", n=2, species="electron")
        mrk["charge"]    = np.array([1, -1])
        mrk["r"][:]      = 7.6
        mrk["phi"][:]    = 90
        mrk["z"][:]      = 0
        mrk["zeta"][:]   = 2
        mrk["energy"][:] = 10e6
        mrk["pitch"]     = np.array([0.4, 0.9])
        for tag in [PhysTest.tag_orbfol_go, PhysTest.tag_orbfol_gcf,
                    PhysTest.tag_orbfol_gca]:
            init("gc", **mrk, desc=tag)

    def run_orbitfollowing(self):
        """Run orbit-following test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_orbfol_go):
            warnings.warn("Results already present: Test orbit-following")
            return
        for tag in [PhysTest.tag_orbfol_go, PhysTest.tag_orbfol_gcf,
                    PhysTest.tag_orbfol_gca]:
            self._activateinputs(tag)
            self._runascot(tag)

    def check_orbitfollowing(self):
        """Check test.
        """
        run_go  = self.ascot.data[PhysTest.tag_orbfol_go]
        run_gcf = self.ascot.data[PhysTest.tag_orbfol_gcf]
        run_gca = self.ascot.data[PhysTest.tag_orbfol_gca]

        # Initialize plots
        fig = a5plt.figuredoublecolumn(3/2)
        gs = GridSpec(3, 4, figure=fig)
        h1a = fig.add_subplot(gs[0,0])
        h2a = fig.add_subplot(gs[1,0])
        h3a = fig.add_subplot(gs[2,0])
        h4a = fig.add_subplot(gs[:,1])
        h1b = fig.add_subplot(gs[0,2])
        h2b = fig.add_subplot(gs[1,2])
        h3b = fig.add_subplot(gs[2,2])
        h4b = fig.add_subplot(gs[:,3])

        h4a.set_aspect("equal", adjustable="box")
        h4b.set_aspect("equal", adjustable="box")

        h1a.set_xlim(0, 5)
        h2a.set_xlim(0, 5)
        h3a.set_xlim(0, 5)
        h1b.set_xlim(0, 5)
        h2b.set_xlim(0, 5)
        h3b.set_xlim(0, 5)

        h4a.set_xlim(5, 8)
        h4a.set_ylim(-2, 2)
        h4b.set_xlim(5, 8)
        h4b.set_ylim(-2, 2)

        h1a.set_xticklabels([])
        h2a.set_xticklabels([])
        h1b.set_xticklabels([])
        h2b.set_xticklabels([])

        h4a.set_yticklabels([])
        h4b.set_xlabel("R [m]")
        h4b.set_ylabel("z [m]")
        h4b.yaxis.set_label_position("right")
        h4b.yaxis.tick_right()

        h1a.set_ylabel(r"$(E-E_0)/E_0$")
        h2a.set_ylabel(r"$(\mu-\mu_0)/\mu_0$")
        h3a.set_ylabel(r"$(P-P_0)/P_0$")

        h3a.set_xlabel("Time [µs]")
        h3b.set_xlabel("Time [µs]")

        def plotreldiff(ax1, ax2, ax3, t, q1, q2, q3, **kwargs):
            """Plot relative diffenrence of given quantities on given axes.
            """
            ax1.plot(t, ( q1 - q1[0] ) / q1[0], **kwargs)
            ax2.plot(t, ( q2 - q2[0] ) / q2[0], **kwargs)
            ax3.plot(t, ( q3 - q3[0] ) / q3[0], **kwargs)

        def fails(t, q, eps, qnt, otype, mode):
            """Check if the change in time of a given quantity is below given
            tolerance
            """
            err = np.polyfit(t, (q - q[0]) / q[0], 1)[0]
            msg = "Rate of change in %6s (%s/%11s): %e Tolerance: %e" \
                % (qnt, otype, mode, err, eps)
            if np.abs(err) > eps:
                msg += " (FAILED)"
                print(msg)
                return True
            print(msg)
            return False

        # Numerical values
        self.ascot.input_init(run=run_go.get_qid(), bfield=True)
        tgo1, ego1, mugo1, pgo1, rgo1, zgo1 = run_go.getorbit(
            "mileage", "ekin", "mu", "ptor", "r", "z", ids=1)
        tgo2, ego2, mugo2, pgo2, rgo2, zgo2 = run_go.getorbit(
            "mileage", "ekin", "mu", "ptor", "r", "z", ids=2)
        tgcf1, egcf1, mugcf1, pgcf1, rgcf1, zgcf1 = run_gcf.getorbit(
            "mileage", "ekin", "mu", "ptor", "r", "z", ids=1)
        tgcf2, egcf2, mugcf2, pgcf2, rgcf2, zgcf2 = run_gcf.getorbit(
            "mileage", "ekin", "mu", "ptor", "r", "z", ids=2)
        tgca1, egca1, mugca1, pgca1, rgca1, zgca1 = run_gca.getorbit(
            "mileage", "ekin", "mu", "ptor", "r", "z", ids=1)
        tgca2, egca2, mugca2, pgca2, rgca2, zgca2 = run_gca.getorbit(
            "mileage", "ekin", "mu", "ptor", "r", "z", ids=2)
        self.ascot.input_free()

        # Plot
        plotreldiff(h1a, h2a, h3a, tgo1.to("µs"),  ego1,  mugo1,  pgo1,
                    color="C0")
        plotreldiff(h1b, h2b, h3b, tgo2.to("µs"),  ego2,  mugo2,  pgo2,
                    color="C0")
        plotreldiff(h1a, h2a, h3a, tgcf1.to("µs"), egcf1, mugcf1, pgcf1,
                    color="C1")
        plotreldiff(h1b, h2b, h3b, tgcf2.to("µs"), egcf2, mugcf2, pgcf2,
                    color="C1")
        plotreldiff(h1a, h2a, h3a, tgca1.to("µs"), egca1, mugca1, pgca1,
                    color="C2")
        plotreldiff(h1b, h2b, h3b, tgca2.to("µs"), egca2, mugca2, pgca2,
                    color="C2")

        h4a.plot(rgo1,  zgo1,  color="C0")
        h4a.plot(rgcf1, zgcf1, color="C1")
        h4a.plot(rgca1, zgca1, color="C2")
        h4b.plot(rgo2,  zgo2,  color="C0")
        h4b.plot(rgcf2, zgcf2, color="C1")
        h4b.plot(rgca2, zgca2, color="C2")

        print("Test orbit-following:")
        passed = True

        otype1 = "passing"; otype2 = "trapped"
        simmode1 = "gyro-orbit"; simmode2 = "fixed GC"; simmode3 = "adaptive GC"
        if fails(tgo1,  ego1,   5e-6, "energy", otype1, simmode1):passed = False
        if fails(tgo2,  ego2,   5e-6, "energy", otype2, simmode1):passed = False
        if fails(tgcf1, egcf1,  5e-6, "energy", otype1, simmode2):passed = False
        if fails(tgcf2, egcf2,  5e-6, "energy", otype2, simmode2):passed = False
        if fails(tgca1, egca1,  5e-6, "energy", otype1, simmode3):passed = False
        if fails(tgca2, egca2,  5e-6, "energy", otype2, simmode3):passed = False
        if fails(tgo1,  mugo1,   2e2, "mu",     otype1, simmode1):passed = False
        if fails(tgo2,  mugo2,  2e-1, "mu",     otype2, simmode1):passed = False
        if fails(tgcf1, mugcf1, 1e-9, "mu",     otype1, simmode2):passed = False
        if fails(tgcf2, mugcf2, 1e-9, "mu",     otype2, simmode2):passed = False
        if fails(tgca1, mugca1, 1e-9, "mu",     otype1, simmode3):passed = False
        if fails(tgca2, mugca2, 1e-9, "mu",     otype2, simmode3):passed = False
        if fails(tgo1,  pgo1,   5e-3, "Ptor",   otype1, simmode1):passed = False
        if fails(tgo2,  pgo2,   5e-5, "Ptor",   otype2, simmode1):passed = False
        if fails(tgcf1, pgcf1,  5e-2, "Ptor",   otype1, simmode2):passed = False
        if fails(tgcf2, pgcf2,  5e-2, "Ptor",   otype2, simmode2):passed = False
        if fails(tgca1, pgca1,  5e-2, "Ptor",   otype1, simmode3):passed = False
        if fails(tgca2, pgca2,  5e-2, "Ptor",   otype2, simmode3):passed = False

        return passed

    def init_gctransform(self):
        """Initialize data for the guiding-center transformation test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_gctransform_go):
            warnings.warn("Inputs already present: Test GC transform")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "SIM_MODE" : 2, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "FIXEDSTEP_USERDEFINED" : 1e-10, "ENDCOND_SIMTIMELIM" : 1,
            "ENDCOND_LIM_SIMTIME" : 3e-5, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
            "ORBITWRITE_INTERVAL" : 4e-10, "ORBITWRITE_NPOINT" : 75002
        })
        init("opt", **opt, desc=PhysTest.tag_gctransform_gc)
        init("opt", **opt, desc=PhysTest.tag_gctransform_first)
        opt.update({"DISABLE_FIRSTORDER_GCTRANS" : 1})
        init("opt", **opt, desc=PhysTest.tag_gctransform_zeroth)

        opt.update({"DISABLE_FIRSTORDER_GCTRANS" : 0, "SIM_MODE" : 1 })
        init("opt", **opt, desc=PhysTest.tag_gctransform_go)
        opt.update({"RECORD_MODE" : 1 })
        init("opt", **opt, desc=PhysTest.tag_gctransform_go2gc)

        # Magnetic field is just some tokamak
        init("bfield_analytical_iter_circular",
             desc=PhysTest.tag_gctransform_go)
        init("bfield_analytical_iter_circular",
             desc=PhysTest.tag_gctransform_gc)
        init("bfield_analytical_iter_circular",
             desc=PhysTest.tag_gctransform_go2gc)
        init("bfield_analytical_iter_circular",
             desc=PhysTest.tag_gctransform_zeroth)
        init("bfield_analytical_iter_circular",
             desc=PhysTest.tag_gctransform_first)

        # Use single alpha particle in tests
        mrk = Marker.generate("gc", n=1, species="alpha")
        mrk["r"][:]      = 7.6
        mrk["phi"][:]    = 90
        mrk["z"][:]      = 0
        mrk["zeta"][:]   = 2
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.4
        init("gc", **mrk, desc=PhysTest.tag_gctransform_go)
        init("gc", **mrk, desc=PhysTest.tag_gctransform_gc)
        init("gc", **mrk, desc=PhysTest.tag_gctransform_go2gc)

    def run_gctransform(self):
        """Run GC transform test.
        """
        init = self.ascot.data.create_input
        if hasattr(self.ascot.data, PhysTest.tag_gctransform_go):
            warnings.warn("Results already present: Test GC transform")
            return
        for tag in [PhysTest.tag_gctransform_go, PhysTest.tag_gctransform_gc,
                    PhysTest.tag_gctransform_go2gc]:
            self._activateinputs(tag)
            self._runascot(tag)

        # Create new marker input from results
        nrep = 10
        mrk = Marker.generate("prt", n=nrep, species="alpha")
        run = self.ascot.data[PhysTest.tag_gctransform_go]

        # Pick initial coordinates along the orbit trajectory
        dt   = 20
        idx = np.s_[:nrep*dt:dt]
        time, r, phi, z, vr, vphi, vz = \
            run.getorbit("time", "r", "phi", "z", "vr", "vphi", "vz")
        mrk.update({
            "time" : time[idx], "r" : r[idx], "phi" : phi[idx], "z" : z[idx],
            "vr" : vr[idx], "vphi" : vphi[idx], "vz" : vz[idx]
        })
        init("prt", **mrk, desc=PhysTest.tag_gctransform_zeroth)
        init("prt", **mrk, desc=PhysTest.tag_gctransform_first)

        for tag in [PhysTest.tag_gctransform_zeroth,
                    PhysTest.tag_gctransform_first]:
            self._activateinputs(tag)
            self._runascot(tag)

    def check_gctransform(self):
        """Check test.
        """
        run_go     = self.ascot.data[PhysTest.tag_gctransform_go]
        run_gc     = self.ascot.data[PhysTest.tag_gctransform_gc]
        run_go2gc  = self.ascot.data[PhysTest.tag_gctransform_go2gc]
        run_zeroth = self.ascot.data[PhysTest.tag_gctransform_zeroth]
        run_first  = self.ascot.data[PhysTest.tag_gctransform_first]

        # Initialize plots
        fig = a5plt.figuredoublecolumn(3/2)
        gs  = GridSpec(3, 3, figure=fig)
        h1a = fig.add_subplot(gs[0,0])
        h1b = fig.add_subplot(gs[1,0])
        h1c = fig.add_subplot(gs[2,0])
        h2  = fig.add_subplot(gs[:,1])
        h3  = fig.add_subplot(gs[:,2])

        h1a.set_xlim(0, 30)
        h1b.set_xlim(0, 30)
        h1c.set_xlim(0, 30)
        h1a.set_xticklabels([])
        h1b.set_xticklabels([])

        h2.set_xlim(6.2, 6.7)
        h2.set_ylim(1.0, 1.6)
        h2.set_aspect("equal", adjustable="box")
        h3.set_xlim(6.2, 6.7)
        h3.set_ylim(1.0, 1.6)
        h3.set_aspect("equal", adjustable="box")
        h3.set_yticklabels([])

        h1c.set_xlabel("Time [µs]")
        h1a.set_ylabel(r"$\mu$")
        h1b.set_ylabel(r"$E_\mathrm{kin}$")
        h1c.set_ylabel(r"$P_\mathrm{ctor}$")
        h2.set_xlabel("R [m]")
        h3.set_xlabel("R [m]")
        h3.set_ylabel("z [m]")
        h3.yaxis.set_label_position("right")
        h3.yaxis.tick_right()

        # Get data and plot
        self.ascot.input_init(run=run_go.get_qid(), bfield=True)
        rgo, zgo, tgo, ego, mugo, pgo = run_go.getorbit(
            "r", "z", "mileage", "ekin", "mu", "ptor")
        rgc, zgc, tgc, egc, mugc, pgc = run_gc.getorbit(
            "r", "z", "mileage", "ekin", "mu", "ptor")
        rgo2gc, zgo2gc, tgo2gc, ego2gc, mugo2gc, pgo2gc = run_go2gc.getorbit(
            "r", "z", "mileage", "ekin", "mu", "ptor")
        self.ascot.input_free()

        h1a.plot(tgo.to("µs"), mugo)
        h1a.plot(tgo2gc.to("µs"), mugo2gc)
        h1a.plot(tgc.to("µs"), mugc)

        h1b.plot(tgo.to("µs"), ego)
        h1b.plot(tgo2gc.to("µs"), ego2gc)
        h1b.plot(tgc.to("µs"), egc)

        h1c.plot(tgo.to("µs"), pgo)
        h1c.plot(tgo2gc.to("µs"), pgo2gc)
        h1c.plot(tgc.to("µs"), pgc)

        h2.plot(rgo, zgo)
        h2.plot(rgo2gc, zgo2gc)
        h2.plot(rgc, zgc)

        nrep = run_zeroth.getstate("ids").size
        for i in range(nrep):
            r, z = run_zeroth.getorbit("r", "z", ids=i)
            h3.plot(r, z, color="C1")
        for i in range(nrep):
            r, z = run_first.getorbit("r", "z", ids=i)
            h3.plot(r.v+0.01, z, color="C2")
        h3.plot(rgo2gc, zgo2gc, color="black")

        # Verify results by interpolating the orbits at fixed intervals and then
        # calculating sum-of-squares of the difference between go and gc and
        # go2gc and gc. The latter should be smaller if the guiding center
        # transformation works.
        t = np.linspace(0, 1e-6, 1000)
        mugo    = np.interp(t, tgo,    mugo)
        mugc    = np.interp(t, tgc,    mugc)
        mugo2gc = np.interp(t, tgo2gc[1:], mugo2gc[1:])
        ego     = np.interp(t, tgo,    ego)
        egc     = np.interp(t, tgc,    egc)
        ego2gc  = np.interp(t, tgo2gc[1:], ego2gc[1:])
        pgo     = np.interp(t, tgo,    pgo)
        pgc     = np.interp(t, tgc,    pgc)
        pgo2gc  = np.interp(t, tgo2gc[1:], pgo2gc[1:])

        print("Test GC transformation:")
        passed = True

        mugo = np.mean(mugo); mugc = np.mean(mugc); mugo2gc = np.mean(mugo2gc)
        ego  = np.mean(ego);  egc  = np.mean(egc);  ego2gc  = np.mean(ego2gc)
        pgo  = np.mean(pgo);  pgc  = np.mean(pgc);  pgo2gc  = np.mean(pgo2gc)

        print("  Mean value  GO        GC        GO2GC")
        err = np.abs(mugc/mugo2gc-1)
        print("  mu          %1.3e %1.3e %1.3e" % (mugo, mugc, mugo2gc))
        if err > 1e-4: print("Error %e (FAILED)" % err); passed = False
        err = np.abs(egc/ego2gc-1)
        print("  Energy      %1.3e %1.3e %1.3e" % (ego, egc, ego2gc))
        if err > 2e-4: print("Error %e (FAILED)" % err); passed = False
        err = np.abs(pgc/pgo2gc-1)
        print("  Pctor       %1.3e %1.3e %1.3e" % (-pgo, -pgc, -pgo2gc))
        if err > 1e-4: print("Error %e (FAILED)" % err); passed = False

        return passed

    def init_ccoll(self):
        """Initialize data for the Coulomb collision test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_ccoll_thermalgo):
            warnings.warn("Inputs already present: Test Coulomb collision")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "ENDCOND_SIMTIMELIM" : 1, "ENDCOND_LIM_SIMTIME" : 2e-2,
            "ENABLE_ORBIT_FOLLOWING" : 1, "ENABLE_COULOMB_COLLISIONS" : 1,
            "ENABLE_DIST_5D" : 1,
            "DIST_MIN_R"    : 4,  "DIST_MAX_R"    : 10, "DIST_NBIN_R"      : 1,
            "DIST_MIN_PHI"  : 0,  "DIST_MAX_PHI"  : 360, "DIST_NBIN_PHI"   : 1,
            "DIST_MIN_Z"    : -5, "DIST_MAX_Z"    : 5, "DIST_NBIN_Z"       : 1,
            "DIST_MIN_TIME" : 0,  "DIST_MAX_TIME" : 2e-2, "DIST_NBIN_TIME" : 2,
            "DIST_MIN_PPA"  : -2.5e-21, "DIST_MAX_PPA" : 2.5e-21,
            "DIST_NBIN_PPA" : 140, "DIST_MIN_PPE" : 0, "DIST_MAX_PPE" : 2.5e-21,
            "DIST_NBIN_PPE" : 80
        })

        opt.update({"SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
                    "FIXEDSTEP_USERDEFINED" : 1e-8})
        init("opt", **opt, desc=PhysTest.tag_ccoll_thermalgo)
        opt.update({"SIM_MODE" : 2, "FIXEDSTEP_USE_USERDEFINED" : 1,
                    "FIXEDSTEP_USERDEFINED" : 2e-8})
        init("opt", **opt, desc=PhysTest.tag_ccoll_thermalgcf)
        opt.update({"ENABLE_ADAPTIVE" : 1, "ADAPTIVE_TOL_ORBIT" : 1e-6,
                    "ADAPTIVE_TOL_CCOL" : 1e-2, "ADAPTIVE_MAX_DRHO" : 0.1,
                    "ADAPTIVE_MAX_DPHI" : 10, "FIXEDSTEP_USERDEFINED" : 1e-8})
        init("opt", **opt, desc=PhysTest.tag_ccoll_thermalgca)

        opt = Opt.get_default()
        opt.update({
            "ENDCOND_ENERGYLIM" : 1, "ENDCOND_MIN_ENERGY" : 50e3,
            "ENDCOND_MIN_THERMAL" : 0, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_COULOMB_COLLISIONS" : 1, "ENABLE_DIST_5D" : 1,
            "DIST_MIN_R" : 4, "DIST_MAX_R" : 10, "DIST_NBIN_R" : 1,
            "DIST_MIN_PHI" : 0, "DIST_MAX_PHI" : 360, "DIST_NBIN_PHI" : 1,
            "DIST_MIN_Z" : -5, "DIST_MAX_Z" : 5, "DIST_NBIN_Z" : 1,
            "DIST_MIN_TIME" : 0, "DIST_MAX_TIME" : 1e0, "DIST_NBIN_TIME" : 1,
            "DIST_MIN_PPA" : -1.3e-19, "DIST_MAX_PPA" : 1.3e-19,
            "DIST_NBIN_PPA" : 200, "DIST_MIN_PPE" : 0, "DIST_MAX_PPE" : 1.3e-19,
            "DIST_NBIN_PPE" : 100
        })

        opt.update({"SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
                    "FIXEDSTEP_USERDEFINED" : 2e-9})
        init("opt", **opt, desc=PhysTest.tag_ccoll_slowinggo)
        opt.update({"SIM_MODE" : 2, "FIXEDSTEP_USE_USERDEFINED" : 1,
                    "FIXEDSTEP_USERDEFINED" : 3e-8})
        init("opt", **opt, desc=PhysTest.tag_ccoll_slowinggcf)
        opt.update({"ENABLE_ADAPTIVE" : 1, "ADAPTIVE_TOL_ORBIT" : 1e-6,
                    "ADAPTIVE_TOL_CCOL" : 1e-2, "ADAPTIVE_MAX_DRHO" : 0.1,
                    "ADAPTIVE_MAX_DPHI" : 10, "FIXEDSTEP_USERDEFINED" : 1e-8})
        init("opt", **opt, desc=PhysTest.tag_ccoll_slowinggca)

        # Magnetic field is just some tokamak and plasma is uniform
        for tag in [PhysTest.tag_ccoll_thermalgo, PhysTest.tag_ccoll_thermalgcf,
                    PhysTest.tag_ccoll_thermalgca, PhysTest.tag_ccoll_slowinggo,
                    PhysTest.tag_ccoll_slowinggcf, PhysTest.tag_ccoll_slowinggca
                    ]:
            init("bfield_analytical_iter_circular", desc=tag)
            init("plasma_flat", density=1e20, temperature=1e3, desc=tag)

        mrk = Marker.generate("gc", n=20, species="proton")
        pol = 2 * np.pi * np.random.rand(mrk["n"],)
        mrk["r"][:]      = 6.2 + 0.8 * np.cos(pol)
        mrk["phi"][:]    = 90
        mrk["z"][:]      = 0.8 * np.sin(pol)
        mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)
        mrk["energy"][:] = 1.e3
        mrk["pitch"][:]  = 0.5
        init("gc", **mrk, desc=PhysTest.tag_ccoll_thermalgo)
        init("gc", **mrk, desc=PhysTest.tag_ccoll_thermalgcf)
        init("gc", **mrk, desc=PhysTest.tag_ccoll_thermalgca)

        mrk = Marker.generate("gc", n=200, species="alpha")
        pol = 2 * np.pi * np.random.rand(mrk["n"],)
        mrk["r"][:]      = 6.2 + 0.8 * np.cos(pol)
        mrk["phi"][:]    = 90
        mrk["z"][:]      = 0.8 * np.sin(pol)
        mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
        init("gc", **mrk, desc=PhysTest.tag_ccoll_slowinggo)
        init("gc", **mrk, desc=PhysTest.tag_ccoll_slowinggcf)
        init("gc", **mrk, desc=PhysTest.tag_ccoll_slowinggca)

    def run_ccoll(self):
        """Run Coulmb collision test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_ccoll_thermalgo):
            warnings.warn("Results already present: Test Coulomb collision")
            return
        for tag in [PhysTest.tag_ccoll_thermalgo, PhysTest.tag_ccoll_thermalgcf,
                    PhysTest.tag_ccoll_thermalgca, PhysTest.tag_ccoll_slowinggo,
                    PhysTest.tag_ccoll_slowinggcf, PhysTest.tag_ccoll_slowinggca
                    ]:
            self._activateinputs(tag)
            self._runascot(tag)

    def check_ccoll(self):
        """Check Coulomb collision test.
        """
        run_tgo  = self.ascot.data[PhysTest.tag_ccoll_thermalgo]
        run_tgcf = self.ascot.data[PhysTest.tag_ccoll_thermalgcf]
        run_tgca = self.ascot.data[PhysTest.tag_ccoll_thermalgca]
        run_sgo  = self.ascot.data[PhysTest.tag_ccoll_slowinggo]
        run_sgcf = self.ascot.data[PhysTest.tag_ccoll_slowinggcf]
        run_sgca = self.ascot.data[PhysTest.tag_ccoll_slowinggca]

        fig = a5plt.figuredoublecolumn()
        gs = GridSpec(2, 2, figure=fig)
        h1 = fig.add_subplot(gs[0,0])
        h2 = fig.add_subplot(gs[0,1])
        h3 = fig.add_subplot(gs[1,0])
        h4 = fig.add_subplot(gs[1,1])

        # Analytical results
        ne      = 1e20 / unyt.m**3
        Te      = 1e3 * unyt.eV
        m_e     = unyt.me
        m_p     = unyt.mp
        m_a     = 4.003*unyt.amu
        eps_0   = unyt.eps_0
        alphaZ  = 2
        clog    = 16
        E0      = 3.5e6*unyt.eV
        Emin    = 50 * Te

        # Thermal distribution with correct normalization
        thdist  = run_tgo.getdist("5d", exi=True)
        egridth = thdist.abscissa_edges("ekin")
        simtime = np.diff(thdist.abscissa_edges("time")[-2:])
        thermal = 2*np.sqrt(egridth/np.pi) * np.power(Te, -3.0/2) \
            * np.exp(-egridth/Te) * simtime

        vth   = np.sqrt(2*Te / m_e)
        vcrit = vth * np.power( (3.0*np.sqrt(np.pi)/4.0) * (m_e / m_p) , 1/3.0)
        Ecrit = (0.5 * m_a * vcrit * vcrit).to("eV")
        ts    = ( 3 * np.sqrt( (2*np.pi * Te)**3 / m_e ) * eps_0**2
                  * m_a / ( alphaZ**2 * unyt.e**4 * ne * clog ) ).to("s")

        egridsd   = np.linspace(Te, 1.1*E0, 100)
        heaviside = np.logical_and(egridsd <= E0, egridsd >= Emin)
        slowing   = heaviside * ts / ( ( 1 + np.power(Ecrit/egridsd, 3.0/2) )
                                       * 2 * egridsd )

        # ts is slowing down rate which gives the slowing down time as
        # t_sd = ts*log(v_0 / v_th) = 0.5*ts*log(E_0/E_th)
        slowingdowntime = 0.5 * ts * np.log( E0 / Emin )

        th_ekin  = np.zeros((3,))
        th_pitch = np.zeros((3,))
        for i, run in enumerate([run_tgo, run_tgcf, run_tgca]):
            dist  = run.getdist("5d", exi=True)
            edist = dist.integrate(
                r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[1:],
                charge=np.s_[:], pitch=np.s_[:], copy=True)
            xdist = dist.integrate(
                r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[1:],
                charge=np.s_[:], ekin=np.s_[:], copy=True)
            edist.plot(axes=h1)
            xdist.plot(axes=h2)
            th_pitch[i] = np.mean(run.getstate("pitch", state="end"))
            th_ekin[i]  = np.mean(run.getstate("ekin",  state="end"))

        sd_time  = np.zeros((3,))
        sd_pitch = np.zeros((3,))
        for i, run in enumerate([run_sgo, run_sgcf, run_sgca]):
            dist = run.getdist("5d", exi=True, ekin_edges=egridsd)
            edist = dist.integrate(
                r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[:],
                charge=np.s_[:], pitch=np.s_[:], copy=True)
            xdist = dist.integrate(
                r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[:],
                charge=np.s_[:], ekin=np.s_[:], copy=True)
            edist.plot(axes=h3)
            xdist.plot(axes=h4)
            sd_pitch[i] = np.mean(run.getstate("pitch", state="end"))
            sd_time[i]  = np.mean(run.getstate("time",  state="end"))

        # These are multiplied with marker number to get correct normalization
        Nmrkth = run_tgo.getstate("ids").size
        Nmrksd = run_sgo.getstate("ids").size
        h1.plot(egridth, thermal*Nmrkth, color="black")
        h3.plot(egridsd, slowing*Nmrksd, color="black")

        print("Test Coulomb collisions:")
        passed = True

        print("  Thermal final energy and pitch")
        print("  GO            %1.1e      %1.2f" % (th_ekin[0], th_pitch[0]))
        print("  GCF           %1.1e      %1.2f" % (th_ekin[1], th_pitch[1]))
        print("  GCA           %1.1e      %1.2f" % (th_ekin[2], th_pitch[2]))
        print("  Expected      %1.1e      %1.2f" % (Te, 0.0))
        if np.amax(np.abs(Te.v - th_ekin)) > 2e3 or \
           np.amax(np.abs(0.0 - th_pitch)) > 0.5:
            print("  (Failed)")
            passed = False
        print("")
        print("  Slowing-down final time  and  pitch")
        print("  GO            %1.1e      %1.2f" % (sd_time[0], sd_pitch[0]))
        print("  GCF           %1.1e      %1.2f" % (sd_time[1], sd_pitch[1]))
        print("  GCA           %1.1e      %1.2f" % (sd_time[2], sd_pitch[2]))
        print("  Expected      %1.1e      %1.2f" % (slowingdowntime, 0.0))
        if np.amax(np.abs(slowingdowntime.v  - sd_time))  > 2e-3 or \
           np.amax(np.abs(0.0 - th_pitch)) > 0.2:
            print("  (Failed)")
            passed = False

        dist = run_tgo.getdist("5d")
        dist.integrate(
            r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[:],
            charge=np.s_[:])
        dist.plot()

        return passed

    def init_classical(self):
        """Initialize data for the classical transport test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_classical_go + "0"):
            warnings.warn("Inputs already present: Test classical transport")
            return
        init = self.ascot.data.create_input

        # Options
        optgo = Opt.get_default()
        optgo.update({
            "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "FIXEDSTEP_USERDEFINED" : 1e-10, "ENDCOND_SIMTIMELIM" : 1,
            "ENDCOND_LIM_SIMTIME" : 1e-5, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_COULOMB_COLLISIONS" : 1
        })
        optgcf = copy.deepcopy(optgo)
        optgcf.update({
            "SIM_MODE" : 2, "FIXEDSTEP_USERDEFINED" : 1e-9
        })
        optgca = copy.deepcopy(optgcf)
        optgca.update({
            "ENABLE_ADAPTIVE" : 1, "ADAPTIVE_MAX_DRHO" : 0.1,
            "ADAPTIVE_TOL_ORBIT" : 1e-8, "ADAPTIVE_TOL_CCOL" : 1e-1,
            "ADAPTIVE_MAX_DPHI" : 10, "FIXEDSTEP_USERDEFINED" : 1e-8
        })

        # Marker input consists of protons
        mrk = Marker.generate("gc", n=1000, species="proton")
        mrk["r"][:]      = 5.0
        mrk["phi"][:]    = 0.0
        mrk["z"][:]      = 0.0
        mrk["energy"][:] = 1e3
        mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
        mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)

        # Plasma consisting of electrons only to avoid proton-proton collisions
        pls = init("plasma_flat", density=1e22, temperature=1e3, dryrun=True)
        pls["idensity"][:] = 1

        for i in range(6):
            init("opt", **optgo,  desc=PhysTest.tag_classical_go  + str(i))
            init("opt", **optgcf, desc=PhysTest.tag_classical_gcf + str(i))
            init("opt", **optgca, desc=PhysTest.tag_classical_gca + str(i))
            for tag in [PhysTest.tag_classical_go, PhysTest.tag_classical_gcf,
                        PhysTest.tag_classical_gca]:
                # Transport is scanned as a function of magnetic field strength
                d = {"bxyz" : np.array([1.0, 0.0, 0.0]), "rhoval" : 0.5,
                     "jacobian" : np.array([0,0,0,0,0,0,0,0,0])}
                d["bxyz"][0] = np.sqrt(1.0/np.linspace(0.01, 1.0, 6))[i]
                init("B_TC", **d, desc=tag + str(i))
                init("gc", **mrk, desc=tag + str(i))
                init("plasma_1D", **pls, desc=tag + str(i))

    def run_classical(self):
        """Run classical transport test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_classical_go+"0"):
            warnings.warn("Results already present: Test classical transport")
            return
        for tag in [PhysTest.tag_classical_go, PhysTest.tag_classical_gcf,
                    PhysTest.tag_classical_gca]:
            i = 0
            while tag + str(i) in self.ascot.data.bfield.ls(show=False):
                self._activateinputs(tag+str(i))
                self._runascot(tag+str(i))
                i += 1

    def check_classical(self):
        """Check classical transport test.
        """
        nscan = 0
        while PhysTest.tag_classical_go + str(nscan) in \
              self.ascot.data.bfield.ls(show=False):
            nscan += 1

        # Initialize plots
        fig = a5plt.figuredoublecolumn(3/2)
        ax = fig.add_subplot(1,1,1)

        # Numerical values
        ndim = 2 # Diffusion happens on 2D plane
        bnorm = np.zeros((nscan,)) * unyt.T
        Dgo   = np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgcf  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgca  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
        for i in range(nscan):
            run_go  = self.ascot.data[PhysTest.tag_classical_go  + str(i)]
            run_gcf = self.ascot.data[PhysTest.tag_classical_gcf + str(i)]
            run_gca = self.ascot.data[PhysTest.tag_classical_gca + str(i)]

            yi, zi, ti = run_go.getstate("y", "z", "mileage", state="ini")
            yf, zf, tf = run_go.getstate("y", "z", "mileage", state="end")
            Dgo[i] = np.mean( ( (yi - yf)**2 + (zi - zf)**2 ) / (tf - ti) ) \
                / (2*ndim)
            yi, zi, ti = run_gcf.getstate("y", "z", "mileage", state="ini")
            yf, zf, tf = run_gcf.getstate("y", "z", "mileage", state="end")
            Dgcf[i] = np.mean( ( (yi - yf)**2 + (zi - zf)**2 ) / (tf - ti) ) \
                / (2*ndim)
            yi, zi, ti = run_gca.getstate("y", "z", "mileage", state="ini")
            yf, zf, tf = run_gca.getstate("y", "z", "mileage", state="end")
            Dgca[i] = np.mean( ( (yi - yf)**2 + (zi - zf)**2 ) / (tf - ti) ) \
                / (2*ndim)

            bnorm[i] = np.sqrt(np.sum(run_go.bfield.read()["bxyz"]**2))

        # Analytical
        clog = 13.4
        ekin = run_go.getstate("ekin")[0]
        self.ascot.input_init(run=run_go.get_qid(), bfield=True, plasma=True)
        ne, Te = self.ascot.input_eval(
            6.2*unyt.m, 0*unyt.deg, 0*unyt.m, 0*unyt.s, "ne", "te")
        self.ascot.input_free()
        collfreq = physlib.collfreq_ie(unyt.mp, unyt.e, ne, Te, clog)
        rhog = physlib.gyrolength(unyt.mp, 1*unyt.e, ekin, 0.0, bnorm).to("m")
        Dana = collfreq * rhog**2 / 2

        # Plotting
        ax.scatter(1/bnorm**2, Dgo)
        ax.scatter(1/bnorm**2, Dgcf)
        ax.scatter(1/bnorm**2, Dgca)
        ax.plot(1/bnorm**2, Dana, color="black")
        ax.set_xlabel(r"$1/B^{2}$ [1/T$^2$]")
        ax.set_ylabel(r"Diffusion [m$^2$/s]")

        print("Test classical transport:")
        passed = True
        k0 = np.polyfit(1/bnorm**2, Dana, 1)[0]
        k1 = np.polyfit(1/bnorm**2, Dgo,  1)[0]
        k2 = np.polyfit(1/bnorm**2, Dgcf, 1)[0]
        k3 = np.polyfit(1/bnorm**2, Dgca, 1)[0]

        print(" slope (expected): %1.3f" % k0)
        f = ""
        if(np.abs(k0-k1) > 1e-2): passed=False; f = "(FAILED)"
        print("               GO: %1.3f %s" % (k1, f))
        f = ""
        if(np.abs(k0-k1) > 1e-2): passed=False; f = "(FAILED)"
        print("              GCF: %1.3f %s" % (k2, f))
        f = ""
        if(np.abs(k0-k1) > 1e-2): passed=False; f = "(FAILED)"
        print("              GCA: %1.3f %s" % (k3, f))

        return passed

    def init_neoclassical(self):
        """Initialize data for the neoclassical transport test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_neoclassical_go+"0"):
            warnings.warn("Inputs already present: Test neoclass. transport")
            return
        init = self.ascot.data.create_input

        # Ion densities to be scanned
        ni = np.array([8.9e17, 8.9e18, 3.0e19, 7.0e20, 2.0e21, 2.0e22])

        # Options (some parameters are changed between the scans)
        optgo = Opt.get_default()
        optgo.update({
            "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "ENDCOND_SIMTIMELIM" : 1, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_COULOMB_COLLISIONS" : 1
        })
        optgcf = copy.deepcopy(optgo)
        optgcf.update({
            "SIM_MODE" : 2
        })
        optgca = copy.deepcopy(optgcf)
        optgca.update({
            "ENABLE_ADAPTIVE" : 1, "ADAPTIVE_MAX_DRHO" : 0.1,
            "ADAPTIVE_TOL_ORBIT" : 1e-8, "ADAPTIVE_TOL_CCOL" : 1e-1,
            "ADAPTIVE_MAX_DPHI" : 10, "FIXEDSTEP_USERDEFINED" : 1e-10
        })

        # Marker input consists of electrons
        mrk = Marker.generate("gc", n=100, species="electron")
        mrk["r"][:]      = 7.2
        mrk["phi"][:]    = 0.0
        mrk["z"][:]      = 0.0
        mrk["energy"][:] = 1e3
        mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
        mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)

        # Plasma consisting of protons only to avoid e-e collisions
        pls = init("plasma_flat", temperature=1e3, dryrun=True)
        pls["edensity"][:] = 1

        simtime = np.array([10e-4, 8e-4, 6e-4, 4e-4, 2e-4, 1e-4])
        for i in range(ni.size):

            # Adjust simulation time and time step as the density changes
            # (otherwise simulations with low density would take very long)
            optgo.update({
                "ENDCOND_LIM_SIMTIME" : simtime[i],
                "FIXEDSTEP_USERDEFINED" : 0.2e-9
            })
            optgcf.update({
                "ENDCOND_LIM_SIMTIME" : simtime[i],
                "FIXEDSTEP_USERDEFINED" : 0.2e-9
            })
            optgca.update({
                "ENDCOND_LIM_SIMTIME" : simtime[i]
            })
            init("opt", **optgo,  desc=PhysTest.tag_neoclassical_go  + str(i))
            init("opt", **optgcf, desc=PhysTest.tag_neoclassical_gcf + str(i))
            init("opt", **optgca, desc=PhysTest.tag_neoclassical_gca + str(i))

            pls["idensity"][:] = ni[i]
            for tag in [PhysTest.tag_neoclassical_go,
                        PhysTest.tag_neoclassical_gcf,
                        PhysTest.tag_neoclassical_gca]:
                init("gc", **mrk, desc=tag + str(i))
                init("plasma_1D", **pls, desc=tag + str(i))
                init("bfield_analytical_iter_circular", desc=tag + str(i))

    def run_neoclassical(self):
        """Run neoclassical transport test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_neoclassical_go+"0"):
            warnings.warn("Results already present: Test neoclass. transport")
            return
        for tag in [PhysTest.tag_neoclassical_go, PhysTest.tag_neoclassical_gcf,
                    PhysTest.tag_neoclassical_gca]:
            i = 0
            while tag + str(i) in self.ascot.data.bfield.ls(show=False):
                self._activateinputs(tag+str(i))
                self._runascot(tag+str(i))
                i += 1

    def check_neoclassical(self):
        """Check neoclassical transport test.
        """
        nscan = 0
        items = self.ascot.data.bfield.ls(show=False)
        while PhysTest.tag_neoclassical_go + str(nscan) in items:
            nscan += 1

        # Evaluate Ti and R_omp(rho) as these are needed
        run_go = self.ascot.data[PhysTest.tag_neoclassical_go  + "0"]
        self.ascot.input_init(run=run_go.get_qid(), bfield=True, plasma=True)
        Ti         = self.ascot.input_eval(
            6.2*unyt.m, 0*unyt.deg, 0*unyt.m, 0*unyt.s, "ti1")
        rhoomp     = np.linspace(0, 1, 100)
        romp, zomp = self.ascot.input_rhotheta2rz(
            rhoomp, 0*unyt.rad, 0*unyt.rad, 0*unyt.s)
        self.ascot.input_free()

        # Initialize plots
        fig = a5plt.figuredoublecolumn(3/2)
        ax = fig.add_subplot(1,1,1)
        ax.set_xscale("log")
        ax.set_yscale("log")

        # Numerical values
        ni    = np.zeros((nscan,)) / unyt.m**3
        Dgo   = np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgoerr= np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgcf  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgcferr= np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgca  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
        Dgcaerr= np.zeros((nscan,)) * unyt.m**2 / unyt.s
        for i in range(nscan):
            run_go  = self.ascot.data[PhysTest.tag_neoclassical_go  + str(i)]
            run_gcf = self.ascot.data[PhysTest.tag_neoclassical_gcf + str(i)]
            run_gca = self.ascot.data[PhysTest.tag_neoclassical_gca + str(i)]

            ri, ti = run_go.getstate("rho", "mileage", state="ini")
            rf, tf = run_go.getstate("rho", "mileage", state="end")
            ri = np.interp(ri, rhoomp, romp) * unyt.m
            rf = np.interp(rf, rhoomp, romp) * unyt.m
            Dgo[i] = 0.5 * np.mean( (rf - ri)**2 / (tf - ti) )
            Dgoerr[i] = np.sqrt( (0.5 * np.var( (rf - ri)**2 / (tf - ti) )) \
                                 / ri.size )

            ri, ti = run_gcf.getstate("rho", "mileage", state="ini")
            rf, tf = run_gcf.getstate("rho", "mileage", state="end")
            ri = np.interp(ri, rhoomp, romp) * unyt.m
            rf = np.interp(rf, rhoomp, romp) * unyt.m
            Dgcf[i] = 0.5 * np.mean( (rf - ri)**2 / (tf - ti) )
            Dgcferr[i] = np.sqrt( (0.5 * np.var( (rf - ri)**2 / (tf - ti) )) \
                                  / ri.size )

            ri, ti = run_gca.getstate("rho", "mileage", state="ini")
            rf, tf = run_gca.getstate("rho", "mileage", state="end")
            ri = np.interp(ri, rhoomp, romp) * unyt.m
            rf = np.interp(rf, rhoomp, romp) * unyt.m
            Dgca[i] = 0.5 * np.mean( (rf - ri)**2 / (tf - ti) )
            Dgcaerr[i] = np.sqrt( (0.5 * np.var( (rf - ri)**2 / (tf - ti) )) \
                                  / ri.size )

            ni[i] = run_go.plasma.read()["idensity"][0, 0]

        r0    = run_go.getstate("r")[0]
        axisr = run_go.bfield.read()["raxis"][0] * unyt.m
        eps   = (r0 - axisr) / axisr
        qfac  = 1.7 # Safety factor at r0 was verified numerically
        bnorm = run_go.bfield.read()["bphi0"][0] * unyt.T

        ekin   = run_go.getstate("ekin")[0]
        omegat = physlib.bouncefrequency(unyt.me, ekin, r0, axisr, qfac)
        rhog   = physlib.gyrolength(unyt.me, 1*unyt.e, ekin, 0.0, bnorm).to("m")

        clog     = 15
        density  = np.power( 10, np.linspace(np.log10(ni[0]) - 1,
                                             np.log10(ni[-1]) + 1, 50) )
        density /= unyt.m**3
        collfreq = physlib.collfreq_ie(unyt.mp, unyt.e, density, Ti, clog) \
            * ( unyt.mp / unyt.me )
        veff = collfreq / omegat
        # Add intermediate values needed for plotting a continuous curve
        veff = np.append(veff, [1, np.power(eps, 3.0/2.0)])
        veff.sort()

        Dps = qfac**2 * veff * omegat * rhog**2 / 2
        Dp  = 0.5 * qfac**2 * omegat * rhog**2 * np.ones(veff.shape)
        Db  = np.power(eps, -3.0/2.0) * Dps

        # x coordinate for plotting the numerical coefficients
        collfreq = physlib.collfreq_ie(unyt.mp, unyt.e, ni, Ti, clog) \
            * ( unyt.mp / unyt.me )
        veff_x = collfreq / omegat

        i1 = np.nonzero(veff==np.power(eps, 3.0/2.0))[0][0]
        i2 = np.nonzero(veff==1)[0][0]
        ax.plot(veff[i2:],     Dps[i2:],    color="black")
        ax.plot(veff[i1:i2+1], Dp[i1:i2+1], color="black")
        ax.plot(veff[:i1+1],   Db[:i1+1],   color="black")

        ax.errorbar(veff_x, Dgo, yerr=Dgoerr,  linestyle="none", marker="*")
        ax.errorbar(veff_x*1.01, Dgcf, yerr=Dgcferr, linestyle="none",
                    marker="o")
        ax.errorbar(veff_x*1.09, Dgca, yerr=Dgcaerr, linestyle="none",
                    marker="^")

        ax.set_xlabel(r"Effective collisionality $\nu^*$")
        ax.set_ylabel(r"Diffusion [m$^2$/s]")

        ax.plot([1,1], [1e-6, 1e-1], color="gray")
        ax.plot([eps**(3.0/2), eps**(3.0/2)], [1e-6, 1e-1], color="gray")
        ax.text(10**(-1.9), 10**(-5.7), r"$\nu^*=\epsilon^{3/2}$", fontsize=10,
           bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
        ax.text(10**(-0.25), 10**(-5.7), r"$\nu^*=1$", fontsize=10,
               bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
        ax.text(10**(-2.5), 10**(-3.7), r"$D_{B}$",  fontsize=10)
        ax.text(10**(-0.8), 10**(-2.7), r"$D_{P}$",  fontsize=10)
        ax.text(10**(0.4),  10**(-1.8), r"$D_{PS}$", fontsize=10)

        ax.set_xlim(1e-4, 1e2)
        ax.set_ylim(1e-6, 1e-1)

        f = "(FAILED)"
        print("Test neoclassical transport:")
        print("                   GO       GCF      GCA      analytical")
        passed = True
        idx = veff_x <= np.power(eps, 3.0/2.0)
        k0 = np.polyfit(veff[:i1+1], Db[:i1+1], 1)[0]
        k1 = np.polyfit(veff_x[idx], Dgo[idx],  1)[0]
        k2 = np.polyfit(veff_x[idx], Dgcf[idx], 1)[0]
        k3 = np.polyfit(veff_x[idx], Dgca[idx], 1)[0]
        f = ""
        if np.amax(np.abs(np.array([k1,k2,k3]) - k0)) > 2e-2:
            f = "(FAILED)"
            passed=False
        print("  Banana regime    %1.2e %1.2e %1.2e %1.2e %s"
              % (k1, k2, k3, k0, f))

        idx = np.logical_and(veff_x >= np.power(eps, 3.0/2.0), veff_x <=1)
        #k0 = np.polyfit(veff[i1:i2+1], Dp[i1:i2+1], 1)[0]
        k0 = 0.0
        k1 = np.polyfit(veff_x[idx],   Dgo[idx],    1)[0]
        k2 = np.polyfit(veff_x[idx],   Dgcf[idx],   1)[0]
        k3 = np.polyfit(veff_x[idx],   Dgca[idx],   1)[0]
        f = ""
        if np.amax(np.abs(np.array([k1,k2,k3]) - k0)) > 3e-3:
            f = "(FAILED)"
            passed=False
        print("  Plateau regime   %1.2e %1.2e %1.2e %1.2e %s"
              % (k1, k2, k3, k0, f))

        idx = veff_x >=1
        k0 = np.polyfit(veff[i2:],   Dps[i2:],  1)[0]
        k1 = np.polyfit(veff_x[idx], Dgo[idx],  1)[0]
        k2 = np.polyfit(veff_x[idx], Dgcf[idx], 1)[0]
        k3 = np.polyfit(veff_x[idx], Dgca[idx], 1)[0]
        f = ""
        if np.amax(np.abs(np.array([k1,k2,k3]) - k0)) > 5e-4:
            f = "(FAILED)"
            passed=False
        print("  Pfirsch-Schlüter %1.2e %1.2e %1.2e %1.2e %s"
              % (k1, k2, k3, k0, f))

        return passed

    def init_boozer(self):
        """Initialize data for the Boozer transformation test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_boozer):
            warnings.warn("Inputs already present: Test Boozer transformation")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "SIM_MODE" : 4, "ENDCOND_SIMTIMELIM" : 1,
            "ENDCOND_MAX_MILEAGE" : 1e2/3e8, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
            "ORBITWRITE_INTERVAL" : 1e-1/3e8, "ORBITWRITE_NPOINT" : 10002
        })

        # Use field line markers
        mrk = Marker.generate("fl", n=1)
        mrk["r"][:]      = 8.0
        mrk["phi"][:]    = 0
        mrk["z"][:]      = 0
        mrk["pitch"][:]  = 1.0

        # Magnetic field is just some tokamak with Boozer data
        b = init("bfield_analytical_iter_circular", dryrun=True)
        b.update({"rmin" : 4, "rmax" : 8.5, "nr" : 120, "zmin" : -4,
                  "zmax" : 4, "nz" : 200})

        bphi = [5.3, -5.3, 5.3, -5.3]
        psimult = [200, 200, -200, -200]
        for i in range(4):
            b.update({"bphi0" : bphi[i], "psimult" : psimult[i]})
            out = B_2DS.convert_B_GS(**b)
            qid = init("B_2DS", **out, desc=PhysTest.tag_boozer+str(i))
            qid = self.ascot.data.bfield[qid].get_qid()
            self.ascot.input_init(bfield=qid)
            qid = init("boozer_tokamak", desc=PhysTest.tag_boozer+str(i),
                       nint=100000)
            self.ascot.input_free(bfield=True)

            init("fl", **mrk, desc=PhysTest.tag_boozer+str(i))
            init("opt", **opt, desc=PhysTest.tag_boozer+str(i))

        # This perturbation is used only in the post-processing
        mhd = {"nmode" : 1, "nmodes" : np.array([2]), "mmodes" : np.array([3]),
               "amplitude" : np.array([1.0]), "omega" : np.array([1.0]),
               "phase" : np.array([0.0]), "nrho" : 100, "rhomin" : 0.0,
               "rhomax" : 1.0}
        rhogrid = np.linspace(mhd["rhomin"], mhd["rhomax"], mhd["nrho"])
        alpha   = np.exp( -(rhogrid-0.85)**2/0.1 )
        phi     = alpha*0
        mhd["phi"]   = np.tile(phi, (mhd["nmode"],1)).T
        mhd["alpha"] = np.tile(alpha, (mhd["nmode"],1)).T
        init("MHD_STAT", **mhd, desc=PhysTest.tag_boozer+"0")

    def run_boozer(self):
        """Run Boozer transformation test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_boozer + "0"):
            warnings.warn("Results already present: Test Boozer transformation")
            return
        for i in range(4):
            self._activateinputs(PhysTest.tag_boozer + str(i))
            self._runascot(PhysTest.tag_boozer + str(i))

    def check_boozer(self):
        """Check Boozer transformation test.
        """
        run1 = self.ascot.data[PhysTest.tag_boozer+"0"]
        run2 = self.ascot.data[PhysTest.tag_boozer+"1"]
        run3 = self.ascot.data[PhysTest.tag_boozer+"2"]
        run4 = self.ascot.data[PhysTest.tag_boozer+"3"]

        # Initialize plots
        fig = a5plt.figuredoublecolumn(3/2)
        ax1 = fig.add_subplot(3,2,1)
        ax2 = fig.add_subplot(3,2,3)
        ax3 = fig.add_subplot(3,2,5)
        ax4 = fig.add_subplot(3,2,2)
        ax5 = fig.add_subplot(3,2,4)

        ax3.set_xlabel("Boozer theta [rad]")
        ax1.set_ylabel("Error in Bpol [T]")
        ax2.set_ylabel("Error in Bphi [T]")
        ax3.set_ylabel("Jacobian")
        ax4.set_ylabel("q")
        ax5.set_ylabel("q")
        ax5.set_xlabel("Orbit parameter s")

        colors = ["C0", "C1", "C2", "C3"]
        ip_out   = [None] * 4
        bphi_out = [None] * 4
        bpol_err = [None] * 4
        bphi_err = [None] * 4
        jac_err  = [None] * 4
        q_err    = [None] * 4
        for i, run in enumerate([run1, run2, run3, run4]):
            self.ascot.input_init(run=run.get_qid(), bfield=True, boozer=True,
                                  mhd=True)
            r, phi, z, t, pol = run.getorbit("r", "phi", "z", "time", "theta")
            theta, zeta, alpha, jacb2 = self.ascot.input_eval(
                r, phi, z, t, "theta", "zeta", "alphaeig", "bjacxb2")

            # Results are plotted as a function of poloidal angle so find the
            # "discontinuity" points at 0/2pi.
            idx = np.nonzero(np.abs(np.diff(theta)) > np.pi)[0]
            theta = theta.v
            zeta  = zeta.v

            # Numerical safety factor from field lines and field data
            dz = np.diff(zeta)
            dz[dz > np.pi]  = dz[dz > np.pi]  - 2*np.pi
            dz[dz < -np.pi] = dz[dz < -np.pi] + 2*np.pi
            dt = np.diff(theta)
            dt[dt > np.pi]  = dt[dt > np.pi]  - 2*np.pi
            dt[dt < -np.pi] = dt[dt < -np.pi] + 2*np.pi
            qfac = dz/dt

            rho = run.getstate("rho", state="ini")
            q, I, g = self.ascot.input_eval_safetyfactor(rho)

            # Evaluate gradients and field components so we can compare those
            a, b, c = self.ascot.input_eval(
                r, phi, z, t, "dpsidr (bzr)", "dpsidphi (bzr)", "dpsidz (bzr)")
            gradpsi = np.array([a, b/r, c]).T
            a, b, c = self.ascot.input_eval(
                r, phi, z, t, "dthetadr", "dthetadphi", "dthetadz")
            gradtheta = np.array([a, b/r, c]).T
            a, b, c = self.ascot.input_eval(
                r, phi, z, t, "dzetadr", "dzetadphi", "dzetadz")
            gradzeta = np.array([a, b/r, c]).T

            br, bphi, bz = self.ascot.input_eval(
                r, phi, z, t, "br", "bphi", "bz")
            self.ascot.input_free()

            # Magnetic field vector from Boozer coordinates
            bvec = -( np.mean(qfac) * np.cross(gradpsi, gradtheta)
                      - np.cross(gradpsi, gradzeta) )

            idx = np.nonzero(np.abs(np.diff(theta)) > np.pi)[0]
            dbpol = np.sqrt((bvec[:,0] - br.v)**2 + (bvec[:,2] - bz.v)**2)
            dbphi = bvec[:,1] - bphi.v
            bnorm = br**2 + bz**2 + bphi**2

            # Jacobian times B^2
            jacb2_bzr = I + g*q

            ip_out[i]   = 1 - 2 * (bz[0] > 0)
            bphi_out[i] = 1 - 2 * (bphi[0] < 0)
            bpol_err[i] = np.amax(np.abs(dbpol[1:idx[-1]]))
            bphi_err[i] = np.amax(np.abs(dbphi[1:idx[-1]]))
            jac_err[i]  = np.amax(np.abs(jacb2-jacb2_bzr))
            q_err[i]    = np.amax(np.abs(q-qfac))

            j0 = 1
            j = idx[-1]
            ax1.plot(theta[j0:j], dbpol[j0:j], color=colors[i])
            ax2.plot(theta[j0:j], dbphi[j0:j], color=colors[i])
            ax3.plot(theta[j0:j], np.abs(jacb2[j0:j]), color=colors[i])
            ax3.plot(theta[np.array([j0,j])], np.abs([jacb2_bzr, jacb2_bzr]),
                     color=colors[i], ls='--')
            if qfac[0] > 0:
                ax4.plot(qfac, color=colors[i])
                ax4.plot([0,qfac.size-1], [q,q], color=colors[i], ls='--')
            else:
                ax5.plot(qfac, color=colors[i])
                ax5.plot([0,qfac.size-1], [q,q], color=colors[i], ls='--')

        passed = True
        print("Test Boozer coordinate mapping:")
        print("Bphi Ip | Delta Bpol  Delta Bphi  Delta Jacobian  Delta q")
        for i in range(4):
            fail = ""
            if bpol_err[i] > 1e-15 or bphi_err[i] > 1e-3 or jac_err[i] > 5e-3\
               or q_err[i] > 1e-4:
                fail = "(FAILED)"
                passed = False
            print(" %2d  %2d    %.3e   %.3e   %.3e       %.3e %s" %
                  (ip_out[i], bphi_out[i], bpol_err[i],
                   bphi_err[i], jac_err[i], q_err[i], fail))
        return passed

    def init_mhd(self):
        """Initialize data for the MHD test.

        https://link.springer.com/article/10.1007/s41614-018-0022-9
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_mhd_go):
            warnings.warn("Inputs already present: Test MHD")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_SIMTIMELIM" : 1,
            "ENDCOND_LIM_SIMTIME" : 1e-5, "ENABLE_ORBIT_FOLLOWING" : 1,
            "ENABLE_MHD" : 1, "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
            "ORBITWRITE_INTERVAL" : 1e-9, "ORBITWRITE_NPOINT" : 10000
        })
        init("opt", **opt, desc=PhysTest.tag_mhd_go)
        opt.update({
            "SIM_MODE" : 2, "FIXEDSTEP_USERDEFINED" : 1e-12
        })
        init("opt", **opt, desc=PhysTest.tag_mhd_gcf)
        opt.update({
            "ENABLE_ADAPTIVE" : 1, "FIXEDSTEP_USERDEFINED" : 1e-11,
            "ADAPTIVE_TOL_ORBIT" : 1e-11, "ADAPTIVE_MAX_DRHO" : 0.1,
            "ADAPTIVE_MAX_DPHI" : 10
        })
        init("opt", **opt, desc=PhysTest.tag_mhd_gca)

        # Use field line markers
        mrk = Marker.generate("gc", n=2, species="electron")
        mrk["r"][:]      = np.array([7.0, 8.0])
        mrk["phi"][:]    = 0
        mrk["z"][:]      = 0
        mrk["pitch"][:]  = np.array([0.4, 0.9])
        mrk["energy"][:] = 1e6

        # ITER-like field with boozer data and field line markers
        for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
                    PhysTest.tag_mhd_gca]:
            qid = init("bfield_analytical_iter_circular", splines=True,
                       desc=tag)
            init("gc", **mrk, desc=tag)

        qid = self.ascot.data.bfield[qid].get_qid()
        self.ascot.input_init(bfield=qid)
        bzr = init("boozer_tokamak", nint=100000, dryrun=True)
        for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
                    PhysTest.tag_mhd_gca]:
            init("Boozer", **bzr, desc=tag)

        mhd = {
            "nmode":1, "nmodes":np.array([2]), "mmodes":np.array([3]),
            "amplitude":np.array([1e-2]), "omega":np.array([1e6]),
            "phase":np.array([np.pi/4]), "nrho":100, "rhomin":0.1, "rhomax":0.99
        }

        rhogrid = np.linspace(mhd["rhomin"], mhd["rhomax"], mhd["nrho"])
        mhd["alpha"] = np.tile(
            np.exp( -(rhogrid-0.85)**2/0.1 ), (mhd["nmode"],1)).T
        mhd = self.ascot.data.create_input(
            "mhd consistent potentials", which="Phi", mhd=mhd, dryrun=True)
        self.ascot.input_free(bfield=True)
        for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
                    PhysTest.tag_mhd_gca]:
            self.ascot.data.create_input("MHD_STAT", **mhd, desc=tag)

    def run_mhd(self):
        """Run MHD test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_mhd_go):
            warnings.warn("Results already present: Test MHD")
            return
        for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
                    PhysTest.tag_mhd_gca]:
            self._activateinputs(tag)
            self._runascot(tag)

    def check_mhd(self):
        """Check MHD test.
        """
        run_go  = self.ascot.data[PhysTest.tag_mhd_go]
        run_gcf = self.ascot.data[PhysTest.tag_mhd_gcf]
        run_gca = self.ascot.data[PhysTest.tag_mhd_gca]

        mhd = run_go.mhd.read()
        self.ascot.input_init(run=run_go.get_qid(),
                              bfield=True, boozer=True, mhd=True)

        # Initialize plots
        fig = a5plt.figuredoublecolumn(3/2)
        ax1a = fig.add_subplot(3,1,1)
        ax1b = ax1a.twinx()
        ax2a = fig.add_subplot(3,1,2)
        ax2b = ax2a.twinx()
        ax3a = fig.add_subplot(3,1,3)
        ax3b = ax3a.twinx()

        ax3a.set_xlabel("Mileage [s]")
        ax2a.set_ylabel("Relative error")

        err1 = [None] * 3
        err2 = [None] * 3
        c = ["C0", "C1", "C2"]
        for i, run in enumerate([run_go, run_gcf, run_gca]):
            ekin, charge, ptor, bphi, r, phi, z, t, ids = run.getorbit(
                "ekin", "charge", "ptor", "bphi", "r", "phi", "z",
                "time", "ids")
            Phi, alpha, br0, bphi0, bz0, er, ephi, ez = self.ascot.input_eval(
                r, phi, z, t, "phieig", "alphaeig", "br", "bphi", "bz",
                "mhd_er", "mhd_ephi", "mhd_ez")
            epar  = ( er * br0 + ephi * bphi0 + ez * bz0 ) / \
                np.sqrt(br0**2 + bphi0**2 + bz0**2)
            enorm = np.sqrt(er**2 + ephi**2 + ez**2)

            ctor = ptor + r * charge * alpha * bphi
            H = ekin + Phi * charge
            P = (mhd["omega"]/unyt.s) * ctor / mhd["nmodes"]
            K = H - P

            H.convert_to_mks()
            P.convert_to_mks()
            K.convert_to_mks()

            idx = ids == 1
            ax1a.plot(t[idx], epar[idx], color=c[i])
            ax1b.plot(t[idx], enorm[idx], ls="--", color=c[i])
            err_e1 = np.amax(np.abs((epar/enorm)[idx]))
            err = K[ids==1]/K[ids==1][0] - 1
            ax2b.plot(t[idx], H[idx] - H[idx][0], ls="--", color=c[i])
            ax2b.plot(t[idx], -P[idx] + P[idx][0], ls=":", color=c[i])
            ax2a.plot(t[idx], K[idx] - K[idx][0], color=c[i])
            err1[i] = np.amax(np.abs(err))

            idx = ids == 2
            ax1a.plot(t[idx], epar[idx], color=c[i])
            ax1b.plot(t[idx], enorm[idx], ls="--", color=c[i])
            err_e2 = np.amax(np.abs((epar/enorm)[idx]))
            err = K[idx]/K[idx][0] - 1
            ax3b.plot(t[idx], H[idx] - H[idx][0], ls="--", color=c[i])
            ax3b.plot(t[idx], -P[idx] + P[idx][0], ls=":", color=c[i])
            ax3a.plot(t[idx], K[idx] - K[idx][0], color=c[i])
            err2[i] = np.amax(np.abs(err))

        self.ascot.input_free()

        passed = True
        print("Test MHD:")
        fail = ""
        if err_e1 > 2e-6 or err_e2 > 2e-5:
            fail = "(FAILED)"
            passed = False
        print("Error in Epar: %e %e %s" % (err_e1, err_e2, fail))

        fail = ""
        if err1[0] > 1e-6 or err2[0] > 1e-3:
            fail = "(FAILED)"
            passed = False
        print("Error in H - P (GO):  %e %e %s" % (err1[0], err2[0], fail))

        fail = ""
        if err1[1] > 1e-7 or err2[1] > 1e-3:
            fail = "(FAILED)"
            passed = False
        print("Error in H - P (GCF): %e %e %s" % (err1[1], err2[1], fail))

        fail = ""
        if err1[2] > 1e-9 or err2[2] > 1e-6:
            fail = "(FAILED)"
            passed = False
        print("Error in H - P (GCA): %e %e %s" % (err1[2], err2[2], fail))
        return passed

    def init_atomic(self):
        """Initialize data for the atomic reaction test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_atomic_cx):
            warnings.warn("Inputs already present: Test atomic reaction")
            return
        init = self.ascot.data.create_input

        # Options
        opt = Opt.get_default()
        opt.update({
            "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
            "ENDCOND_SIMTIMELIM" : 1, "ENDCOND_LIM_SIMTIME" : 2e-2,
            "ENABLE_ORBIT_FOLLOWING" : 1, "ENABLE_ATOMIC" : 1
        })
        opt.update({
            "FIXEDSTEP_USERDEFINED" : 1e-9, "ENDCOND_NEUTRALIZED" : 1
        })
        init("opt", **opt, desc=PhysTest.tag_atomic_cx)
        opt.update({
            "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_NEUTRALIZED" : 0,
            "ENDCOND_IONIZED" : 1, "ENDCOND_LIM_SIMTIME" : 2e-6
        })
        init("opt", **opt, desc=PhysTest.tag_atomic_ionz)

        mrk = Marker.generate("gc", n=1000, species="deuterium")
        mrk["r"][:]      = 6.2
        mrk["phi"][:]    = 0.0
        mrk["z"][:]      = 0.0
        mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)
        mrk["energy"][:] = 1e5
        mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
        init("gc", **mrk, desc=PhysTest.tag_atomic_cx)

        mrk = Marker.generate("gc", n=1000, species="deuterium")
        mrk["charge"][:] = 0
        mrk["r"][:]      = 6.2
        mrk["phi"][:]    = 0.0
        mrk["z"][:]      = 0.0
        mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)
        mrk["energy"][:] = 1e5
        mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
        init("gc", **mrk, desc=PhysTest.tag_atomic_ionz)

        for tag in [PhysTest.tag_atomic_cx, PhysTest.tag_atomic_ionz]:
            # Uniform magnetic field
            btc = {"bxyz":np.array([0,1.0,0]), "jacobian":np.zeros((3,3)),
                   "rhoval":0.5}
            init("B_TC", **btc, desc=tag)

            # Plasma and neutral data
            init("plasma_flat", anum=2, znum=1, mass=2.0135532,
                 density=1e20, temperature=1e3, desc=tag)
            init("neutral_flat", anum=2, znum=1, density=1e16, temperature=1e3,
                 desc=tag)

            # ADAS data
            init("import_adas", desc=tag)

    def run_atomic(self):
        """Run atomic reaction test.
        """
        if hasattr(self.ascot.data, PhysTest.tag_atomic_cx):
            warnings.warn("Results already present: Test atomic reaction")
            return
        for tag in [PhysTest.tag_atomic_cx, PhysTest.tag_atomic_ionz]:
            self._activateinputs(tag)
            self._runascot(tag)

    def check_atomic(self):
        """Check atomic reaction test.
        """
        run_cx   = self.ascot.data[PhysTest.tag_atomic_cx]
        run_ionz = self.ascot.data[PhysTest.tag_atomic_ionz]

        # Calculate mean free paths from the cross sections
        ma, anum, znum, r, phi, z, t, va = run_cx.getstate(
            "mass", "anum", "znum", "r", "phi", "z", "time", "vnorm",
            state="ini", mode="prt", ids=1)
        self.ascot.input_init(run=run_cx.get_qid(), bfield=True, plasma=True,
                              neutral=True, asigma=True)
        sigmacx  = self.ascot.input_eval_atomicsigma(
            ma, anum[0], znum[0], r, phi, z, t, va, ion=0, reaction=6)
        sigmabms = self.ascot.input_eval_atomicsigma(
            ma, anum[0], znum[0], r, phi, z, t, va, ion=0, reaction=7)
        self.ascot.input_free()
        mfp_cx0  = va/(sigmacx * 1e16/unyt.m**3)
        mfp_bms0 = va/(sigmabms * 1e20/unyt.m**3)

        # Mean free paths from the simulation
        vnorm1 = run_cx.getstate("vnorm", state="ini", mode="prt")
        vnorm2, mil = run_cx.getstate("vnorm", "mileage", state="end",
                                      mode="prt")
        mfp_cx = np.mean(mil*vnorm2)

        vnorm1 = run_ionz.getstate("vnorm", state="ini", mode="prt")
        vnorm2, mil = run_ionz.getstate("vnorm", "mileage", state="end",
                                        mode="prt")
        mfp_bms = np.mean(mil*vnorm2)

        passed = True
        print("Test atomic:")
        fail = ""
        if np.abs(mfp_cx/mfp_cx0 - 1) > 0.1:
            fail = "(FAILED)"
            passed = False
        print("Mean free path CX:  %.3e m (numerical) %.3e m (analytical) %s" %
              (mfp_cx, mfp_cx0[0,0], fail))
        fail = ""
        if np.abs(mfp_bms/mfp_bms0 - 1) > 0.1:
            fail = "(FAILED)"
            passed = False
        print("Mean free path BMS: %.3e m (numerical) %.3e m (analytical) %s" %
              (mfp_bms, mfp_bms0[0,0], fail))
        return passed

    def init_afsi(self):
        """Initialize data for AFSI test.
        """
        if hasattr(self.ascot.data.bfield, PhysTest.tag_afsi_thermal):
            warnings.warn("Inputs already present: Test AFSI")
            return

        bgs = self.ascot.data.create_input("bfield analytical iter circular",
                                           dryrun=True)

        # DT-plasma
        nrho  = 4
        rho   = np.array([0, 1, 1+1e-1, 10])
        edens = 2e20 * np.ones((nrho, 1))
        etemp = 1e4  * np.ones((nrho, 1))
        idens = 1e20 * np.ones((nrho, 2))
        itemp = 1e4  * np.ones((nrho, 1))

        edens[rho>1]   = 1
        idens[rho>1,:] = 1

        pls = {
            "nrho" : nrho, "nion" : 2, "rho" : rho,
            "anum" : np.array([2, 3]), "znum" : np.array([1, 1]),
            "mass" : np.array([2.014, 3.016]), "charge" : np.array([1, 1]),
            "edensity" : edens, "etemperature" : etemp,
            "idensity" : idens, "itemperature" : itemp}

        for tag in [PhysTest.tag_afsi_thermal, PhysTest.tag_afsi_beamthermal,
                    PhysTest.tag_afsi_beambeam]:
            self.ascot.data.create_input("B_GS", **bgs, desc=tag)
            self.ascot.data.create_input("plasma_1D", **pls, desc=tag)

    def run_afsi(self):
        """Run AFSI tests.
        """
        if hasattr(self.ascot.data, PhysTest.tag_afsi_thermal):
            warnings.warn("Results already present: Test AFSI")
            return

        rmin =  5.9; rmax = 6.1; nr = 1
        zmin = -0.1; zmax = 0.1; nz = 1
        self._activateinputs(PhysTest.tag_afsi_thermal)
        self.ascot.afsi.thermal(
            "DT_He4n",
            rmin, rmax, nr, zmin, zmax, nz,
            minphi=0, maxphi=2*np.pi, nphi=1, nmc=10**6,
            minppara=-1.0e-19, maxppara=1.0e-19, nppara=400,
            minpperp=0, maxpperp=1.0e-19, npperp=200)
        self.ascot.data.active.set_desc(PhysTest.tag_afsi_thermal)

        dist = self.ascot.data[PhysTest.tag_afsi_thermal].getdist("prod1")
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"), dist.abscissa("pperp"),
                               indexing="ij")
        vol = np.diff(dist.abscissa_edges("r")[:2]) \
            * np.diff(dist.abscissa_edges("z")[:2]) \
            * np.diff(dist.abscissa_edges("phi")[:2]) \
            * np.diff(dist.abscissa_edges("time")[:2]) \
            * np.diff(dist.abscissa_edges("charge")[:2])

        def maxwellian(mass):
            mass *= unyt.amu
            T = 1e4*unyt.eV
            nprt = 1e20*2*np.pi*6*0.2*0.2
            return nprt * ppe * np.exp( -(ppe**2 + ppa**2) / (2*mass*T) ) \
                / ( np.power(mass*T, 3/2.0) * vol * np.sqrt(2*np.pi) )

        dist._distribution[0,0,0,:,:,0,0] = maxwellian(2.014)

        self._activateinputs(PhysTest.tag_afsi_beamthermal)
        self.ascot.afsi.beamthermal("DT_He4n", dist, nmc=10**6)
        self.ascot.data.active.set_desc(PhysTest.tag_afsi_beamthermal)

        dist = self.ascot.data[PhysTest.tag_afsi_thermal].getdist("prod1")
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"), dist.abscissa("pperp"),
                               indexing="ij")
        dist._distribution[0,0,0,:,:,0,0] = maxwellian(2.014)

        dist1 = self.ascot.data[PhysTest.tag_afsi_thermal].getdist("prod1")
        dist1._distribution[0,0,0,:,:,0,0] = maxwellian(3.016)

        self._activateinputs(PhysTest.tag_afsi_beambeam)
        self.ascot.afsi.beambeam("DT_He4n", dist, dist1, nmc=10**6)
        self.ascot.data.active.set_desc(PhysTest.tag_afsi_beambeam)

    def check_afsi(self):
        """Check AFSI tests.
        """
        alphadist = self.ascot.data[PhysTest.tag_afsi_thermal].getdist(
            "prod1", exi=True, ekin_edges=np.linspace(0, 4.5e6, 100),
            pitch_edges=100)

        alphadist.integrate(
            time=np.s_[:], charge=np.s_[:], r=np.s_[:], z=np.s_[:],
            phi=np.s_[:])
        ekindist = alphadist.integrate(copy=True, pitch=np.s_[:])
        xidist   = alphadist.integrate(copy=True, ekin=np.s_[:])

        ekin  = alphadist.abscissa("ekin")
        pitch = alphadist.abscissa("pitch")

        emean = 3.537e6 + 1e4 * 2
        gauss = np.exp(-(ekin.v - emean)**2 / (4*1e4*emean*4.001/(4.001+1.008)))

        pdens = 1e20*unyt.m**(-3)*1e20*unyt.m**(-3)*1.1e-16*unyt.cm**3
        vol = 6*unyt.m*np.pi*2*0.2*unyt.m*0.2*unyt.m

        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)

        ekindist.plot(axes=ax1)
        ax1.plot(ekin, pdens * vol * gauss / np.trapz(gauss,ekin))
        xidist.plot(axes=ax2)
        ax2.plot(pitch, 0.5 * pdens * vol * np.ones(pitch.shape))

        eana  = pdens * vol * gauss / np.trapz(gauss,ekin)
        xiana = 0.5 * pdens * vol * np.ones(pitch.shape)
        ekin1 = ekindist.distribution()
        xi1   = xidist.distribution()

        alphadist = self.ascot.data[PhysTest.tag_afsi_beamthermal].getdist(
            "prod1", exi=True, ekin_edges=np.linspace(0, 4.5e6, 100),
            pitch_edges=100)

        alphadist.integrate(
            time=np.s_[:], charge=np.s_[:], r=np.s_[:], z=np.s_[:],
            phi=np.s_[:])
        ekindist = alphadist.integrate(copy=True, pitch=np.s_[:])
        xidist   = alphadist.integrate(copy=True, ekin=np.s_[:])

        ekindist.plot(axes=ax1)
        xidist.plot(axes=ax2)
        ekin2 = ekindist.distribution()
        xi2   = xidist.distribution()

        alphadist = self.ascot.data[PhysTest.tag_afsi_beambeam].getdist(
            "prod1", exi=True, ekin_edges=np.linspace(0, 4.5e6, 100),
            pitch_edges=100)

        alphadist.integrate(
            time=np.s_[:], charge=np.s_[:], r=np.s_[:], z=np.s_[:],
            phi=np.s_[:])
        ekindist = alphadist.integrate(copy=True, pitch=np.s_[:])
        xidist   = alphadist.integrate(copy=True, ekin=np.s_[:])
        ekin3 = ekindist.distribution()
        xi3   = xidist.distribution()

        ekindist.plot(axes=ax1)
        xidist.plot(axes=ax2)

        passed = True
        print("Test AFSI:")
        err1 = np.square(np.subtract(eana, ekin1)).mean()
        err2 = np.square(np.subtract(xiana, xi1)).mean()
        if err1.v > 1e23 or err2.v > 2e35:
            passed = False
            print("Thermal %e %e (FAILED)" % (err1.v, err2.v))
        else:
            print("Thermal %e %e" % (err1.v, err2.v))

        err1 = np.square(np.subtract(eana, ekin2)).mean()
        err2 = np.square(np.subtract(xiana, xi2)).mean()
        if err1.v > 1e23 or err2.v > 2e35:
            passed = False
            print("Beam-Thermal %e %e (FAILED)" % (err1.v, err2.v))
        else:
            print("Beam-Thermal %e %e" % (err1.v, err2.v))

        err1 = np.square(np.subtract(eana, ekin3)).mean()
        err2 = np.square(np.subtract(xiana, xi3)).mean()
        if err1.v > 1e23 or err2.v > 2e35:
            passed = False
            print("Beam-Beam %e %e (FAILED)" % (err1.v, err2.v))
        else:
            print("Beam-Beam %e %e" % (err1.v, err2.v))

        return passed

    def init_bbnbi(self):
        """Initialize data for the BBNBI deposition test.
        """
        if hasattr(self.ascot.data.options, PhysTest.tag_bbnbi):
            warnings.warn("Inputs already present: Test BBNBI")
            return
        init = self.ascot.data.create_input
        bdata = {"bxyz":np.array([0.0,1.0,0.0]), "jacobian":1.0*np.zeros((3,3)),
                 "rhoval":0.5}
        init("B_TC", **bdata, desc=PhysTest.tag_bbnbi)

        # Hydrogen plasma
        nrho  = 4
        rho   = np.array([0, 1, 1+1e-8, 10])
        edens = 1e19 * np.ones((nrho, 1))
        etemp = 2e3  * np.ones((nrho, 1))
        idens = 1e19 * np.ones((nrho, 1))
        itemp = 2e3  * np.ones((nrho, 1))
        edens[rho>=1]   = 1
        idens[rho>=1,:] = 1

        pls = {
            "nrho" : nrho, "nion" : 1, "rho" : rho,
            "anum" : np.array([1]), "znum" : np.array([1]),
            "mass" : np.array([1.014]), "charge" : np.array([1]),
            "edensity" : edens, "etemperature" : etemp,
            "idensity" : idens, "itemperature" : itemp}
        init("plasma_1D", **pls, desc=PhysTest.tag_bbnbi)

        opt = Opt.get_default()
        opt.update({
            "ENABLE_DIST_5D":0,
            "DIST_MIN_R":4.3,       "DIST_MAX_R":8.3,      "DIST_NBIN_R":50,
            "DIST_MIN_PHI":0,       "DIST_MAX_PHI":360,    "DIST_NBIN_PHI":1,
            "DIST_MIN_Z":-2.0,      "DIST_MAX_Z":2.0,      "DIST_NBIN_Z":50,
            "DIST_MIN_PPA":-1.e-19, "DIST_MAX_PPA":1.e-19, "DIST_NBIN_PPA":100,
            "DIST_MIN_PPE":0,       "DIST_MAX_PPE":1.e-19, "DIST_NBIN_PPE":50,
            "DIST_MIN_TIME":0,      "DIST_MAX_TIME":1.0,   "DIST_NBIN_TIME":1,
            "DIST_MIN_CHARGE":0.5,  "DIST_MAX_CHARGE":1.5, "DIST_NBIN_CHARGE":1,
        })
        init("opt", **opt, desc=PhysTest.tag_bbnbi)

        halo_f = 0.15
        div_v  = 10e-3
        div_h  =  5e-3
        halo_v = 30e-3
        halo_h = 20e-3
        inj = Injector(
            ids=1, anum=1, znum=1, mass=1.0,
            energy=4e4, efrac=np.array([1.0,0.0,0.0]), power=1.0,
            divh=div_h, divv=div_v, divhalofrac=halo_f,
            divhaloh=halo_h, divhalov=halo_v, nbeamlet=1,
            beamletx=np.array([10.0]), beamlety=np.array([0.0]),
            beamletz=np.array([0.0]),
            beamletdx=np.array([-1.0]), beamletdy=np.array([0.0]),
            beamletdz=np.array([0.0]))
        init("nbi", **{"ninj":1, "injectors":[inj]}, desc=PhysTest.tag_bbnbi)

    def run_bbnbi(self):
        """Run BBNBI.
        """
        if hasattr(self.ascot.data, PhysTest.tag_bbnbi):
            warnings.warn("Results already present: Test BBNBI")
            return
        self._activateinputs(PhysTest.tag_bbnbi)
        self.ascot.data.nbi[PhysTest.tag_bbnbi].activate()
        subprocess.call(
            ["./../../build/bbnbi5", "--in=testascot.h5",
             "--d="+PhysTest.tag_bbnbi],
            stdout=subprocess.DEVNULL)
        self.ascot = Ascot(self.ascot.file_getpath())

    def check_bbnbi(self):
        """Check the BBNBI deposition test.
        """
        run  = self.ascot.data[PhysTest.tag_bbnbi]
        r, z, x, y, s = run.getstate("r", "z", "x", "y", "mileage", mode="prt")

        omega_x = np.mod(np.arctan2(y.v, x.v-10) + 2*np.pi, 2*np.pi)-np.pi
        omega_y = np.mod(np.arctan2(z.v, r.v-10) + 2*np.pi, 2*np.pi)-np.pi

        inj = run.nbi.read()["injectors"][0]
        div_h  = inj.divh
        div_v  = inj.divv
        halo_h = inj.divhaloh
        halo_v = inj.divhalov
        halo_f = inj.divhalofrac

        # The plasma parameters and this value are taken from "Neutral beam
        # stopping and emission in fusion plasmas I: deuterium beams",
        # H Anderson et al 2000 Plasma Phys. Control. Fusion 42 781, p. 10
        mfp = 1.0/(1e19*1.1e-13)

        def gaussian(omega, omegadiv):
            return np.exp( -(omega/omegadiv)**2 ) / omegadiv

        def profile(omega, div, halo_div):
            return halo_f * gaussian(omega, halo_div) / np.sqrt(np.pi) \
                + (1.0 - halo_f) * gaussian(omega, div) / np.sqrt(np.pi)

        def ionization(distance, mfp):
            mfp *= 1e6
            return np.exp(-distance/mfp) / mfp

        dist_edges = np.linspace(0,7,100)
        distance = ( dist_edges[1:] + dist_edges[:-1] ) / 2

        omega_edges = np.linspace(-0.04*np.pi/2, 0.04*np.pi/2, 181)
        omega = ( omega_edges[1:] + omega_edges[:-1] ) / 2

        mrk1 = np.histogram(omega_x, omega_edges, density=True)[0]
        mrk2 = np.histogram(omega_y, omega_edges, density=True)[0]
        mrk3 = np.histogram(s.v*1e6,  dist_edges, density=True)[0]

        ana1 = profile(omega, div_h, halo_h)
        ana2 = profile(omega, div_v, halo_v)
        ana3 = ionization(distance, mfp)

        p = scipy.optimize.curve_fit(profile, omega, mrk1, [div_h, halo_h])[0]
        fit1 = profile(omega, p[0], p[1])
        div_h_fit  = p[0]
        halo_h_fit = p[1]

        p = scipy.optimize.curve_fit(profile, omega, mrk2, [div_v, halo_v])[0]
        fit2 = profile(omega, p[0], p[1])
        div_v_fit  = p[0]
        halo_v_fit = p[1]

        p = scipy.optimize.curve_fit(ionization, distance, mrk3, [mfp])[0]
        fit3 = ionization(distance, p[0])
        mfp_fit = p[0]

        fig = plt.figure(figsize=(5,5))
        gs = GridSpec(3, 1, figure=fig, hspace=0.6)
        ax1 = fig.add_subplot(gs[0,0])
        a5plt.mesh1d(omega_edges, mrk1, axes=ax1)
        ax1.plot(omega, ana1, color="black")
        ax1.plot(omega, fit1, color="C3", ls="--")
        ax1.set_xlabel(r"Horizontal divergence $\omega$ [rad]")

        ax2 = fig.add_subplot(gs[1,0])
        a5plt.mesh1d(omega_edges, mrk2, axes=ax2)
        ax2.plot(omega, ana2, color="black")
        ax2.plot(omega, fit2, color="C3", ls="--")
        ax2.set_xlabel(r"Vertical divergence $\omega$ [rad]")

        ax3 = fig.add_subplot(gs[2,0])
        a5plt.mesh1d(dist_edges, mrk3, axes=ax3)
        ax3.plot(distance, ana3, color="black")
        ax3.plot(distance, fit3, color="C3", ls="--")
        ax3.set_xlabel(r"Time until ionized [$\mu$s]")

        ax2.set_ylabel("Probability density")

        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color="C0", lw=2, label="Markers"),
            Line2D([0], [0], color="black", lw=2, label="Analytical"),
            Line2D([0], [0], color="C3", lw=2, ls="--", label="Fit"),
        ]
        ax3.legend(handles=legend_elements, loc="upper right", frameon=False)

        passed = True
        if np.abs(div_h  - div_h_fit)  > 1e-3: passed = False
        if np.abs(div_v  - div_v_fit)  > 5e-3: passed = False
        if np.abs(halo_h - halo_h_fit) > 5e-3: passed = False
        if np.abs(halo_v - halo_v_fit) > 5e-3: passed = False
        if np.abs(mfp    - mfp_fit)    > 1e-7: passed = False

        failed = "" if passed else "(FAILED)"

        print("Test BBNBI:")
        print("          div_h     div_v     halo_h    halo_v    mfp")
        print("  Actual: %.3e %.3e %.3e %.3e %.3e" \
              % (div_h[0], div_v[0], halo_h[0], halo_v[0], mfp) )
        print("  Fitted: %.3e %.3e %.3e %.3e %.3e %s" \
              % (div_h_fit, div_v_fit, halo_h_fit, halo_v_fit, mfp_fit, failed))

        return passed

    def init_biosaw(self):
        """Dummy function as BioSaw test requires no initialization.
        """
        pass

    def run_biosaw(self):
        """Dummy function as BioSaw test requires no simulation.
        """
        pass

    def check_biosaw(self):
        """Check and run BioSaw test results.
        """
        # Define coil geometry (Circular coil centered at x=r0,y=0,z=0)
        nseg = 200
        r0   = 2.0
        rad  = 1.0
        coilxyz = np.zeros((nseg,3))
        coilxyz[:,0] = r0  + rad * np.cos(np.linspace(0, 2*np.pi, nseg))
        coilxyz[:,1] = 0.0 + rad * np.sin(np.linspace(0, 2*np.pi, nseg))
        coilxyz[:,2] = 0.0

        # Coordinates where field is evaluated
        npnt = 100
        r   = np.array([r0])
        phi = np.array([0.0])
        z   = np.linspace(-10.0, 10.0, npnt)

        # Calculate analytical and numerical fields
        Icoil = 1.0
        br, bphi, bz = self.ascot.biosaw.calculate(
            r, phi, z, coilxyz, current=Icoil)
        br   = np.squeeze(br)
        bphi = np.squeeze(bphi)
        bz   = np.squeeze(bz)

        br0   = np.zeros((npnt,))
        bphi0 = np.zeros((npnt,))
        bz0   = 0.5 * unyt.mu_0.v * rad**2 * Icoil \
            / np.power(rad**2 + z**2, 3.0/2)

        fig = plt.figure()
        gs = GridSpec(4, 1, figure=fig, hspace=0.9, top=0.95)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0])

        ax1.plot(bz, color="C0")
        ax1.plot(bz0, color="black", ls="--")
        ax2.plot(br**2  + bphi**2, color="C0")
        ax2.plot(br0**2 + bphi0**2, color="black", ls="--")

        passed = True
        print("Test BioSaw:")
        if any(br**2 + bphi**2 > 1e-40) or any(np.abs(bz - bz0) > 1e-10):
            print("- Field on axis of a cylindrical loop (FAILED)")
            passed = False
        else:
            print("- Field on axis of a cylindrical loop")

        # Tokamak field test
        coilxyz = np.zeros((nseg,3))
        coilxyz[:,0] = r0  + rad * np.cos(np.linspace(0, 2*np.pi, nseg))
        coilxyz[:,1] = 0.0
        coilxyz[:,2] = 0.0 + rad * np.sin(np.linspace(0, 2*np.pi, nseg))

        npnt = 100
        r   = np.linspace(r0-rad*0.8, r0+rad*0.8, npnt)
        phi = np.linspace(0, 2*np.pi, npnt+1)[:-1]
        z   = np.array([0.0])

        Icoil = -1.0
        _, bphi, _ = self.ascot.biosaw.calculate(
            r, phi, z, coilxyz, current=Icoil, revolve=1)

        brad  = np.squeeze(bphi[:,0,0])
        bphi  = np.squeeze(bphi[int(npnt/2),:,0])
        brad0 = -unyt.mu_0.v * npnt * Icoil / ( 2 * np.pi * r )
        bphi0 = bphi[0] * np.ones((npnt,))

        if any(np.abs(bphi - bphi0) > 1e-10) or \
           any(np.abs(brad - brad0) > 1e-7):
            print("- Tokamak field (FAILED)")
            passed = False
        else:
            print("- Tokamak field")

        ax3 = fig.add_subplot(gs[2,0])
        ax4 = fig.add_subplot(gs[3,0])
        ax3.plot(r, brad, color="C0")
        ax3.plot(r, brad0, color="black", ls="--")
        ax4.plot(phi, bphi, color="C0")
        ax4.plot(phi, bphi0, color="black", ls="--")

        ax1.set_xlabel("$z$ [m]")
        ax1.set_ylabel("$B_z$ [T]")
        ax2.set_xlabel("$z$ [m]")
        ax2.set_ylabel("$|B_{xy}|$ [T]")

        ax3.set_xlabel("$r$ [m]")
        ax3.set_ylabel("$B_phi$ [m]")
        ax4.set_xlabel("$\phi$ [rad]")
        ax4.set_ylabel("$B_phi$ [m]")
        return passed

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

if __name__ == '__main__':
    test = PhysTest()
    # Atomic does not work standalone yet
    failed = test.execute(
        init=True, run=True, check=True,
        tests=["elementary", "orbitfollowing", "gctransform", "ccoll",
               "classical", "neoclassical", "boozer", "mhd", "afsi",
               "bbnbi", "biosaw"])#, "atomic"])
    plt.show(block=False)
    if failed: raise Exception("Verification failed")
