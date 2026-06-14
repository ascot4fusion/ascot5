"""Verify that the code reproduces gyromotion and guiding center drifts
accurately.

The tests in this file are done with energetic electrons to verify that the
orbit-following physics work in the relativistic regime.
"""
import numpy as np
import unyt

import pytest

from a5py import Ascot, SimulationOptions, physlib
from a5py import plotting as a5plt


def test_drifts_nodrifts(ascot, inspect, plot):
    """Verify that the code reproduces correct gyroradius and gyrofrequency in
    a uniform magnetic field.
    """
    if not inspect:
        parameters = SimulationOptions(
                simulation_mode=1,
                enable_adaptive=False,
                timestep=1e-11,
                enable_orbit_following=True,
                collect_orbit=True,
                number_of_points_per_marker=202,
                interval=1e-11,
                activate_simulation_time_limits=True,
                max_mileage=2e-9,
            )

        ascot.data.create_guidingcentermarker(
            species="electron",
            ids=np.array([1, 2]),
            r=5*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=100e6*unyt.eV,
            pitch=0.5,
            gyroangle=0.0*unyt.rad,
            charge=np.array([1, -1])*unyt.e,
            )

        ascot.data.create_bfieldcartesian(
            bxyz=np.array([5., 0., 0.])*unyt.T,
            jacobian=np.full((3, 3), 0.)*unyt.T/unyt.m,
            axisrz=np.array([1., 0.])*unyt.m,
            rhoval=1.,
            )

        ascot.data.bfield.active.stage()
        ascot.data.marker.active.stage()
        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    def plotarrow(ax, ini, end, txt):
        """Plot line with an arrow
        """
        ax.annotate("", xytext=ini, xy=end,
                    arrowprops={"arrowstyle":"->", "color":"black"})
        ax.plot([ini[0], end[0]], [ini[1], end[1]], color="black")
        ax.annotate(txt, ini + .005)

    m, q, pitch, ekin, x0, y0 = \
        run.getstate("mass", "charge", "pitch", "ekin", "y", "z")
    x, y, time = run.getorbit("y", "z", "time", filter=[1])

    gyro_radius, gyro_frequency = {}, {}
    bnorm = np.sqrt(np.sum(run.bfield.bxyz**2)) * unyt.T
    gamma = physlib.Quantity.registry["gamma"].compute(mass=m, energy=ekin)
    gyro_radius["analytical"] = (
        physlib.Quantity.registry["gyrolength"].compute(
            mass=m, charge=q, gamma=gamma, pitch=pitch, bnorm=bnorm,
            ).to("m")
    )
    gyro_frequency["analytical"] = (
        physlib.Quantity.registry["gyrofrequency"].compute(
            mass=m, charge=q, bnorm=bnorm, gamma=gamma,
            ).to("rad/s")
    )
    gyro_radius["numerical"] = (
        np.mean( np.sqrt( (x - x0[0])**2 + (y - y0[0])**2 ) )
    )
    gyro_frequency["numerical"] = np.sum( np.sqrt( np.diff(x)**2 + np.diff(y)**2 ) ) \
        * unyt.rad / (gyro_radius["numerical"] * time[-1])

    if plot:
        fig, ax = a5plt.subplots()

        ax.set_title("Gyro-motion")
        #h2a.set_title("ExB drift")
        #h3a.set_title("Grad-B drift")

        ax.set_ylabel("z [m]")
        ax.set_xlabel("y [m]")

        ax.set_aspect("equal", adjustable="box")
        orbx = x - x0[0]
        orby = y - y0[0]

        anax   = gyro_radius["numerical"][0] * np.sin(np.linspace(0, 2*np.pi, 360))
        anay   = gyro_radius["numerical"][0] * np.cos(np.linspace(0, 2*np.pi, 360))
        fangle = (
            np.arctan2(orby[0], orbx[0]) * unyt.rad + gyro_frequency["analytical"][0] * time[-1]
        )
        fx = gyro_radius["numerical"][0] * np.cos(-fangle) # Minus sign because this we are
        fy = gyro_radius["numerical"][0] * np.sin(-fangle) # in LHS coordinate system.

        ax.plot(orbx, orby, linewidth=3)
        ax.plot(anax, anay, linestyle="--", color="black", alpha=0.7)
        ax.scatter(orbx[-1], orby[-1], marker="o", color="C0", zorder=3)
        ax.scatter(fx, fy, marker="x", color="black", zorder=4)
        a5plt.show()

        ## E x B drift
        # yi_go            = run_exbgo.getstate("z", state="ini")
        # yf_go, deltat_go = run_exbgo.getstate("z", "mileage", state="end")
        # x_go1, y_go1     = run_exbgo.getorbit("y", "z", ids=1)
        # x_go2, y_go2     = run_exbgo.getorbit("y", "z", ids=2)

        # yi_gc            = run_exbgc.getstate("z", state="ini")
        # yf_gc, deltat_gc = run_exbgc.getstate("z", "mileage", state="end")
        # x_gc1, y_gc1     = run_exbgc.getorbit("y", "z", ids=1)
        # x_gc2, y_gc2     = run_exbgc.getorbit("y", "z", ids=2)

        # # Analytical values
        # bvec  = run_exbgo.bfield.read()["bxyz"].ravel()
        # evec  = run_exbgo.efield.read()["exyz"].ravel()
        # v_ExB = np.cross(evec, bvec) / np.inner(bvec, bvec) * unyt.m / unyt.s
        # time  = run_exbgo.options.read()["ENDCOND_LIM_SIMTIME"] * unyt.s

        # # Numerical values
        # vgo_ExB = (yf_go - yi_go) / deltat_go
        # vgc_ExB = (yf_gc - yi_gc) / deltat_gc

        # # Plot
        # h2a.plot(x_go1, y_go1)
        # h2b.plot(x_go2, y_go2)
        # h2a.plot(x_gc1, y_gc1)
        # h2b.plot(x_gc2, y_gc2)

        # ini = np.array([x_gc1[0].v + .01, y_gc1[0].v])
        # end = np.array([x_gc1[0].v + .01, (y_gc1[0] + v_ExB[2]*deltat_gc[1]).v])
        # plotarrow(h2a, ini, end, r"$e^+$")
        # ini = np.array([x_gc2[0].v + .01, y_gc2[0].v])
        # end = np.array([x_gc2[0].v + .01, (y_gc2[0] + v_ExB[2]*deltat_gc[0]).v])
        # plotarrow(h2b, ini, end, r"$e^-$")

        # ## gradB drift
        # pitch, ekin, q, m = run_gradbgo.getstate(
        #     "pitch", "ekin", "charge", "mass", ids=1)
        # xi_go            = run_gradbgo.getstate("y", state="ini")
        # xf_go, deltat_go = run_gradbgo.getstate("y", "mileage", state="end")
        # x_go1, y_go1     = run_gradbgo.getorbit("y", "z", ids=1)
        # x_go2, y_go2     = run_gradbgo.getorbit("y", "z", ids=2)

        # xi_gc            = run_gradbgc.getstate("y", state="ini")
        # xf_gc, deltat_gc = run_gradbgc.getstate("y", "mileage", state="end")
        # x_gc1, y_gc1     = run_gradbgc.getorbit("y", "z", ids=1)
        # x_gc2, y_gc2     = run_gradbgc.getorbit("y", "z", ids=2)

        # # Analytical values
        # gamma   = physlib.gamma_energy(m, ekin)
        # vnorm   = physlib.vnorm_gamma(gamma)
        # bvec    = run_gradbgo.bfield.read()["bxyz"].ravel()
        # gradb   = run_gradbgo.bfield.read()["jacobian"][0,:]
        # bnorm   = np.sqrt(np.sum(bvec**2)) * unyt.T
        # v_gradb = gamma * m * ( (1.0 - pitch**2) * vnorm**2 / ( q * bnorm ) ) \
        #     * 0.5 * np.cross(bvec, gradb) / bnorm**2 * unyt.T**2 / unyt.m
        # v_gradb.convert_to_units("m/s")
        # time    = run_gradbgo.options.read()["ENDCOND_LIM_SIMTIME"] * unyt.s

        # # Numerical values
        # vgo_gradb = (xf_go - xi_go) / deltat_go
        # vgc_gradb = (xf_gc - xi_gc) / deltat_gc

        # Plot
        #h3a.plot(x_go1, y_go1)
        #h3b.plot(x_go2, y_go2)
        #h3a.plot(x_gc1, y_gc1)
        #h3b.plot(x_gc2, y_gc2)


@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_drifts_exb(ascot, method, inspect, plot):

    if not inspect:
        parameters = SimulationOptions(
            simulation_mode=1,
            enable_orbit_following=True,
            collect_orbit=True,
            activate_simulation_time_limits=True,
            max_mileage=1e-7,
            )

        if method == "go":
            parameters.simulation.simulation_mode = 1
            parameters.simulation.timestep=1e-10
            parameters.orbit.number_of_points_per_marker=10002
            parameters.orbit.interval=1e-10
        else:
            parameters.simulation.simulation_mode = 2
            parameters.simulation.timestep = 1e-9
            parameters.orbit.number_of_points_per_marker=1002
            parameters.orbit.interval=1e-9
            if method == "gc-fixedstep":
                parameters.simulation.enable_adaptive = False

        ascot.data.create_guidingcentermarker(
            species="electron",
            ids=np.array([1, 2]),
            r=5*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=100e6*unyt.eV,
            pitch=0.5,
            gyroangle=0.0*unyt.rad,
            charge=np.array([1, -1])*unyt.e,
            )

        ascot.data.create_bfieldcartesian(
            bxyz=np.array([5., 0., 0.])*unyt.T,
            jacobian=np.full((3, 3), 0.)*unyt.T/unyt.m,
            axisrz=np.array([1., 0.])*unyt.m,
            rhoval=1.,
            )

        ascot.data.create_efieldcartesian(
            exyz=np.array([0, 1e6, 0])*unyt.V/unyt.m,
            )

        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    x1, y1, z1 = run.getorbit("x", "y", "z", mode="particle", filter=[1])
    x2, y2, z2 = run.getorbit("x", "y", "z", mode="particle", filter=[2])
    import matplotlib.pyplot as plt
    plt.plot(y1,z1)
    plt.plot(y2,z2)
    plt.show()

def test_drifts_gradb(ascot, method, inspect, plot):
    if not inspect:
        parameters = SimulationOptions(
            simulation_mode=1,
            enable_adaptive=False,
            use_explicit_fixedstep=True,
            explicit_fixedstep=1e-10,
            enable_orbit_following=True,
            collect_orbit=True,
            number_of_points_per_marker=10002,
            interval=1e-11,
            activate_simulation_time_limits=True,
            max_mileage=1e-7,
            activate_real_time_limit=True,
            max_real_time=5.0,
        )

        ascot.data.create_guidingcentermarker(
            species="electron",
            ids=np.array([1, 2]),
            r=5*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=100e6*unyt.eV,
            pitch=0.5,
            gyroangle=0.0*unyt.rad,
            charge=np.array([1, -1])*unyt.e,
            )

        jacobian = np.array(
            [
                [0., 0., 0.1],
                [0., 0., 0.],
                [0., 0., 0.],
            ]
            )*unyt.T/unyt.m
        ascot.data.create_bfieldcartesian(
            bxyz=np.array([5., 0., 0.])*unyt.T,
            jacobian=jacobian,
            axisrz=np.array([1., 0.])*unyt.m,
            rhoval=1.,
            )

        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    #print(run.getorbit("x", "y", "z", mode="guidingcenter", filter=[1]))
    x1, y1, z1 = run.getorbit("x", "y", "z", mode="particle", filter=[1])
    x2, y2, z2 = run.getorbit("x", "y", "z", mode="particle", filter=[2])
    #print(z1.size, z2.size)
    import matplotlib.pyplot as plt
    plt.plot(y1,z1)
    plt.plot(y2,z2)
    plt.show()


    # def init_elementary(self):
    #     """Initialize data for the elementary test.

    #     This test covers elementary results such as the gyromotion and drifts by
    #     verifying that ASCOT5 reproduces:

    #     1. correct gyroradius and gyrofrequency in an uniform magnetic field.
    #     2. correct E X B drift in uniform electromagnetic field.
    #     3. correct gradB drift in a magnetic field with constant gradient.

    #     Test 1 is done with GO mode and 2 and 3 with both GO and GC modes,
    #     latter using the fixed step scheme. Tests are done in Cartesian
    #     electromagnetic field (B_TC and E_TC) and without collisions.
    #     The test particles are an energetic electron and a positron so
    #     the tests also verify that ASCOT5 is valid in the relativistic regime.
    #     """
    #     if hasattr(self.ascot.data.options, PhysTest.tag_elementary_gyro):
    #         warnings.warn("Inputs already present: Test elementary")
    #         return
    #     init = self.ascot.data.create_input

    #     # Options
    #     opt = Opt.get_default()
    #     opt.update({
    #         "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
    #         "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_SIMTIMELIM" : 1,
    #         "ENDCOND_LIM_SIMTIME" : 2e-9, "ENABLE_ORBIT_FOLLOWING" : 1,
    #         "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
    #         "ORBITWRITE_INTERVAL" : 1e-11, "ORBITWRITE_NPOINT" : 202
    #     })
    #     init("opt", **opt, desc=PhysTest.tag_elementary_gyro)
    #     opt.update({
    #         "FIXEDSTEP_USERDEFINED" : 1e-10, "ENDCOND_LIM_SIMTIME" : 1e-7,
    #         "ORBITWRITE_NPOINT" : 10002
    #     })
    #     init("opt", **opt, desc=PhysTest.tag_elementary_exbgo)
    #     init("opt", **opt, desc=PhysTest.tag_elementary_gradbgo)
    #     opt.update({
    #         "SIM_MODE" : 2, "FIXEDSTEP_USERDEFINED" : 1e-9,
    #         "ORBITWRITE_INTERVAL" : 1e-9, "ORBITWRITE_NPOINT" : 102
    #     })
    #     init("opt", **opt, desc=PhysTest.tag_elementary_exbgc)
    #     init("opt", **opt, desc=PhysTest.tag_elementary_gradbgc)

    #     # Magnetic field
    #     d = {"bxyz" : np.array([5, 0, 0]),
    #          "jacobian" : np.array([0,0,0,0,0,0,0,0,0]),
    #          "rhoval" : 1.5}
    #     init("B_TC", **d, desc=PhysTest.tag_elementary_gyro)
    #     init("B_TC", **d, desc=PhysTest.tag_elementary_exbgo)
    #     init("B_TC", **d, desc=PhysTest.tag_elementary_exbgc)

    #     d.update({"jacobian" : np.array([0,0,0.1,0,0,0,0,0,0])})
    #     init("B_TC", **d, desc=PhysTest.tag_elementary_gradbgo)
    #     init("B_TC", **d, desc=PhysTest.tag_elementary_gradbgc)

    #     # Electric field
    #     d = {"exyz" : np.array([0, 0, 0])}
    #     init("E_TC", **d, desc=PhysTest.tag_elementary_gyro)
    #     init("E_TC", **d, desc=PhysTest.tag_elementary_gradbgo)
    #     init("E_TC", **d, desc=PhysTest.tag_elementary_gradbgc)

    #     d.update({"exyz" : np.array([0, 1e6, 0])})
    #     init("E_TC", **d, desc=PhysTest.tag_elementary_exbgo)
    #     init("E_TC", **d, desc=PhysTest.tag_elementary_exbgc)

    #     # Marker input is an electron and positron
    #     mrk = Marker.generate("gc", n=2, species="electron")
    #     mrk["charge"]    = np.array([1, -1])
    #     mrk["r"][:]      = 5
    #     mrk["phi"][:]    = 90
    #     mrk["z"][:]      = 0
    #     mrk["zeta"][:]   = 0
    #     mrk["energy"][:] = 100e6
    #     mrk["pitch"][:]  = 0.5

    #     for tag in [PhysTest.tag_elementary_gyro,
    #                 PhysTest.tag_elementary_exbgo,
    #                 PhysTest.tag_elementary_exbgc,
    #                 PhysTest.tag_elementary_gradbgo,
    #                 PhysTest.tag_elementary_gradbgc]:
    #         init("gc", **mrk, desc=tag)

    # def run_elementary(self):
    #     """Run elementary test.
    #     """
    #     if hasattr(self.ascot.data, PhysTest.tag_elementary_gyro):
    #         warnings.warn("Results already present: Test elementary")
    #         return
    #     for tag in [PhysTest.tag_elementary_gyro,
    #                 PhysTest.tag_elementary_exbgo,
    #                 PhysTest.tag_elementary_exbgc,
    #                 PhysTest.tag_elementary_gradbgo,
    #                 PhysTest.tag_elementary_gradbgc]:
    #         self._activateinputs(tag)
    #         self._runascot(tag)

    # def check_elementary(self):
    #     """Verify and plot elementary test results.

    #     Returns
    #     -------
    #     passed : bool
    #         True if the test passed.
    #     """
    #     run_gyro    = self.ascot.data[PhysTest.tag_elementary_gyro]
    #     run_exbgo   = self.ascot.data[PhysTest.tag_elementary_exbgo]
    #     run_exbgc   = self.ascot.data[PhysTest.tag_elementary_exbgc]
    #     run_gradbgo = self.ascot.data[PhysTest.tag_elementary_gradbgo]
    #     run_gradbgc = self.ascot.data[PhysTest.tag_elementary_gradbgc]

    #     # Initialize plots
    #     fig = a5plt.figuredoublecolumn()
    #     gs = GridSpec(2, 3, figure=fig)
    #     h1a = fig.add_subplot(gs[0,0])
    #     h2a = fig.add_subplot(gs[0,1])
    #     h2b = fig.add_subplot(gs[1,1])
    #     h3a = fig.add_subplot(gs[0,2])
    #     h3b = fig.add_subplot(gs[1,2])

    #     h1a.set_title("Gyro-motion")
    #     h2a.set_title("ExB drift")
    #     h3a.set_title("Grad-B drift")

    #     h1a.set_ylabel("z [m]")
    #     h1a.set_xlabel("y [m]")

    #     h1a.set_aspect("equal", adjustable="box")
    #     h2a.set_aspect("equal", adjustable="box")
    #     h2b.set_aspect("equal", adjustable="box")
    #     h3a.set_aspect("equal", adjustable="box")
    #     h3b.set_aspect("equal", adjustable="box")

    #     def plotarrow(ax, ini, end, txt):
    #         """Plot line with an arrow
    #         """
    #         ax.annotate("", xytext=ini, xy=end,
    #                     arrowprops={"arrowstyle":"->", "color":"black"})
    #         ax.plot([ini[0], end[0]], [ini[1], end[1]], color="black")
    #         ax.annotate(txt, ini + .005)

    #     ## Gyro test ##
    #     m, q, pitch, ekin, x0, y0 = \
    #         run_gyro.getstate("mass", "charge", "pitch", "ekin", "y", "z")
    #     x, y, time = run_gyro.getorbit("y", "z", "time", ids=1)

    #     # Analytical values
    #     bnorm = np.sqrt(np.sum(run_gyro.bfield.read()["bxyz"]**2)) * unyt.T
    #     larmorrad_ana = physlib.gyrolength(m, q, ekin, pitch, bnorm).to("m")
    #     gyrofreq_ana  = physlib.gyrofrequency(m, q, ekin, bnorm).to("rad/s")

    #     # Numerical values
    #     larmorrad_go = np.mean( np.sqrt( (x - x0[0])**2 + (y - y0[0])**2 ) )
    #     gyrofreq_go  = np.sum( np.sqrt( np.diff(x)**2 + np.diff(y)**2 ) ) \
    #         * unyt.rad / (larmorrad_go * time[-1])

    #     # Plot
    #     orbx = x - x0[0]
    #     orby = y - y0[0]

    #     anax   = larmorrad_ana[0] * np.sin(np.linspace(0, 2*np.pi, 360))
    #     anay   = larmorrad_ana[0] * np.cos(np.linspace(0, 2*np.pi, 360))
    #     fangle = np.arctan2(orby[0], orbx[0]) * unyt.rad \
    #         + gyrofreq_ana[0] * time[-1]
    #     fx = larmorrad_ana[0] * np.cos(-fangle) # Minus sign because this we are
    #     fy = larmorrad_ana[0] * np.sin(-fangle) # in LHS coordinate system.

    #     h1a.plot(orbx, orby, linewidth=3)
    #     h1a.plot(anax, anay, linestyle="--", color="black", alpha=0.7)
    #     h1a.scatter(orbx[-1], orby[-1], marker="o", color="C0", zorder=3)
    #     h1a.scatter(fx, fy, marker="x", color="black", zorder=4)

    #     ## E x B drift
    #     yi_go            = run_exbgo.getstate("z", state="ini")
    #     yf_go, deltat_go = run_exbgo.getstate("z", "mileage", state="end")
    #     x_go1, y_go1     = run_exbgo.getorbit("y", "z", ids=1)
    #     x_go2, y_go2     = run_exbgo.getorbit("y", "z", ids=2)

    #     yi_gc            = run_exbgc.getstate("z", state="ini")
    #     yf_gc, deltat_gc = run_exbgc.getstate("z", "mileage", state="end")
    #     x_gc1, y_gc1     = run_exbgc.getorbit("y", "z", ids=1)
    #     x_gc2, y_gc2     = run_exbgc.getorbit("y", "z", ids=2)

    #     # Analytical values
    #     bvec  = run_exbgo.bfield.read()["bxyz"].ravel()
    #     evec  = run_exbgo.efield.read()["exyz"].ravel()
    #     v_ExB = np.cross(evec, bvec) / np.inner(bvec, bvec) * unyt.m / unyt.s
    #     time  = run_exbgo.options.read()["ENDCOND_LIM_SIMTIME"] * unyt.s

    #     # Numerical values
    #     vgo_ExB = (yf_go - yi_go) / deltat_go
    #     vgc_ExB = (yf_gc - yi_gc) / deltat_gc

    #     # Plot
    #     h2a.plot(x_go1, y_go1)
    #     h2b.plot(x_go2, y_go2)
    #     h2a.plot(x_gc1, y_gc1)
    #     h2b.plot(x_gc2, y_gc2)

    #     ini = np.array([x_gc1[0].v + .01, y_gc1[0].v])
    #     end = np.array([x_gc1[0].v + .01, (y_gc1[0] + v_ExB[2]*deltat_gc[1]).v])
    #     plotarrow(h2a, ini, end, r"$e^+$")
    #     ini = np.array([x_gc2[0].v + .01, y_gc2[0].v])
    #     end = np.array([x_gc2[0].v + .01, (y_gc2[0] + v_ExB[2]*deltat_gc[0]).v])
    #     plotarrow(h2b, ini, end, r"$e^-$")

    #     ## gradB drift
    #     pitch, ekin, q, m = run_gradbgo.getstate(
    #         "pitch", "ekin", "charge", "mass", ids=1)
    #     xi_go            = run_gradbgo.getstate("y", state="ini")
    #     xf_go, deltat_go = run_gradbgo.getstate("y", "mileage", state="end")
    #     x_go1, y_go1     = run_gradbgo.getorbit("y", "z", ids=1)
    #     x_go2, y_go2     = run_gradbgo.getorbit("y", "z", ids=2)

    #     xi_gc            = run_gradbgc.getstate("y", state="ini")
    #     xf_gc, deltat_gc = run_gradbgc.getstate("y", "mileage", state="end")
    #     x_gc1, y_gc1     = run_gradbgc.getorbit("y", "z", ids=1)
    #     x_gc2, y_gc2     = run_gradbgc.getorbit("y", "z", ids=2)

    #     # Analytical values
    #     gamma   = physlib.gamma_energy(m, ekin)
    #     vnorm   = physlib.vnorm_gamma(gamma)
    #     bvec    = run_gradbgo.bfield.read()["bxyz"].ravel()
    #     gradb   = run_gradbgo.bfield.read()["jacobian"][0,:]
    #     bnorm   = np.sqrt(np.sum(bvec**2)) * unyt.T
    #     v_gradb = gamma * m * ( (1.0 - pitch**2) * vnorm**2 / ( q * bnorm ) ) \
    #         * 0.5 * np.cross(bvec, gradb) / bnorm**2 * unyt.T**2 / unyt.m
    #     v_gradb.convert_to_units("m/s")
    #     time    = run_gradbgo.options.read()["ENDCOND_LIM_SIMTIME"] * unyt.s

    #     # Numerical values
    #     vgo_gradb = (xf_go - xi_go) / deltat_go
    #     vgc_gradb = (xf_gc - xi_gc) / deltat_gc

    #     # Plot
    #     h3a.plot(x_go1, y_go1)
    #     h3b.plot(x_go2, y_go2)
    #     h3a.plot(x_gc1, y_gc1)
    #     h3b.plot(x_gc2, y_gc2)

    #     ini = np.array([x_gc1[0].v, y_gc1[0].v + .01])
    #     end = np.array([x_gc1[0] + v_gradb[1]*deltat_gc[0], y_gc1[0].v + .01])
    #     plotarrow(h3a, ini, end, r"$e^+$")
    #     ini = np.array([x_gc2[0].v, y_gc2[0] + .01])
    #     end = np.array([x_gc2[0] - v_gradb[1]*deltat_gc[0], y_gc2[0].v + .01])
    #     plotarrow(h3b, ini, end, r"$e^-$")

    #     print("Test elementary:")
    #     passed = True

    #     print("  Gyroradius %e, Gyrofrequency %e (gyro-orbit)" \
    #           % (larmorrad_go, gyrofreq_go))
    #     err = np.abs(larmorrad_go - larmorrad_ana[0]) / larmorrad_ana[0]
    #     if err > 5e-8: print("Error: %e (FAILED)" % err); passed = False

    #     print("  Gyroradius %e, Gyrofrequency %e (expected)" \
    #           % (larmorrad_ana[0], gyrofreq_ana[0]))
    #     err = np.abs(gyrofreq_go - gyrofreq_ana[0]) / gyrofreq_ana[0]
    #     if err > 5e-4: print("Error: %e (FAILED)" % err); passed = False

    #     print("  ExB drift for positrons %e electrons %e (gyro-orbit)" \
    #           % (vgo_ExB[0], vgo_ExB[1]) )
    #     print("  ExB drift for positrons %e electrons %e (guiding center)" \
    #           % (vgc_ExB[0], vgc_ExB[1]) )
    #     print("  ExB drift for positrons %e electrons %e (expected)" \
    #           % (v_ExB[2], v_ExB[2]) )
    #     err = np.amax([np.abs( ( vgo_ExB[0] - v_ExB[2] ) / v_ExB[2] ),
    #                    np.abs( ( vgo_ExB[1] - v_ExB[2] ) / v_ExB[2] ),
    #                    np.abs( ( vgc_ExB[0] - v_ExB[2] ) / v_ExB[2] ),
    #                    np.abs( ( vgc_ExB[1] - v_ExB[2] ) / v_ExB[2] )])
    #     if err > 5e-5: print("Error: %e (FAILED)" % err); passed = False

    #     print("  Grad-B drift for positrons %e electrons %e (gyro-orbit)" \
    #           % (vgo_gradb[0], vgo_gradb[1]) )
    #     print("  Grad-B drift for positrons %e electrons %e (guiding center)" \
    #           % (vgc_gradb[0], vgc_gradb[1]) )
    #     print("  Grad-B drift for positrons %e electrons %e (expected)" \
    #           % (v_gradb[1], -v_gradb[1]) )
    #     err = np.amax([np.abs( ( vgo_gradb[0] - v_gradb[1] ) / v_gradb[1] ),
    #                    np.abs( ( vgo_gradb[1] + v_gradb[1] ) / v_gradb[1] ),
    #                    np.abs( ( vgc_gradb[0] - v_gradb[1] ) / v_gradb[1] ),
    #                    np.abs( ( vgc_gradb[1] + v_gradb[1] ) / v_gradb[1] )])
    #     if err > 1e-3: print("Error: %e (FAILED)" % err); passed = False

    #     return passed
