"""Verify the guiding center transformation.

The transformation is verified by comparing the guiding-center orbit to the
"actual" guiding center orbit obtained by performing guiding center
transformation at each point along the particle orbit. The test also shows that
including the first order terms in the momentum space (in the transformation)
imporves accuracy.
"""
import numpy as np
import unyt

import pytest

from a5py import Ascot, SimulationOptions, physlib
from a5py import plotting as a5plt
from a5py.templates import PremadeMagneticField

def test_gctransform_compare_actual(ascot, inspect, plot):
    if not inspect:
        parameters = SimulationOptions.from_dict(
            physics={
                "enable_orbit_following": True,
            },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 1e-7,
            },
            orbit={
                "collect": "interval",
            },
        )
        parameters = SimulationOptions(
                simulation_mode=2,
                enable_adaptive=False,
                timestep=1e-10,
                enable_orbit_following=True,
                collect_orbit=True,
                number_of_points_per_marker=75002,
                interval=4e-10,
                activate_simulation_time_limits=True,
                max_mileage=3e-5,
            )

        ascot.data.create_guidingcentermarker(
            species="alpha",
            r=7.6*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=3.5e6*unyt.eV,
            pitch=0.4,
            gyroangle=2.0*unyt.rad,
            charge=2*unyt.e,
            )
        PremadeMagneticField(ascot, field="iter-circular").create_input()



def test_gctransform_verify_firstorder(ascot, inspect, plot):
    if not inspect:
        PremadeMagneticField(ascot, field="iter-circular").create_input()

        ascot.data.create_guidingcentermarker(
            species="alpha",
            r=7.6*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=3.5e6*unyt.eV,
            pitch=0.4,
            gyroangle=2.0*unyt.rad,
            charge=2*unyt.e,
            )

        nmarker = 10
        parameters = SimulationOptions(
            simulation_mode=1,
            timestep=1e-11,
            enable_orbit_following=True,
            collect_orbit=True,
            number_of_points_per_marker=nmarker,
            interval=2e-10,
            activate_simulation_time_limits=True,
            max_mileage=2e-9,
            )

        ascot.data.bfield.active.stage()
        ascot.data.marker.active.stage()
        run = ascot.simulate(params=parameters)
        r, z, phi, vr, vphi, vz = run.getorbit(
            "r", "z", "phi", "vr", "vphi", "vz",
            )

        ascot.data.create_guidingcentermarker(
            species="alpha",
            r=r,
            z=z,
            phi=phi,
            vr=vr,
            vz=vz,
            vphi=vphi,
            charge=2*unyt.e,
            )

        parameters = SimulationOptions(
            simulation_mode=2,
            timestep=1e-9,
            enable_adaptive=False,
            enable_orbit_following=True,
            collect_orbit=True,
            number_of_points_per_marker=10**4,
            interval=3e-9,
            activate_simulation_time_limits=True,
            max_mileage=3e-5,
            )

        ascot.data.bfield.active.stage()
        ascot.data.marker.active.stage()
        run = ascot.simulate(params=parameters)

        parameters.physics.disable_first_order_gctransformation = 1
        ascot.data.bfield.active.stage()
        ascot.data.marker.active.stage()
        run = ascot.simulate(params=parameters)
#     def check_gctransform(self):
#         """Check test.
#         """
#         run_go     = self.ascot.data[PhysTest.tag_gctransform_go]
#         run_gc     = self.ascot.data[PhysTest.tag_gctransform_gc]
#         run_go2gc  = self.ascot.data[PhysTest.tag_gctransform_go2gc]
#         run_zeroth = self.ascot.data[PhysTest.tag_gctransform_zeroth]
#         run_first  = self.ascot.data[PhysTest.tag_gctransform_first]

#         # Initialize plots
#         fig = a5plt.figuredoublecolumn(3/2)
#         gs  = GridSpec(3, 3, figure=fig)
#         h1a = fig.add_subplot(gs[0,0])
#         h1b = fig.add_subplot(gs[1,0])
#         h1c = fig.add_subplot(gs[2,0])
#         h2  = fig.add_subplot(gs[:,1])
#         h3  = fig.add_subplot(gs[:,2])

#         h1a.set_xlim(0, 30)
#         h1b.set_xlim(0, 30)
#         h1c.set_xlim(0, 30)
#         h1a.set_xticklabels([])
#         h1b.set_xticklabels([])

#         h2.set_xlim(6.2, 6.7)
#         h2.set_ylim(1.0, 1.6)
#         h2.set_aspect("equal", adjustable="box")
#         h3.set_xlim(6.2, 6.7)
#         h3.set_ylim(1.0, 1.6)
#         h3.set_aspect("equal", adjustable="box")
#         h3.set_yticklabels([])

#         h1c.set_xlabel("Time [µs]")
#         h1a.set_ylabel(r"$\mu$")
#         h1b.set_ylabel(r"$E_\mathrm{kin}$")
#         h1c.set_ylabel(r"$P_\mathrm{ctor}$")
#         h2.set_xlabel("R [m]")
#         h3.set_xlabel("R [m]")
#         h3.set_ylabel("z [m]")
#         h3.yaxis.set_label_position("right")
#         h3.yaxis.tick_right()

#         # Get data and plot
#         self.ascot.input_init(run=run_go.get_qid(), bfield=True)
#         rgo, zgo, tgo, ego, mugo, pgo = run_go.getorbit(
#             "r", "z", "mileage", "ekin", "mu", "ptor")
#         rgc, zgc, tgc, egc, mugc, pgc = run_gc.getorbit(
#             "r", "z", "mileage", "ekin", "mu", "ptor")
#         rgo2gc, zgo2gc, tgo2gc, ego2gc, mugo2gc, pgo2gc = run_go2gc.getorbit(
#             "r", "z", "mileage", "ekin", "mu", "ptor")
#         self.ascot.input_free()

#         h1a.plot(tgo.to("µs"), mugo)
#         h1a.plot(tgo2gc.to("µs"), mugo2gc)
#         h1a.plot(tgc.to("µs"), mugc)

#         h1b.plot(tgo.to("µs"), ego)
#         h1b.plot(tgo2gc.to("µs"), ego2gc)
#         h1b.plot(tgc.to("µs"), egc)

#         h1c.plot(tgo.to("µs"), pgo)
#         h1c.plot(tgo2gc.to("µs"), pgo2gc)
#         h1c.plot(tgc.to("µs"), pgc)

#         h2.plot(rgo, zgo)
#         h2.plot(rgo2gc, zgo2gc)
#         h2.plot(rgc, zgc)

#         nrep = run_zeroth.getstate("ids").size
#         for i in range(nrep):
#             r, z = run_zeroth.getorbit("r", "z", ids=i)
#             h3.plot(r, z, color="C1")
#         for i in range(nrep):
#             r, z = run_first.getorbit("r", "z", ids=i)
#             h3.plot(r.v+0.01, z, color="C2")
#         h3.plot(rgo2gc, zgo2gc, color="black")

#         # Verify results by interpolating the orbits at fixed intervals and then
#         # calculating sum-of-squares of the difference between go and gc and
#         # go2gc and gc. The latter should be smaller if the guiding center
#         # transformation works.
#         t = np.linspace(0, 1e-6, 1000) * unyt.s
#         mugo    = np.interp(t, tgo,    mugo)
#         mugc    = np.interp(t, tgc,    mugc)
#         mugo2gc = np.interp(t, tgo2gc[1:], mugo2gc[1:])
#         ego     = np.interp(t, tgo,    ego)
#         egc     = np.interp(t, tgc,    egc)
#         ego2gc  = np.interp(t, tgo2gc[1:], ego2gc[1:])
#         pgo     = np.interp(t, tgo,    pgo)
#         pgc     = np.interp(t, tgc,    pgc)
#         pgo2gc  = np.interp(t, tgo2gc[1:], pgo2gc[1:])

#         print("Test GC transformation:")
#         passed = True

#         mugo = np.mean(mugo); mugc = np.mean(mugc); mugo2gc = np.mean(mugo2gc)
#         ego  = np.mean(ego);  egc  = np.mean(egc);  ego2gc  = np.mean(ego2gc)
#         pgo  = np.mean(pgo);  pgc  = np.mean(pgc);  pgo2gc  = np.mean(pgo2gc)

#         print("  Mean value  GO        GC        GO2GC")
#         err = np.abs(mugc/mugo2gc-1)
#         print("  mu          %1.3e %1.3e %1.3e" % (mugo, mugc, mugo2gc))
#         if err > 1e-4: print("Error %e (FAILED)" % err); passed = False
#         err = np.abs(egc/ego2gc-1)
#         print("  Energy      %1.3e %1.3e %1.3e" % (ego, egc, ego2gc))
#         if err > 2e-4: print("Error %e (FAILED)" % err); passed = False
#         err = np.abs(pgc/pgo2gc-1)
#         print("  Pctor       %1.3e %1.3e %1.3e" % (-pgo, -pgc, -pgo2gc))
#         if err > 1e-4: print("Error %e (FAILED)" % err); passed = False

#         return passed