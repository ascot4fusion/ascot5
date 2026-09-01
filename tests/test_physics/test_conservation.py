"""Verify that the constants of motion in orbit-following are conserved."""
import numpy as np
import unyt

import pytest
from matplotlib.gridspec import GridSpec

from a5py import Ascot, SimulationOptions, physlib
from a5py import plotting as a5plt
from a5py.templates import PremadeMagneticField

@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_conservation(ascot, method, inspect, plot):
    if not inspect:
        parameters = SimulationOptions.from_dict(
            physics={
                "enable_orbit_following": True,
            },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 5e-6,
            },
            orbit={
                "collect": "interval",
            },
        )
        if method == "go":
            parameters.simulation.mode = "gyro-orbit"
            parameters.simulation.timestep = 1e-11
            parameters.orbit.buffer_size = 50002
            parameters.orbit.interval = 1e-10
        else:
            parameters.simulation.mode = "guiding-center"
            parameters.simulation.timestep = 1e-12
            parameters.orbit.buffer_size = 502
            parameters.orbit.interval = 1e-8
            if method == "gc-fixedstep":
                parameters.simulation.enable_adaptive = False
                parameters.simulation.adaptive_tolerance_orbit = 1e-14

        ascot.data.create_guidingcentermarker(
            species="electron",
            r=7.6*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=10e6*unyt.eV,
            pitch=np.array([0.4, 0.9]),
            gyroangle=0.0*unyt.rad,
            charge=np.array([1., -1.])*unyt.e,
            )

        PremadeMagneticField(ascot, field="iter-circular").create_input()
        ascot.data.create_efieldcartesian(
            exyz=np.array([0, 0, 0])*unyt.V/unyt.m,
            )

        run = ascot.simulate(params=parameters)
    else:
        run = ascot.data.active

    t1, e1, mu1, p1, r1, z1 = run.getorbit(
        "mileage", "ekin", "mu", "ptor", "r", "z", filter=[1])
    t2, e2, mu2, p2, r2, z2 = run.getorbit(
        "mileage", "ekin", "mu", "ptor", "r", "z", filter=[2])

    if plot:
        fig = a5plt.figure()
        gs = GridSpec(3, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0])
        ax3 = fig.add_subplot(gs[2,0])
        ax4 = fig.add_subplot(gs[:,1])

        def plot_difference(t, e, mu, p, color):
            ax1.plot(t, e - e[0], color=color)
            ax2.plot(t, mu - mu[0], color=color)
            ax3.plot(t, p - p[0], color=color)

        plot_difference(t1, e1, mu1, p1, "C0")
        plot_difference(t2, e2, mu2, p2, "C1")

        ax4.plot(r1, z1, color="C0")
        ax4.plot(r2, z2, color="C1")

        ax3.set_xlabel("Mileage [m]")
        ax4.set_xlabel("$r$ [m]")
        ax4.set_ylabel("$z$ [m]")

        a5plt.show()

#     def check_orbitfollowing(self):
#         """Check test.
#         """
#         run_go  = self.ascot.data[PhysTest.tag_orbfol_go]
#         run_gcf = self.ascot.data[PhysTest.tag_orbfol_gcf]
#         run_gca = self.ascot.data[PhysTest.tag_orbfol_gca]

#         # Initialize plots
#         fig = a5plt.figuredoublecolumn(3/2)
#         gs = GridSpec(3, 4, figure=fig)
#         h1a = fig.add_subplot(gs[0,0])
#         h2a = fig.add_subplot(gs[1,0])
#         h3a = fig.add_subplot(gs[2,0])
#         h4a = fig.add_subplot(gs[:,1])
#         h1b = fig.add_subplot(gs[0,2])
#         h2b = fig.add_subplot(gs[1,2])
#         h3b = fig.add_subplot(gs[2,2])
#         h4b = fig.add_subplot(gs[:,3])

#         h4a.set_aspect("equal", adjustable="box")
#         h4b.set_aspect("equal", adjustable="box")

#         h1a.set_xlim(0, 5)
#         h2a.set_xlim(0, 5)
#         h3a.set_xlim(0, 5)
#         h1b.set_xlim(0, 5)
#         h2b.set_xlim(0, 5)
#         h3b.set_xlim(0, 5)

#         h4a.set_xlim(5, 8)
#         h4a.set_ylim(-2, 2)
#         h4b.set_xlim(5, 8)
#         h4b.set_ylim(-2, 2)

#         h1a.set_xticklabels([])
#         h2a.set_xticklabels([])
#         h1b.set_xticklabels([])
#         h2b.set_xticklabels([])

#         h4a.set_yticklabels([])
#         h4b.set_xlabel("R [m]")
#         h4b.set_ylabel("z [m]")
#         h4b.yaxis.set_label_position("right")
#         h4b.yaxis.tick_right()

#         h1a.set_ylabel(r"$(E-E_0)/E_0$")
#         h2a.set_ylabel(r"$(\mu-\mu_0)/\mu_0$")
#         h3a.set_ylabel(r"$(P-P_0)/P_0$")

#         h3a.set_xlabel("Time [µs]")
#         h3b.set_xlabel("Time [µs]")

#         def plotreldiff(ax1, ax2, ax3, t, q1, q2, q3, **kwargs):
#             """Plot relative diffenrence of given quantities on given axes.
#             """
#             ax1.plot(t, ( q1 - q1[0] ) / q1[0], **kwargs)
#             ax2.plot(t, ( q2 - q2[0] ) / q2[0], **kwargs)
#             ax3.plot(t, ( q3 - q3[0] ) / q3[0], **kwargs)

#         def fails(t, q, eps, qnt, otype, mode):
#             """Check if the change in time of a given quantity is below given
#             tolerance
#             """
#             err = _linfit(t, (q - q[0]) / q[0])
#             msg = "Rate of change in %6s (%s/%11s): %e Tolerance: %e" \
#                 % (qnt, otype, mode, err, eps)
#             if np.abs(err) > eps:
#                 msg += " (FAILED)"
#                 print(msg)
#                 return True
#             print(msg)
#             return False

#         # Numerical values
#         self.ascot.input_init(run=run_go.get_qid(), bfield=True)
#         tgo1, ego1, mugo1, pgo1, rgo1, zgo1 = run_go.getorbit(
#             "mileage", "ekin", "mu", "ptor", "r", "z", ids=1)
#         tgo2, ego2, mugo2, pgo2, rgo2, zgo2 = run_go.getorbit(
#             "mileage", "ekin", "mu", "ptor", "r", "z", ids=2)
#         tgcf1, egcf1, mugcf1, pgcf1, rgcf1, zgcf1 = run_gcf.getorbit(
#             "mileage", "ekin", "mu", "ptor", "r", "z", ids=1)
#         tgcf2, egcf2, mugcf2, pgcf2, rgcf2, zgcf2 = run_gcf.getorbit(
#             "mileage", "ekin", "mu", "ptor", "r", "z", ids=2)
#         tgca1, egca1, mugca1, pgca1, rgca1, zgca1 = run_gca.getorbit(
#             "mileage", "ekin", "mu", "ptor", "r", "z", ids=1)
#         tgca2, egca2, mugca2, pgca2, rgca2, zgca2 = run_gca.getorbit(
#             "mileage", "ekin", "mu", "ptor", "r", "z", ids=2)
#         self.ascot.input_free()

#         # Plot
#         plotreldiff(h1a, h2a, h3a, tgo1.to("µs"),  ego1,  mugo1,  pgo1,
#                     color="C0")
#         plotreldiff(h1b, h2b, h3b, tgo2.to("µs"),  ego2,  mugo2,  pgo2,
#                     color="C0")
#         plotreldiff(h1a, h2a, h3a, tgcf1.to("µs"), egcf1, mugcf1, pgcf1,
#                     color="C1")
#         plotreldiff(h1b, h2b, h3b, tgcf2.to("µs"), egcf2, mugcf2, pgcf2,
#                     color="C1")
#         plotreldiff(h1a, h2a, h3a, tgca1.to("µs"), egca1, mugca1, pgca1,
#                     color="C2")
#         plotreldiff(h1b, h2b, h3b, tgca2.to("µs"), egca2, mugca2, pgca2,
#                     color="C2")

#         h4a.plot(rgo1,  zgo1,  color="C0")
#         h4a.plot(rgcf1, zgcf1, color="C1")
#         h4a.plot(rgca1, zgca1, color="C2")
#         h4b.plot(rgo2,  zgo2,  color="C0")
#         h4b.plot(rgcf2, zgcf2, color="C1")
#         h4b.plot(rgca2, zgca2, color="C2")

#         print("Test orbit-following:")
#         passed = True

#         otype1 = "passing"; otype2 = "trapped"
#         simmode1 = "gyro-orbit"; simmode2 = "fixed GC"; simmode3 = "adaptive GC"
#         if fails(tgo1,  ego1,   5e-6, "energy", otype1, simmode1):passed = False
#         if fails(tgo2,  ego2,   5e-6, "energy", otype2, simmode1):passed = False
#         if fails(tgcf1, egcf1,  5e-6, "energy", otype1, simmode2):passed = False
#         if fails(tgcf2, egcf2,  5e-6, "energy", otype2, simmode2):passed = False
#         if fails(tgca1, egca1,  5e-6, "energy", otype1, simmode3):passed = False
#         if fails(tgca2, egca2,  5e-6, "energy", otype2, simmode3):passed = False
#         if fails(tgo1,  mugo1,   2e2, "mu",     otype1, simmode1):passed = False
#         if fails(tgo2,  mugo2,  2e-1, "mu",     otype2, simmode1):passed = False
#         if fails(tgcf1, mugcf1, 1e-9, "mu",     otype1, simmode2):passed = False
#         if fails(tgcf2, mugcf2, 1e-9, "mu",     otype2, simmode2):passed = False
#         if fails(tgca1, mugca1, 1e-9, "mu",     otype1, simmode3):passed = False
#         if fails(tgca2, mugca2, 1e-9, "mu",     otype2, simmode3):passed = False
#         if fails(tgo1,  pgo1,   5e-3, "Ptor",   otype1, simmode1):passed = False
#         if fails(tgo2,  pgo2,   5e-5, "Ptor",   otype2, simmode1):passed = False
#         if fails(tgcf1, pgcf1,  5e-2, "Ptor",   otype1, simmode2):passed = False
#         if fails(tgcf2, pgcf2,  5e-2, "Ptor",   otype2, simmode2):passed = False
#         if fails(tgca1, pgca1,  5e-2, "Ptor",   otype1, simmode3):passed = False
#         if fails(tgca2, pgca2,  5e-2, "Ptor",   otype2, simmode3):passed = False

#         return passed