"""Test the Coulomb collision operator."""
import numpy as np
import unyt

import pytest

from a5py import Ascot, SimulationOptions, physlib
from a5py import plotting as a5plt
from a5py.templates import PremadeMagneticField

@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_maxwellian(ascot: Ascot, method, inspect, plot):
    if not inspect:
        PremadeMagneticField(ascot, field="iter-circular").create_input()
        ascot.data.create_efieldcartesian(
            exyz=np.array([0, 0, 0])*unyt.V/unyt.m,
            )
        ascot.data.create_plasmalinear1d(
            species=["H"],
            rhogrid=np.array([0, 2]),
            ni=1e20*unyt.m**-3,
            Ti=1e3*unyt.eV,
            )
#         opt.update({
#             "ENABLE_DIST_5D" : 1,
#             "DIST_MIN_R"    : 4,  "DIST_MAX_R"    : 10, "DIST_NBIN_R"      : 1,
#             "DIST_MIN_PHI"  : 0,  "DIST_MAX_PHI"  : 360, "DIST_NBIN_PHI"   : 1,
#             "DIST_MIN_Z"    : -5, "DIST_MAX_Z"    : 5, "DIST_NBIN_Z"       : 1,
#             "DIST_MIN_TIME" : 0,  "DIST_MAX_TIME" : 2e-2, "DIST_NBIN_TIME" : 2,
#             "DIST_MIN_PPA"  : -2.5e-21, "DIST_MAX_PPA" : 2.5e-21,
#             "DIST_NBIN_PPA" : 140, "DIST_MIN_PPE" : 0, "DIST_MAX_PPE" : 2.5e-21,
#             "DIST_NBIN_PPE" : 80
#         })
        pol, rminor = 2*np.pi * np.random.rand(20,), 0.8
        ascot.data.create_guidingcentermarker(
            species="proton",
            r=(6.2 + rminor * np.cos(pol))*unyt.m,
            z=(0.0 + rminor * np.sin(pol))*unyt.m,
            ekin=1e3*unyt.eV,
            pitch=0.4,
            )

        parameters = SimulationOptions.from_dict(
            physics={
                "enable_orbit_following": True,
                "enable_coulomb_collisions": True,
                },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 2e-2,#2
            },
            histograms=[
                {
                    "dimensions": [
                        ("ekin", 0.*unyt.eV, 1e4*unyt.eV, 100),
                        ("pitch", -1, 1, 10),
                        ],
                    "charge_interval": (1, 1,),
                },
            ]
            )

        if method == "go":
            parameters.simulation.mode = "gyro-orbit"
            parameters.simulation.timestep=1e-8
        else:
            parameters.simulation.mode = "guiding-center"
            parameters.simulation.timestep = 2e-8
            parameters.simulation.enable_adaptive = False
            if method == "gc-adaptivestep":
                parameters.simulation.enable_adaptive = True
                parameters.simulation.adaptive_tolerance_orbit = 1e-6
                parameters.simulation.adaptive_tolerance_collisions = 1e-2

        run = ascot.simulate(params=parameters)

        # Analytical results
        plasma = ascot.data.plasma.active
        ne, Te = plasma.ne[0], plasma.Te[0]
        m_e = unyt.me
        m_p     = unyt.mp
        m_a     = 4.003*unyt.amu
        eps_0   = unyt.eps_0
        #alphaZ  = 2
        #clog    = 16
        #E0      = 3.5e6*unyt.eV
        #Emin    = 50 * Te

        # Thermal distribution with correct normalization
        dist0 = run.getdist(0, ("pitch", "ekin"))
        time = run.options.endconditions.max_mileage * unyt.s
        weight = np.sum(ascot.data.marker.active.weight)
        thermal = (
            2*np.sqrt(dist0.dimensions["ekin"]/np.pi) * np.power(Te, -3.0/2)
            * np.exp(-dist0.dimensions["ekin"]/Te) * time * weight
        ).to("1/eV")

        #ekin0 = run.getstate("ekin", state="ini")
        #ekin1 = run.getstate("ekin", state="end")
        #print(run.getdist(0, ("ekin", "pitch")))
        #dist = run.getdist(0, ("pitch", "ekin"))
        dist0.integrate({"pitch": None})
        dist1 = run.getdist(0, ("pitch", "ekin"))
        dist1.integrate({"ekin": None})
        _, ax = a5plt.subplots()
        ax.plot(dist0.dimensions["ekin"].to("eV"), dist0.distribution.to("1/eV"))
        ax.plot(dist0.dimensions["ekin"].to("eV"), thermal)
        _, ax = a5plt.subplots()
        ax.plot(dist1.dimensions["pitch"], dist1.histogram)
        a5plt.show()
        #plt.hist(ekin0)
        #plt.hist(ekin1)
        #plt.show()


@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_slowingdown(method, ascot: Ascot, inspect, plot):
    if not inspect:
        PremadeMagneticField(ascot, field="iter-circular").create_input()
        ascot.data.create_plasmalinear1d(
            species=["H"],
            rhogrid=np.array([0, 2]),
            ni=1e20*unyt.m**-3,
            Ti=1e3*unyt.eV,
            )

        nmrk = 200
        pol, rminor = 2*np.pi * np.random.rand(nmrk,), 0.8
        ascot.data.create_guidingcentermarker(
            species="alpha",
            r=(6.2 + rminor * np.cos(pol))*unyt.m,
            z=(0.0 + rminor * np.sin(pol))*unyt.m,
            ekin=3.5e3*unyt.eV,
            pitch=1.0 - 2.0 * np.random.rand(nmrk,),
            charge=2*unyt.e,
            )

        parameters = SimulationOptions()
        parameters.physics.enable_orbit_following = True
        parameters.physics.enable_coulomb_collisions = True
        parameters.endconditions.activate_energy_limits = True
        parameters.endconditions.min_energy = 1e3*unyt.eV
        parameters.endconditions.local_thermal_limit = 0.

        if method == "go":
            parameters.simulation.simulation_mode = 1
            parameters.simulation.timestep=2e-9
        else:
            parameters.simulation.simulation_mode = 2
            parameters.simulation.timestep = 3e-8
            parameters.simulation.enable_adaptive = False
            if method == "gc-adaptivestep":
                parameters.simulation.enable_adaptive = True
                parameters.simulation.adaptive_tolerance_orbit = 1e-6
                parameters.simulation.adaptive_tolerance_collisions = 1e-2

        ascot.data.bfield.active.stage()
        ascot.data.plasma.active.stage()
        ascot.data.marker.active.stage()
        run = ascot.simulate(params=parameters)

#         opt = Opt.get_default()
#         opt.update({
#             "ENABLE_DIST_5D" : 1,
#             "DIST_MIN_R" : 4, "DIST_MAX_R" : 10, "DIST_NBIN_R" : 1,
#             "DIST_MIN_PHI" : 0, "DIST_MAX_PHI" : 360, "DIST_NBIN_PHI" : 1,
#             "DIST_MIN_Z" : -5, "DIST_MAX_Z" : 5, "DIST_NBIN_Z" : 1,
#             "DIST_MIN_TIME" : 0, "DIST_MAX_TIME" : 1e0, "DIST_NBIN_TIME" : 1,
#             "DIST_MIN_PPA" : -1.3e-19, "DIST_MAX_PPA" : 1.3e-19,
#             "DIST_NBIN_PPA" : 200, "DIST_MIN_PPE" : 0, "DIST_MAX_PPE" : 1.3e-19,
#             "DIST_NBIN_PPE" : 100
#         })





#     def check_ccoll(self):
#         """Check Coulomb collision test.
#         """
#         run_tgo  = self.ascot.data[PhysTest.tag_ccoll_thermalgo]
#         run_tgcf = self.ascot.data[PhysTest.tag_ccoll_thermalgcf]
#         run_tgca = self.ascot.data[PhysTest.tag_ccoll_thermalgca]
#         run_sgo  = self.ascot.data[PhysTest.tag_ccoll_slowinggo]
#         run_sgcf = self.ascot.data[PhysTest.tag_ccoll_slowinggcf]
#         run_sgca = self.ascot.data[PhysTest.tag_ccoll_slowinggca]

#         fig = a5plt.figuredoublecolumn()
#         gs = GridSpec(2, 2, figure=fig)
#         h1 = fig.add_subplot(gs[0,0])
#         h2 = fig.add_subplot(gs[0,1])
#         h3 = fig.add_subplot(gs[1,0])
#         h4 = fig.add_subplot(gs[1,1])

#         # Analytical results
#         ne      = 1e20 / unyt.m**3
#         Te      = 1e3 * unyt.eV
#         m_e     = unyt.me
#         m_p     = unyt.mp
#         m_a     = 4.003*unyt.amu
#         eps_0   = unyt.eps_0
#         alphaZ  = 2
#         clog    = 16
#         E0      = 3.5e6*unyt.eV
#         Emin    = 50 * Te

#         # Thermal distribution with correct normalization
#         thdist  = run_tgo.getdist("5d", exi=True)
#         egridth = thdist.abscissa_edges("ekin")
#         simtime = np.diff(thdist.abscissa_edges("time")[-2:])
#         thermal = 2*np.sqrt(egridth/np.pi) * np.power(Te, -3.0/2) \
#             * np.exp(-egridth/Te) * simtime

#         vth   = np.sqrt(2*Te / m_e)
#         vcrit = vth * np.power( (3.0*np.sqrt(np.pi)/4.0) * (m_e / m_p) , 1/3.0)
#         Ecrit = (0.5 * m_a * vcrit * vcrit).to("eV")
#         ts    = ( 3 * np.sqrt( (2*np.pi * Te)**3 / m_e ) * eps_0**2
#                   * m_a / ( alphaZ**2 * unyt.e**4 * ne * clog ) ).to("s")

#         egridsd   = np.linspace(Te, 1.1*E0, 100)
#         heaviside = np.logical_and(egridsd <= E0, egridsd >= Emin)
#         slowing   = heaviside * ts / ( ( 1 + np.power(Ecrit/egridsd, 3.0/2) )
#                                        * 2 * egridsd )

#         # ts is slowing down rate which gives the slowing down time as
#         # t_sd = ts*log(v_0 / v_th) = 0.5*ts*log(E_0/E_th)
#         slowingdowntime = 0.5 * ts * np.log( E0 / Emin )

#         th_ekin  = np.zeros((3,))
#         th_pitch = np.zeros((3,))
#         for i, run in enumerate([run_tgo, run_tgcf, run_tgca]):
#             dist  = run.getdist("5d", exi=True)
#             edist = dist.integrate(
#                 r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[1:],
#                 charge=np.s_[:], pitch=np.s_[:], copy=True)
#             xdist = dist.integrate(
#                 r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[1:],
#                 charge=np.s_[:], ekin=np.s_[:], copy=True)
#             edist.plot(axes=h1)
#             xdist.plot(axes=h2)
#             th_pitch[i] = np.mean(run.getstate("pitch", state="end"))
#             th_ekin[i]  = np.mean(run.getstate("ekin",  state="end"))

#         sd_time  = np.zeros((3,))
#         sd_pitch = np.zeros((3,))
#         for i, run in enumerate([run_sgo, run_sgcf, run_sgca]):
#             dist = run.getdist("5d", exi=True, ekin_edges=egridsd)
#             edist = dist.integrate(
#                 r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[:],
#                 charge=np.s_[:], pitch=np.s_[:], copy=True)
#             xdist = dist.integrate(
#                 r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[:],
#                 charge=np.s_[:], ekin=np.s_[:], copy=True)
#             edist.plot(axes=h3)
#             xdist.plot(axes=h4)
#             sd_pitch[i] = np.mean(run.getstate("pitch", state="end"))
#             sd_time[i]  = np.mean(run.getstate("time",  state="end"))

#         # These are multiplied with marker number to get correct normalization
#         Nmrkth = run_tgo.getstate("ids").size
#         Nmrksd = run_sgo.getstate("ids").size
#         h1.plot(egridth, thermal*Nmrkth, color="black")
#         h3.plot(egridsd, slowing*Nmrksd, color="black")

#         print("Test Coulomb collisions:")
#         passed = True

#         print("  Thermal final energy and pitch")
#         print("  GO            %1.1e      %1.2f" % (th_ekin[0], th_pitch[0]))
#         print("  GCF           %1.1e      %1.2f" % (th_ekin[1], th_pitch[1]))
#         print("  GCA           %1.1e      %1.2f" % (th_ekin[2], th_pitch[2]))
#         print("  Expected      %1.1e      %1.2f" % (Te, 0.0))
#         if np.amax(np.abs(Te.v - th_ekin)) > 2e3 or \
#            np.amax(np.abs(0.0 - th_pitch)) > 0.5:
#             print("  (Failed)")
#             passed = False
#         print("")
#         print("  Slowing-down final time  and  pitch")
#         print("  GO            %1.1e      %1.2f" % (sd_time[0], sd_pitch[0]))
#         print("  GCF           %1.1e      %1.2f" % (sd_time[1], sd_pitch[1]))
#         print("  GCA           %1.1e      %1.2f" % (sd_time[2], sd_pitch[2]))
#         print("  Expected      %1.1e      %1.2f" % (slowingdowntime, 0.0))
#         if np.amax(np.abs(slowingdowntime.v  - sd_time))  > 2e-3 or \
#            np.amax(np.abs(0.0 - sd_pitch)) > 0.2:
#             print("  (Failed)")
#             passed = False

#         fig = a5plt.figuredoublecolumn(3/2)
#         ax = fig.add_subplot(1,1,1)
#         dist = run_tgo.getdist("5d")
#         dist.integrate(
#             r=np.s_[:], phi=np.s_[:], z=np.s_[:], time=np.s_[:],
#             charge=np.s_[:])
#         dist.plot(axes=ax)

#         return passed