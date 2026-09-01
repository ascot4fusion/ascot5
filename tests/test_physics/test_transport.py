"""Verify that the reproduces classical and neoclassical transport."""
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
def test_classical(ascot: Ascot, method, inspect, plot):

    if not inspect:
        parameters = SimulationOptions.from_dict(
            physics={
                "enable_orbit_following": True,
                "enable_coulomb_collisions": True,
            },
            endconditions={
                "activate_simulation_time_limits": True,
                "max_mileage": 1e-5,
            },
        )

        if method == "go":
            parameters.simulation.mode = "gyro-orbit"
            parameters.simulation.timestep=1e-10
        else:
            parameters.simulation.mode = "guiding-center"
            parameters.simulation.timestep = 1e-9
            if method == "gc-fixedstep":
                parameters.simulation.enable_adaptive = False

        ascot.data.create_guidingcentermarker(
            species="proton",
            r=5*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=1e3*unyt.eV,
            pitch=1.0 - 2.0 * np.random.rand(1000,),
            )

        # Pure electron plasma to avoid proton-proton collisions
        ascot.data.create_plasmalinear1d(
            species=["H"],
            rhogrid=np.linspace([0, 2]),
            ni=1e0*unyt.m**-3,
            ne=1e22*unyt.m**-3,
            Ti=1e3*unyt.eV,
            )

        b_scan = 1.0 / np.sqrt(np.linspace(0.01, 1.0, 6))
        for bscale in b_scan:
            ascot.data.create_bfieldcartesian(
            bxyz=np.array([bscale, 0., 0.])*unyt.T,
            jacobian=np.full((3, 3), 0.)*unyt.T/unyt.m,
            axisrz=np.array([1., 0.])*unyt.m,
            rhoval=1.,
            )
            ascot.simulate(params=parameters)


@pytest.mark.parametrize(
    "method",
    ["go", "gc-fixedstep", "gc-adaptive"],
)
def test_neoclassical_tokamak(ascot: Ascot, method, inspect, plot):
    if not inspect:
        ascot.data.create_guidingcentermarker(
            species="proton",
            r=7.2*unyt.m,
            z=0.0*unyt.m,
            phi=0.0*unyt.deg,
            ekin=1e3*unyt.eV,
            pitch=1.0 - 2.0 * np.random.rand(100,),
            charge=1.*unyt.e,
            )

        density = np.array([8.9e17, 8.9e18, 3.0e19, 7.0e20, 2.0e21, 2.0e22]) * unyt.m**-3
        max_mileage = np.array([10e-4, 8e-4, 6e-4, 4e-4, 2e-4, 1e-4])
        for tmax, ni in zip(max_mileage, density):
            parameters = SimulationOptions(
                enable_orbit_following=True,
                enable_coulomb_collisions=True,
                activate_simulation_time_limits=True,
                max_mileage=tmax,
                )

            if method == "go":
                parameters.simulation.simulation_mode = 1
                parameters.simulation.timestep=2e-10
            else:
                parameters.simulation.simulation_mode = 2
                parameters.simulation.enable_adaptive = True
                parameters.simulation.timestep = 2e-10
                if method == "gc-fixedstep":
                    parameters.simulation.enable_adaptive = False

            # Pure proton plasma to avoid proton-electron collisions
            ascot.data.create_plasmalinear1d(
                species=["H"],
                rhogrid=np.linspace([0, 2]),
                ni=ni,
                ne=1e0*unyt.m**-3,
                Ti=1e3*unyt.eV,
                )
            ascot.simulate(params=parameters)



#     def check_classical(self):
#         """Check classical transport test.
#         """
#         nscan = 0
#         while PhysTest.tag_classical_go + str(nscan) in \
#               self.ascot.data.bfield.ls(show=False):
#             nscan += 1

#         # Initialize plots
#         fig = a5plt.figuredoublecolumn(3/2)
#         ax = fig.add_subplot(1,1,1)

#         # Numerical values
#         ndim = 2 # Diffusion happens on 2D plane
#         bnorm = np.zeros((nscan,)) * unyt.T
#         Dgo   = np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgcf  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgca  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         for i in range(nscan):
#             run_go  = self.ascot.data[PhysTest.tag_classical_go  + str(i)]
#             run_gcf = self.ascot.data[PhysTest.tag_classical_gcf + str(i)]
#             run_gca = self.ascot.data[PhysTest.tag_classical_gca + str(i)]

#             yi, zi, ti = run_go.getstate("y", "z", "mileage", state="ini")
#             yf, zf, tf = run_go.getstate("y", "z", "mileage", state="end")
#             Dgo[i] = np.mean( ( (yi - yf)**2 + (zi - zf)**2 ) / (tf - ti) ) \
#                 / (2*ndim)
#             yi, zi, ti = run_gcf.getstate("y", "z", "mileage", state="ini")
#             yf, zf, tf = run_gcf.getstate("y", "z", "mileage", state="end")
#             Dgcf[i] = np.mean( ( (yi - yf)**2 + (zi - zf)**2 ) / (tf - ti) ) \
#                 / (2*ndim)
#             yi, zi, ti = run_gca.getstate("y", "z", "mileage", state="ini")
#             yf, zf, tf = run_gca.getstate("y", "z", "mileage", state="end")
#             Dgca[i] = np.mean( ( (yi - yf)**2 + (zi - zf)**2 ) / (tf - ti) ) \
#                 / (2*ndim)

#             bnorm[i] = np.sqrt(np.sum(run_go.bfield.read()["bxyz"]**2))

#         # Analytical
#         clog = 13.4
#         ekin = run_go.getstate("ekin")[0]
#         self.ascot.input_init(run=run_go.get_qid(), bfield=True, plasma=True)
#         ne, Te = self.ascot.input_eval(
#             6.2*unyt.m, 0*unyt.deg, 0*unyt.m, 0*unyt.s, "ne", "te")
#         self.ascot.input_free()
#         collfreq = physlib.collfreq_ie(unyt.mp, unyt.e, ne, Te, clog)
#         rhog = physlib.gyrolength(unyt.mp, 1*unyt.e, ekin, 0.0, bnorm).to("m")
#         Dana = collfreq * rhog**2 / 2

#         # Plotting
#         ax.scatter(1/bnorm**2, Dgo)
#         ax.scatter(1/bnorm**2, Dgcf)
#         ax.scatter(1/bnorm**2, Dgca)
#         ax.plot(1/bnorm**2, Dana, color="black")
#         ax.set_xlabel(r"$1/B^{2}$ [1/T$^2$]")
#         ax.set_ylabel(r"Diffusion [m$^2$/s]")

#         print("Test classical transport:")
#         passed = True
#         k0 = _linfit(1/bnorm**2, Dana)
#         k1 = _linfit(1/bnorm**2, Dgo)
#         k2 = _linfit(1/bnorm**2, Dgcf)
#         k3 = _linfit(1/bnorm**2, Dgca)

#         print(" slope (expected): %1.3f" % k0)
#         f = ""
#         if(np.abs(k0-k1) > 1e-2): passed=False; f = "(FAILED)"
#         print("               GO: %1.3f %s" % (k1, f))
#         f = ""
#         if(np.abs(k0-k1) > 1e-2): passed=False; f = "(FAILED)"
#         print("              GCF: %1.3f %s" % (k2, f))
#         f = ""
#         if(np.abs(k0-k1) > 1e-2): passed=False; f = "(FAILED)"
#         print("              GCA: %1.3f %s" % (k3, f))

#         return passed

#     def init_neoclassical(self):
#         """Initialize data for the neoclassical transport test.
#         """
#         if hasattr(self.ascot.data.options, PhysTest.tag_neoclassical_go+"0"):
#             warnings.warn("Inputs already present: Test neoclass. transport")
#             return
#         init = self.ascot.data.create_input

#         # Ion densities to be scanned
#         ni = np.array([8.9e17, 8.9e18, 3.0e19, 7.0e20, 2.0e21, 2.0e22])

#         

#     def run_neoclassical(self):
#         """Run neoclassical transport test.
#         """
#         if hasattr(self.ascot.data, PhysTest.tag_neoclassical_go+"0"):
#             warnings.warn("Results already present: Test neoclass. transport")
#             return
#         for tag in [PhysTest.tag_neoclassical_go, PhysTest.tag_neoclassical_gcf,
#                     PhysTest.tag_neoclassical_gca]:
#             i = 0
#             while tag + str(i) in self.ascot.data.bfield.ls(show=False):
#                 self._activateinputs(tag+str(i))
#                 self._runascot(tag+str(i))
#                 i += 1

#     def check_neoclassical(self):
#         """Check neoclassical transport test.
#         """
#         nscan = 0
#         items = self.ascot.data.bfield.ls(show=False)
#         while PhysTest.tag_neoclassical_go + str(nscan) in items:
#             nscan += 1

#         # Evaluate Ti and R_omp(rho) as these are needed
#         run_go = self.ascot.data[PhysTest.tag_neoclassical_go  + "0"]
#         self.ascot.input_init(run=run_go.get_qid(), bfield=True, plasma=True)
#         Ti         = self.ascot.input_eval(
#             6.2*unyt.m, 0*unyt.deg, 0*unyt.m, 0*unyt.s, "ti1")
#         rhoomp     = np.linspace(0, 1, 100)
#         romp, zomp = self.ascot.input_rhotheta2rz(
#             rhoomp, 0*unyt.rad, 0*unyt.rad, 0*unyt.s)
#         self.ascot.input_free()

#         # Initialize plots
#         fig = a5plt.figuredoublecolumn(3/2)
#         ax = fig.add_subplot(1,1,1)
#         ax.set_xscale("log")
#         ax.set_yscale("log")

#         # Numerical values
#         ni    = np.zeros((nscan,)) / unyt.m**3
#         Dgo   = np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgoerr= np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgcf  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgcferr= np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgca  = np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         Dgcaerr= np.zeros((nscan,)) * unyt.m**2 / unyt.s
#         for i in range(nscan):
#             run_go  = self.ascot.data[PhysTest.tag_neoclassical_go  + str(i)]
#             run_gcf = self.ascot.data[PhysTest.tag_neoclassical_gcf + str(i)]
#             run_gca = self.ascot.data[PhysTest.tag_neoclassical_gca + str(i)]

#             ri, ti = run_go.getstate("rho", "mileage", state="ini")
#             rf, tf = run_go.getstate("rho", "mileage", state="end")
#             ri = np.interp(ri, rhoomp, romp)
#             rf = np.interp(rf, rhoomp, romp)
#             Dgo[i] = 0.5 * np.mean( (rf - ri)**2 / (tf - ti) )
#             Dgoerr[i] = np.sqrt( (0.5 * np.var( (rf - ri)**2 / (tf - ti) )) \
#                                  / ri.size )

#             ri, ti = run_gcf.getstate("rho", "mileage", state="ini")
#             rf, tf = run_gcf.getstate("rho", "mileage", state="end")
#             ri = np.interp(ri, rhoomp, romp)
#             rf = np.interp(rf, rhoomp, romp)
#             Dgcf[i] = 0.5 * np.mean( (rf - ri)**2 / (tf - ti) )
#             Dgcferr[i] = np.sqrt( (0.5 * np.var( (rf - ri)**2 / (tf - ti) )) \
#                                   / ri.size )

#             ri, ti = run_gca.getstate("rho", "mileage", state="ini")
#             rf, tf = run_gca.getstate("rho", "mileage", state="end")
#             ri = np.interp(ri, rhoomp, romp)
#             rf = np.interp(rf, rhoomp, romp)
#             Dgca[i] = 0.5 * np.mean( (rf - ri)**2 / (tf - ti) )
#             Dgcaerr[i] = np.sqrt( (0.5 * np.var( (rf - ri)**2 / (tf - ti) )) \
#                                   / ri.size )

#             ni[i] = run_go.plasma.read()["idensity"][0, 0]

#         r0    = run_go.getstate("r")[0]
#         axisr = run_go.bfield.read()["raxis"][0] * unyt.m
#         eps   = (r0 - axisr) / axisr
#         qfac  = 1.7 # Safety factor at r0 was verified numerically
#         bnorm = run_go.bfield.read()["bphi0"][0] * unyt.T

#         ekin   = run_go.getstate("ekin")[0]
#         omegat = physlib.bouncefrequency(unyt.me, ekin, r0, axisr, qfac)
#         rhog   = physlib.gyrolength(unyt.me, 1*unyt.e, ekin, 0.0, bnorm).to("m")

#         clog     = 15
#         density  = np.power( 10, np.linspace(np.log10(ni[0]) - 1,
#                                              np.log10(ni[-1]) + 1, 50) )
#         density /= unyt.m**3
#         collfreq = physlib.collfreq_ie(unyt.mp, unyt.e, density, Ti, clog) \
#             * ( unyt.mp / unyt.me )
#         veff = collfreq / omegat
#         # Add intermediate values needed for plotting a continuous curve
#         veff = np.append(veff, [1, np.power(eps, 3.0/2.0)])
#         veff.sort()

#         Dps = qfac**2 * veff * omegat * rhog**2 / 2
#         Dp  = 0.5 * qfac**2 * omegat * rhog**2 * np.ones(veff.shape)
#         Db  = np.power(eps, -3.0/2.0) * Dps

#         # x coordinate for plotting the numerical coefficients
#         collfreq = physlib.collfreq_ie(unyt.mp, unyt.e, ni, Ti, clog) \
#             * ( unyt.mp / unyt.me )
#         veff_x = collfreq / omegat

#         i1 = np.nonzero(veff==np.power(eps, 3.0/2.0))[0][0]
#         i2 = np.nonzero(veff==1)[0][0]
#         ax.plot(veff[i2:],     Dps[i2:],    color="black")
#         ax.plot(veff[i1:i2+1], Dp[i1:i2+1], color="black")
#         ax.plot(veff[:i1+1],   Db[:i1+1],   color="black")

#         ax.errorbar(veff_x, Dgo, yerr=Dgoerr,  linestyle="none", marker="*")
#         ax.errorbar(veff_x*1.01, Dgcf, yerr=Dgcferr, linestyle="none",
#                     marker="o")
#         ax.errorbar(veff_x*1.09, Dgca, yerr=Dgcaerr, linestyle="none",
#                     marker="^")

#         ax.set_xlabel(r"Effective collisionality $\nu^*$")
#         ax.set_ylabel(r"Diffusion [m$^2$/s]")

#         ax.plot([1,1], [1e-6, 1e-1], color="gray")
#         ax.plot([eps**(3.0/2), eps**(3.0/2)], [1e-6, 1e-1], color="gray")
#         ax.text(10**(-1.9), 10**(-5.7), r"$\nu^*=\epsilon^{3/2}$", fontsize=10,
#            bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
#         ax.text(10**(-0.25), 10**(-5.7), r"$\nu^*=1$", fontsize=10,
#                bbox={'facecolor':'white', 'edgecolor':'none', 'pad':0})
#         ax.text(10**(-2.5), 10**(-3.7), r"$D_{B}$",  fontsize=10)
#         ax.text(10**(-0.8), 10**(-2.7), r"$D_{P}$",  fontsize=10)
#         ax.text(10**(0.4),  10**(-1.8), r"$D_{PS}$", fontsize=10)

#         ax.set_xlim(1e-4, 1e2)
#         ax.set_ylim(1e-6, 1e-1)

#         f = "(FAILED)"
#         print("Test neoclassical transport:")
#         print("                   GO       GCF      GCA      analytical")
#         passed = True
#         idx = veff_x <= np.power(eps, 3.0/2.0)
#         k0 = _linfit(veff[:i1+1], Db[:i1+1])
#         k1 = _linfit(veff_x[idx], Dgo[idx])
#         k2 = _linfit(veff_x[idx], Dgcf[idx])
#         k3 = _linfit(veff_x[idx], Dgca[idx])
#         f = ""
#         if np.amax(np.abs(np.array([k1,k2,k3]) - k0)) > 2e-2:
#             f = "(FAILED)"
#             passed=False
#         print("  Banana regime    %1.2e %1.2e %1.2e %1.2e %s"
#               % (k1, k2, k3, k0, f))

#         idx = np.logical_and(veff_x >= np.power(eps, 3.0/2.0), veff_x <=1)
#         k0 = 0.0
#         k1 = _linfit(veff_x[idx], Dgo[idx])
#         k2 = _linfit(veff_x[idx], Dgcf[idx])
#         k3 = _linfit(veff_x[idx], Dgca[idx])
#         f = ""
#         if np.amax(np.abs(np.array([k1,k2,k3]) - k0)) > 3e-3:
#             f = "(FAILED)"
#             passed=False
#         print("  Plateau regime   %1.2e %1.2e %1.2e %1.2e %s"
#               % (k1, k2, k3, k0, f))

#         idx = veff_x >=1
#         k0 = _linfit(veff[i2:],   Dps[i2:])
#         k1 = _linfit(veff_x[idx], Dgo[idx])
#         k2 = _linfit(veff_x[idx], Dgcf[idx])
#         k3 = _linfit(veff_x[idx], Dgca[idx])
#         f = ""
#         if np.amax(np.abs(np.array([k1,k2,k3]) - k0)) > 5e-4:
#             f = "(FAILED)"
#             passed=False
#         print("  Pfirsch-Schlüter %1.2e %1.2e %1.2e %1.2e %s"
#               % (k1, k2, k3, k0, f))

#         return passed