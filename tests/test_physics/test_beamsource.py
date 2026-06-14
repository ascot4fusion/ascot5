# def init_bbnbi(self):
#         """Initialize data for the BBNBI deposition test.
#         """
#         if hasattr(self.ascot.data.options, PhysTest.tag_bbnbi):
#             warnings.warn("Inputs already present: Test BBNBI")
#             return
#         init = self.ascot.data.create_input
#         bdata = {"bxyz":np.array([0.0,1.0,0.0]), "jacobian":1.0*np.zeros((3,3)),
#                  "rhoval":0.5}
#         init("B_TC", **bdata, desc=PhysTest.tag_bbnbi)

#         # Hydrogen plasma
#         nrho  = 4
#         rho   = np.array([0, 1, 1+1e-8, 10])
#         vtor  = np.zeros((nrho,1))
#         edens = 1e19 * np.ones((nrho, 1))
#         etemp = 2e3  * np.ones((nrho, 1))
#         idens = 1e19 * np.ones((nrho, 1))
#         itemp = 2e3  * np.ones((nrho, 1))
#         edens[rho>1]   = 1
#         idens[rho>1,:] = 1

#         pls = {
#             "nrho" : nrho, "nion" : 1, "rho" : rho, "vtor" : vtor,
#             "anum" : np.array([1]), "znum" : np.array([1]),
#             "mass" : np.array([1.014]), "charge" : np.array([1]),
#             "edensity" : edens, "etemperature" : etemp,
#             "idensity" : idens, "itemperature" : itemp}
#         init("plasma_1D", **pls, desc=PhysTest.tag_bbnbi)

#         opt = Opt.get_default()
#         opt.update({
#             "ENABLE_DIST_5D":0,
#             "DIST_MIN_R":4.3,       "DIST_MAX_R":8.3,      "DIST_NBIN_R":50,
#             "DIST_MIN_PHI":0,       "DIST_MAX_PHI":360,    "DIST_NBIN_PHI":1,
#             "DIST_MIN_Z":-2.0,      "DIST_MAX_Z":2.0,      "DIST_NBIN_Z":50,
#             "DIST_MIN_PPA":-1.e-19, "DIST_MAX_PPA":1.e-19, "DIST_NBIN_PPA":100,
#             "DIST_MIN_PPE":0,       "DIST_MAX_PPE":1.e-19, "DIST_NBIN_PPE":50,
#             "DIST_MIN_TIME":0,      "DIST_MAX_TIME":1.0,   "DIST_NBIN_TIME":1,
#             "DIST_MIN_CHARGE":0.5,  "DIST_MAX_CHARGE":1.5, "DIST_NBIN_CHARGE":1,
#         })
#         init("opt", **opt, desc=PhysTest.tag_bbnbi)

#         halo_f = 0.15
#         div_v  = 10e-3
#         div_h  =  5e-3
#         halo_v = 30e-3
#         halo_h = 20e-3
#         inj = Injector(
#             ids=1, anum=1, znum=1, mass=1.0,
#             energy=4e4, efrac=np.array([1.0,0.0,0.0]), power=1.0,
#             divh=div_h, divv=div_v, divhalofrac=halo_f,
#             divhaloh=halo_h, divhalov=halo_v, nbeamlet=1,
#             beamletx=np.array([10.0]), beamlety=np.array([0.0]),
#             beamletz=np.array([0.0]),
#             beamletdx=np.array([-1.0]), beamletdy=np.array([0.0]),
#             beamletdz=np.array([0.0]))
#         init("nbi", **{"ninj":1, "injectors":[inj]}, desc=PhysTest.tag_bbnbi)

#     def run_bbnbi(self):
#         """Run BBNBI.
#         """
#         if hasattr(self.ascot.data, PhysTest.tag_bbnbi):
#             warnings.warn("Results already present: Test BBNBI")
#             return
#         self._activateinputs(PhysTest.tag_bbnbi)
#         self.ascot.data.nbi[PhysTest.tag_bbnbi].activate()
#         subprocess.call(
#             ["./../../build/bbnbi5", "--in=testascot.h5",
#              "--d="+PhysTest.tag_bbnbi],
#             stdout=subprocess.DEVNULL)
#         self.ascot = Ascot(self.ascot.file_getpath())

#     def check_bbnbi(self):
#         """Check the BBNBI deposition test.
#         """
#         run  = self.ascot.data[PhysTest.tag_bbnbi]
#         r, z, x, y, s = run.getstate("r", "z", "x", "y", "mileage", mode="prt")

#         omega_x = np.mod(np.arctan2(y.v, x.v-10) + 2*np.pi, 2*np.pi)-np.pi
#         omega_y = np.mod(np.arctan2(z.v, r.v-10) + 2*np.pi, 2*np.pi)-np.pi

#         inj = run.nbi.read()["injectors"][0]
#         div_h  = inj.divh
#         div_v  = inj.divv
#         halo_h = inj.divhaloh
#         halo_v = inj.divhalov
#         halo_f = inj.divhalofrac

#         # The plasma parameters and this value are taken from "Neutral beam
#         # stopping and emission in fusion plasmas I: deuterium beams",
#         # H Anderson et al 2000 Plasma Phys. Control. Fusion 42 781, p. 10
#         mfp = 1.0/(1e19*1.1e-13)

#         def gaussian(omega, omegadiv):
#             return np.exp( -(omega/omegadiv)**2 ) / omegadiv

#         def profile(omega, div, halo_div):
#             return halo_f * gaussian(omega, halo_div) / np.sqrt(np.pi) \
#                 + (1.0 - halo_f) * gaussian(omega, div) / np.sqrt(np.pi)

#         def ionization(distance, mfp):
#             mfp *= 1e6
#             return np.exp(-distance/mfp) / mfp

#         dist_edges = np.linspace(0,7,100)
#         distance = ( dist_edges[1:] + dist_edges[:-1] ) / 2

#         omega_edges = np.linspace(-0.04*np.pi/2, 0.04*np.pi/2, 181)
#         omega = ( omega_edges[1:] + omega_edges[:-1] ) / 2

#         mrk1 = np.histogram(omega_x, omega_edges, density=True)[0]
#         mrk2 = np.histogram(omega_y, omega_edges, density=True)[0]
#         mrk3 = np.histogram(s.v*1e6,  dist_edges, density=True)[0]

#         ana1 = profile(omega, div_h, halo_h)
#         ana2 = profile(omega, div_v, halo_v)
#         ana3 = ionization(distance, mfp)

#         p = scipy.optimize.curve_fit(profile, omega, mrk1, [div_h, halo_h])[0]
#         fit1 = profile(omega, p[0], p[1])
#         div_h_fit  = p[0]
#         halo_h_fit = p[1]

#         p = scipy.optimize.curve_fit(profile, omega, mrk2, [div_v, halo_v])[0]
#         fit2 = profile(omega, p[0], p[1])
#         div_v_fit  = p[0]
#         halo_v_fit = p[1]

#         p = scipy.optimize.curve_fit(ionization, distance, mrk3, [mfp])[0]
#         fit3 = ionization(distance, p[0])
#         mfp_fit = p[0]

#         fig = plt.figure(figsize=(5,5))
#         gs = GridSpec(3, 1, figure=fig, hspace=0.6)
#         ax1 = fig.add_subplot(gs[0,0])
#         a5plt.mesh1d(omega_edges, mrk1, axes=ax1)
#         ax1.plot(omega, ana1, color="black")
#         ax1.plot(omega, fit1, color="C3", ls="--")
#         ax1.set_xlabel(r"Horizontal divergence $\omega$ [rad]")

#         ax2 = fig.add_subplot(gs[1,0])
#         a5plt.mesh1d(omega_edges, mrk2, axes=ax2)
#         ax2.plot(omega, ana2, color="black")
#         ax2.plot(omega, fit2, color="C3", ls="--")
#         ax2.set_xlabel(r"Vertical divergence $\omega$ [rad]")

#         ax3 = fig.add_subplot(gs[2,0])
#         a5plt.mesh1d(dist_edges, mrk3, axes=ax3)
#         ax3.plot(distance, ana3, color="black")
#         ax3.plot(distance, fit3, color="C3", ls="--")
#         ax3.set_xlabel(r"Time until ionized [$\mu$s]")

#         ax2.set_ylabel("Probability density")

#         from matplotlib.lines import Line2D
#         legend_elements = [
#             Line2D([0], [0], color="C0", lw=2, label="Markers"),
#             Line2D([0], [0], color="black", lw=2, label="Analytical"),
#             Line2D([0], [0], color="C3", lw=2, ls="--", label="Fit"),
#         ]
#         ax3.legend(handles=legend_elements, loc="upper right", frameon=False)

#         passed = True
#         if np.abs(div_h  - div_h_fit)  > 1e-3: passed = False
#         if np.abs(div_v  - div_v_fit)  > 5e-3: passed = False
#         if np.abs(halo_h - halo_h_fit) > 5e-3: passed = False
#         if np.abs(halo_v - halo_v_fit) > 5e-3: passed = False
#         if np.abs(mfp    - mfp_fit)    > 1e-7: passed = False

#         failed = "" if passed else "(FAILED)"

#         print("Test BBNBI:")
#         print("          div_h     div_v     halo_h    halo_v    mfp")
#         print("  Actual: %.3e %.3e %.3e %.3e %.3e" \
#               % (div_h[0], div_v[0], halo_h[0], halo_v[0], mfp) )
#         print("  Fitted: %.3e %.3e %.3e %.3e %.3e %s" \
#               % (div_h_fit, div_v_fit, halo_h_fit, halo_v_fit, mfp_fit, failed))

#         return passed