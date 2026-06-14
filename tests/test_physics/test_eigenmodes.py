# def init_mhd(self):
#         """Initialize data for the MHD test.

#         https://link.springer.com/article/10.1007/s41614-018-0022-9
#         """
#         if hasattr(self.ascot.data.options, PhysTest.tag_mhd_go):
#             warnings.warn("Inputs already present: Test MHD")
#             return
#         init = self.ascot.data.create_input

#         # Options
#         opt = Opt.get_default()
#         opt.update({
#             "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
#             "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_SIMTIMELIM" : 1,
#             "ENDCOND_LIM_SIMTIME" : 1e-5, "ENABLE_ORBIT_FOLLOWING" : 1,
#             "ENABLE_MHD" : 1, "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
#             "ORBITWRITE_INTERVAL" : 1e-9, "ORBITWRITE_NPOINT" : 10000
#         })
#         init("opt", **opt, desc=PhysTest.tag_mhd_go)
#         opt.update({
#             "SIM_MODE" : 2, "ENABLE_ADAPTIVE": 0,
#             "FIXEDSTEP_USERDEFINED" : 1e-11
#         })
#         init("opt", **opt, desc=PhysTest.tag_mhd_gcf)
#         opt.update({
#             "ENABLE_ADAPTIVE" : 1, "FIXEDSTEP_USERDEFINED" : 1e-11,
#             "ADAPTIVE_TOL_ORBIT" : 1e-14, "ADAPTIVE_MAX_DRHO" : 0.1,
#             "ADAPTIVE_MAX_DPHI" : 10
#         })
#         init("opt", **opt, desc=PhysTest.tag_mhd_gca)

#         # Use field line markers
#         mrk = Marker.generate("gc", n=2, species="electron")
#         mrk["r"][:]      = np.array([7.0, 8.0])
#         mrk["phi"][:]    = 0
#         mrk["z"][:]      = 0
#         mrk["pitch"][:]  = np.array([0.4, 0.9])
#         mrk["energy"][:] = 1e6

#         # ITER-like field with boozer data and field line markers
#         for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
#                     PhysTest.tag_mhd_gca]:
#             qid = init("bfield_analytical_iter_circular", splines=True,
#                        desc=tag)
#             init("gc", **mrk, desc=tag)

#         qid = self.ascot.data.bfield[qid].get_qid()
#         self.ascot.input_init(bfield=qid)
#         bzr = init("boozer_tokamak", nint=100000, dryrun=True)
#         for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
#                     PhysTest.tag_mhd_gca]:
#             init("Boozer", **bzr, desc=tag)

#         mhd = {
#             "nmode":1, "nmodes":np.array([2]), "mmodes":np.array([3]),
#             "amplitude":np.array([1e-2]), "omega":np.array([1e6]),
#             "phase":np.array([np.pi/4]), "nrho":100, "rhomin":0.1, "rhomax":0.99
#         }

#         rhogrid = np.linspace(mhd["rhomin"], mhd["rhomax"], mhd["nrho"])
#         mhd["alpha"] = np.tile(
#             np.exp( -(rhogrid-0.85)**2/0.1 ), (mhd["nmode"],1)).T
#         mhd = self.ascot.data.create_input(
#             "mhd consistent potentials", which="Phi", mhd=mhd, dryrun=True)
#         self.ascot.input_free(bfield=True)
#         for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
#                     PhysTest.tag_mhd_gca]:
#             self.ascot.data.create_input("MHD_STAT", **mhd, desc=tag)

#     def run_mhd(self):
#         """Run MHD test.
#         """
#         if hasattr(self.ascot.data, PhysTest.tag_mhd_go):
#             warnings.warn("Results already present: Test MHD")
#             return
#         for tag in [PhysTest.tag_mhd_go, PhysTest.tag_mhd_gcf,
#                     PhysTest.tag_mhd_gca]:
#             self._activateinputs(tag)
#             self._runascot(tag)

#     def check_mhd(self):
#         """Check MHD test.
#         """
#         run_go  = self.ascot.data[PhysTest.tag_mhd_go]
#         run_gcf = self.ascot.data[PhysTest.tag_mhd_gcf]
#         run_gca = self.ascot.data[PhysTest.tag_mhd_gca]

#         mhd = run_go.mhd.read()
#         self.ascot.input_init(run=run_go.get_qid(),
#                               bfield=True, boozer=True, mhd=True)

#         # Initialize plots
#         fig = a5plt.figuredoublecolumn(3/2)
#         ax1a = fig.add_subplot(3,1,1)
#         ax1b = ax1a.twinx()
#         ax2a = fig.add_subplot(3,1,2)
#         ax2b = ax2a.twinx()
#         ax3a = fig.add_subplot(3,1,3)
#         ax3b = ax3a.twinx()

#         ax3a.set_xlabel("Mileage [s]")
#         ax2a.set_ylabel("Relative error")

#         err1 = [None] * 3
#         err2 = [None] * 3
#         c = ["C0", "C1", "C2"]
#         for i, run in enumerate([run_go, run_gcf, run_gca]):
#             ekin, charge, ptor, bphi, r, phi, z, t, ids = run.getorbit(
#                 "ekin", "charge", "ptor", "bphi", "r", "phi", "z",
#                 "time", "ids")
#             Phi, alpha, br0, bphi0, bz0, er, ephi, ez = self.ascot.input_eval(
#                 r, phi, z, t, "phieig", "alphaeig", "br", "bphi", "bz",
#                 "mhd_er", "mhd_ephi", "mhd_ez")
#             epar  = ( er * br0 + ephi * bphi0 + ez * bz0 ) / \
#                 np.sqrt(br0**2 + bphi0**2 + bz0**2)
#             enorm = np.sqrt(er**2 + ephi**2 + ez**2)

#             ctor = ptor + r * charge * alpha * bphi
#             H = ekin + Phi * charge
#             P = (mhd["omega"]/unyt.s) * ctor / mhd["nmodes"]
#             K = H - P

#             H.convert_to_mks()
#             P.convert_to_mks()
#             K.convert_to_mks()

#             idx = ids == 1
#             ax1a.plot(t[idx], epar[idx], color=c[i])
#             ax1b.plot(t[idx], enorm[idx], ls="--", color=c[i])
#             err_e1 = np.amax(np.abs((epar/enorm)[idx]))
#             err = K[ids==1]/K[ids==1][0] - 1
#             ax2b.plot(t[idx], H[idx] - H[idx][0], ls="--", color=c[i])
#             ax2b.plot(t[idx], -P[idx] + P[idx][0], ls=":", color=c[i])
#             ax2a.plot(t[idx], K[idx] - K[idx][0], color=c[i])
#             err1[i] = np.amax(np.abs(err))

#             idx = ids == 2
#             ax1a.plot(t[idx], epar[idx], color=c[i])
#             ax1b.plot(t[idx], enorm[idx], ls="--", color=c[i])
#             err_e2 = np.amax(np.abs((epar/enorm)[idx]))
#             err = K[idx]/K[idx][0] - 1
#             ax3b.plot(t[idx], H[idx] - H[idx][0], ls="--", color=c[i])
#             ax3b.plot(t[idx], -P[idx] + P[idx][0], ls=":", color=c[i])
#             ax3a.plot(t[idx], K[idx] - K[idx][0], color=c[i])
#             err2[i] = np.amax(np.abs(err))

#         self.ascot.input_free()

#         passed = True
#         print("Test MHD:")
#         fail = ""
#         if err_e1 > 2e-5 or err_e2 > 2e-5:
#             fail = "(FAILED)"
#             passed = False
#         print("Error in Epar: %e %e %s" % (err_e1, err_e2, fail))

#         fail = ""
#         if err1[0] > 1e-6 or err2[0] > 1e-3:
#             fail = "(FAILED)"
#             passed = False
#         print("Error in H - P (GO):  %e %e %s" % (err1[0], err2[0], fail))

#         fail = ""
#         if err1[1] > 1e-7 or err2[1] > 1e-3:
#             fail = "(FAILED)"
#             passed = False
#         print("Error in H - P (GCF): %e %e %s" % (err1[1], err2[1], fail))

#         fail = ""
#         if err1[2] > 1e-9 or err2[2] > 1e-6:
#             fail = "(FAILED)"
#             passed = False
#         print("Error in H - P (GCA): %e %e %s" % (err1[2], err2[2], fail))
#         return passed