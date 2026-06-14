# def init_boozer(self):
#         """Initialize data for the Boozer transformation test.
#         """
#         if hasattr(self.ascot.data.options, PhysTest.tag_boozer):
#             warnings.warn("Inputs already present: Test Boozer transformation")
#             return
#         init = self.ascot.data.create_input

#         # Options
#         opt = Opt.get_default()
#         opt.update({
#             "SIM_MODE" : 4, "ENDCOND_SIMTIMELIM" : 1,
#             "ENDCOND_MAX_MILEAGE" : 1e2/3e8, "ENABLE_ORBIT_FOLLOWING" : 1,
#             "ENABLE_ORBITWRITE" : 1, "ORBITWRITE_MODE" : 1,
#             "ORBITWRITE_INTERVAL" : 1e-1/3e8, "ORBITWRITE_NPOINT" : 10002
#         })

#         # Use field line markers
#         mrk = Marker.generate("fl", n=1)
#         mrk["r"][:]      = 8.0
#         mrk["phi"][:]    = 0
#         mrk["z"][:]      = 0
#         mrk["pitch"][:]  = 1.0

#         # Magnetic field is just some tokamak with Boozer data
#         b = init("bfield_analytical_iter_circular", dryrun=True)
#         b.update({"rmin" : 4, "rmax" : 8.5, "nr" : 120, "zmin" : -4,
#                   "zmax" : 4, "nz" : 200})

#         bphi = [5.3, -5.3, 5.3, -5.3]
#         psimult = [200, 200, -200, -200]
#         for i in range(4):
#             b.update({"bphi0" : bphi[i], "psimult" : psimult[i]})
#             out = B_2DS.convert_B_GS(**b)
#             qid = init("B_2DS", **out, desc=PhysTest.tag_boozer+str(i))
#             qid = self.ascot.data.bfield[qid].get_qid()
#             self.ascot.input_init(bfield=qid)
#             qid = init("boozer_tokamak", desc=PhysTest.tag_boozer+str(i),
#                        nint=100000)
#             self.ascot.input_free(bfield=True)

#             init("fl", **mrk, desc=PhysTest.tag_boozer+str(i))
#             init("opt", **opt, desc=PhysTest.tag_boozer+str(i))

#         # This perturbation is used only in the post-processing
#         mhd = {"nmode" : 1, "nmodes" : np.array([2]), "mmodes" : np.array([3]),
#                "amplitude" : np.array([1.0]), "omega" : np.array([1.0]),
#                "phase" : np.array([0.0]), "nrho" : 100, "rhomin" : 0.0,
#                "rhomax" : 1.0}
#         rhogrid = np.linspace(mhd["rhomin"], mhd["rhomax"], mhd["nrho"])
#         alpha   = np.exp( -(rhogrid-0.85)**2/0.1 )
#         phi     = alpha*0
#         mhd["phi"]   = np.tile(phi, (mhd["nmode"],1)).T
#         mhd["alpha"] = np.tile(alpha, (mhd["nmode"],1)).T
#         init("MHD_STAT", **mhd, desc=PhysTest.tag_boozer+"0")

#     def run_boozer(self):
#         """Run Boozer transformation test.
#         """
#         if hasattr(self.ascot.data, PhysTest.tag_boozer + "0"):
#             warnings.warn("Results already present: Test Boozer transformation")
#             return
#         for i in range(4):
#             self._activateinputs(PhysTest.tag_boozer + str(i))
#             self._runascot(PhysTest.tag_boozer + str(i))

#     def check_boozer(self):
#         """Check Boozer transformation test.
#         """
#         run1 = self.ascot.data[PhysTest.tag_boozer+"0"]
#         run2 = self.ascot.data[PhysTest.tag_boozer+"1"]
#         run3 = self.ascot.data[PhysTest.tag_boozer+"2"]
#         run4 = self.ascot.data[PhysTest.tag_boozer+"3"]

#         # Initialize plots
#         fig = a5plt.figuredoublecolumn(3/2)
#         ax1 = fig.add_subplot(3,2,1)
#         ax2 = fig.add_subplot(3,2,3)
#         ax3 = fig.add_subplot(3,2,5)
#         ax4 = fig.add_subplot(3,2,2)
#         ax5 = fig.add_subplot(3,2,4)

#         ax3.set_xlabel("Boozer theta [rad]")
#         ax1.set_ylabel("Error in Bpol [T]")
#         ax2.set_ylabel("Error in Bphi [T]")
#         ax3.set_ylabel("Jacobian")
#         ax4.set_ylabel("q")
#         ax5.set_ylabel("q")
#         ax5.set_xlabel("Orbit parameter s")

#         colors = ["C0", "C1", "C2", "C3"]
#         ip_out   = [None] * 4
#         bphi_out = [None] * 4
#         bpol_err = [None] * 4
#         bphi_err = [None] * 4
#         jac_err  = [None] * 4
#         q_err    = [None] * 4
#         for i, run in enumerate([run1, run2, run3, run4]):
#             self.ascot.input_init(run=run.get_qid(), bfield=True, boozer=True,
#                                   mhd=True)
#             r, phi, z, t, pol = run.getorbit("r", "phi", "z", "time", "theta")
#             theta, zeta, alpha, jacb2 = self.ascot.input_eval(
#                 r, phi, z, t, "theta", "zeta", "alphaeig", "bjacxb2")

#             # Results are plotted as a function of poloidal angle so find the
#             # "discontinuity" points at 0/2pi.
#             idx = np.nonzero(np.abs(np.diff(theta)) > np.pi)[0]
#             theta = theta.v
#             zeta  = zeta.v

#             # Numerical safety factor from field lines and field data
#             dz = np.diff(zeta)
#             dz[dz > np.pi]  = dz[dz > np.pi]  - 2*np.pi
#             dz[dz < -np.pi] = dz[dz < -np.pi] + 2*np.pi
#             dt = np.diff(theta)
#             dt[dt > np.pi]  = dt[dt > np.pi]  - 2*np.pi
#             dt[dt < -np.pi] = dt[dt < -np.pi] + 2*np.pi
#             qfac = dz/dt

#             rho = run.getstate("rho", state="ini")
#             q, I, g = self.ascot.input_eval_safetyfactor(rho)

#             # Evaluate gradients and field components so we can compare those
#             a, b, c = self.ascot.input_eval(
#                 r, phi, z, t, "dpsidr (bzr)", "dpsidphi (bzr)", "dpsidz (bzr)")
#             gradpsi = np.array([a, b/r, c]).T
#             a, b, c = self.ascot.input_eval(
#                 r, phi, z, t, "dthetadr", "dthetadphi", "dthetadz")
#             gradtheta = np.array([a, b/r, c]).T
#             a, b, c = self.ascot.input_eval(
#                 r, phi, z, t, "dzetadr", "dzetadphi", "dzetadz")
#             gradzeta = np.array([a, b/r, c]).T

#             br, bphi, bz = self.ascot.input_eval(
#                 r, phi, z, t, "br", "bphi", "bz")
#             self.ascot.input_free()

#             # Magnetic field vector from Boozer coordinates
#             bvec = -( np.mean(qfac) * np.cross(gradpsi, gradtheta)
#                       - np.cross(gradpsi, gradzeta) )

#             idx = np.nonzero(np.abs(np.diff(theta)) > np.pi)[0]
#             dbpol = np.sqrt((bvec[:,0] - br.v)**2 + (bvec[:,2] - bz.v)**2)
#             dbphi = bvec[:,1] - bphi.v
#             bnorm = br**2 + bz**2 + bphi**2

#             # Jacobian times B^2
#             jacb2_bzr = I + g*q

#             ip_out[i]   = 1 - 2 * (bz[0] > 0)
#             bphi_out[i] = 1 - 2 * (bphi[0] < 0)
#             bpol_err[i] = np.amax(np.abs(dbpol[1:idx[-1]]))
#             bphi_err[i] = np.amax(np.abs(dbphi[1:idx[-1]]))
#             jac_err[i]  = np.amax(np.abs(jacb2-jacb2_bzr.v))
#             q_err[i]    = np.amax(np.abs(q-qfac))

#             j0 = 1
#             j = idx[-1]
#             ax1.plot(theta[j0:j], dbpol[j0:j], color=colors[i])
#             ax2.plot(theta[j0:j], dbphi[j0:j], color=colors[i])
#             ax3.plot(theta[j0:j], np.abs(jacb2[j0:j]), color=colors[i])
#             ax3.plot(theta[np.array([j0,j])], np.abs([jacb2_bzr, jacb2_bzr]),
#                      color=colors[i], ls='--')
#             if qfac[0] > 0:
#                 ax4.plot(qfac, color=colors[i])
#                 ax4.plot([0,qfac.size-1], [q,q], color=colors[i], ls='--')
#             else:
#                 ax5.plot(qfac, color=colors[i])
#                 ax5.plot([0,qfac.size-1], [q,q], color=colors[i], ls='--')

#         passed = True
#         print("Test Boozer coordinate mapping:")
#         print("Bphi Ip | Delta Bpol  Delta Bphi  Delta Jacobian  Delta q")
#         for i in range(4):
#             fail = ""
#             if bpol_err[i] > 1e-15 or bphi_err[i] > 1e-3 or jac_err[i] > 1e-2\
#                or q_err[i] > 1e-4:
#                 fail = "(FAILED)"
#                 passed = False
#             print(" %2d  %2d    %.3e   %.3e   %.3e       %.3e %s" %
#                   (ip_out[i], bphi_out[i], bpol_err[i],
#                    bphi_err[i], jac_err[i], q_err[i], fail))
#         return passed