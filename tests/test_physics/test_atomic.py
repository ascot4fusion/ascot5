# def init_atomic(self):
#         """Initialize data for the atomic reaction test.
#         """
#         if hasattr(self.ascot.data.options, PhysTest.tag_atomic_cx):
#             warnings.warn("Inputs already present: Test atomic reaction")
#             return
#         init = self.ascot.data.create_input

#         # Options
#         opt = Opt.get_default()
#         opt.update({
#             "SIM_MODE" : 1, "FIXEDSTEP_USE_USERDEFINED" : 1,
#             "ENDCOND_SIMTIMELIM" : 1, "ENDCOND_LIM_SIMTIME" : 2e-2,
#             "ENABLE_ORBIT_FOLLOWING" : 1, "ENABLE_ATOMIC" : 1
#         })
#         opt.update({
#             "FIXEDSTEP_USERDEFINED" : 1e-9, "ENDCOND_NEUTRALIZED" : 1
#         })
#         init("opt", **opt, desc=PhysTest.tag_atomic_cx)
#         opt.update({
#             "FIXEDSTEP_USERDEFINED" : 1e-11, "ENDCOND_NEUTRALIZED" : 0,
#             "ENDCOND_IONIZED" : 1, "ENDCOND_LIM_SIMTIME" : 2e-6
#         })
#         init("opt", **opt, desc=PhysTest.tag_atomic_ionz)

#         mrk = Marker.generate("gc", n=1000, species="deuterium")
#         mrk["r"][:]      = 6.2
#         mrk["phi"][:]    = 0.0
#         mrk["z"][:]      = 0.0
#         mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)
#         mrk["energy"][:] = 1e5
#         mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
#         init("gc", **mrk, desc=PhysTest.tag_atomic_cx)

#         mrk = Marker.generate("gc", n=1000, species="deuterium")
#         mrk["charge"][:] = 0
#         mrk["r"][:]      = 6.2
#         mrk["phi"][:]    = 0.0
#         mrk["z"][:]      = 0.0
#         mrk["zeta"][:]   = 2 * np.pi * np.random.rand(mrk["n"],)
#         mrk["energy"][:] = 1e5
#         mrk["pitch"][:]  = 1.0 - 2.0 * np.random.rand(mrk["n"],)
#         init("gc", **mrk, desc=PhysTest.tag_atomic_ionz)

#         for tag in [PhysTest.tag_atomic_cx, PhysTest.tag_atomic_ionz]:
#             # Uniform magnetic field
#             btc = {"bxyz":np.array([0,1.0,0]), "jacobian":np.zeros((3,3)),
#                    "rhoval":0.5}
#             init("B_TC", **btc, desc=tag)

#             # Plasma and neutral data
#             init("plasma_flat", anum=2, znum=1, mass=2.0135532,
#                  density=1e20, temperature=1e3, desc=tag)
#             init("neutral_flat", anum=2, znum=1, density=1e16, temperature=1e3,
#                  desc=tag)

#             # Atomic data
#             try:
#                 init("import_adas")
#             except Exception as err:
#                 print(err)
#                 print("Using analytical model instead.")
#                 init("asigma_chebyshev_cx_hh0", desc=tag)

#     def run_atomic(self):
#         """Run atomic reaction test.
#         """
#         if hasattr(self.ascot.data, PhysTest.tag_atomic_cx):
#             warnings.warn("Results already present: Test atomic reaction")
#             return
#         for tag in [PhysTest.tag_atomic_cx, PhysTest.tag_atomic_ionz]:
#             self._activateinputs(tag)
#             self._runascot(tag)

#     def check_atomic(self):
#         """Check atomic reaction test.
#         """
#         run_cx   = self.ascot.data[PhysTest.tag_atomic_cx]
#         run_ionz = self.ascot.data[PhysTest.tag_atomic_ionz]

#         # Calculate mean free paths from the cross sections
#         ma, anum, znum, r, phi, z, t, va = run_cx.getstate(
#             "mass", "anum", "znum", "r", "phi", "z", "time", "vnorm",
#             state="ini", mode="prt", ids=1)
#         self.ascot.input_init(run=run_cx.get_qid(), bfield=True, plasma=True,
#                               neutral=True, asigma=True)
#         sigmacx  = self.ascot.input_eval_atomiccoefs(
#             ma, anum[0], znum[0], r, phi, z, t, va,
#             reaction="charge-exchange")
#         sigmabms = self.ascot.input_eval_atomiccoefs(
#             ma, anum[0], znum[0], r, phi, z, t, va,
#             reaction="beamstopping")
#         self.ascot.input_free()
#         mfp_cx0  = va / sigmacx
#         mfp_bms0 = va / sigmabms

#         # Mean free paths from the simulation
#         vnorm1 = run_cx.getstate("vnorm", state="ini", mode="prt")
#         vnorm2, mil = run_cx.getstate("vnorm", "mileage", state="end",
#                                       mode="prt")
#         mfp_cx = np.mean(mil*vnorm2)

#         vnorm1 = run_ionz.getstate("vnorm", state="ini", mode="prt")
#         vnorm2, mil = run_ionz.getstate("vnorm", "mileage", state="end",
#                                         mode="prt")
#         mfp_bms = np.mean(mil*vnorm2)

#         passed = True
#         print("Test atomic:")
#         fail = ""
#         if np.abs(mfp_cx/mfp_cx0 - 1) > 0.1:
#             fail = "(FAILED)"
#             passed = False
#         print("Mean free path CX:  %.3e m (numerical) %.3e m (analytical) %s" %
#               (mfp_cx, mfp_cx0[0,0], fail))
#         fail = ""
#         if np.abs(mfp_bms/mfp_bms0 - 1) > 0.1:
#             fail = "(FAILED)"
#             passed = False
#         print("Mean free path BMS: %.3e m (numerical) %.3e m (analytical) %s" %
#               (mfp_bms, mfp_bms0[0,0], fail))
#         return passed
