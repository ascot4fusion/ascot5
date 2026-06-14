# def init_afsi(self):
#         """Initialize data for AFSI test.
#         """
#         if hasattr(self.ascot.data.bfield, PhysTest.tag_afsi_thermal):
#             warnings.warn("Inputs already present: Test AFSI")
#             return

#         bgs = self.ascot.data.create_input("bfield analytical iter circular",
#                                            dryrun=True)

#         # DT-plasma
#         nrho  = 4
#         rho   = np.array([0, 1, 1+1e-1, 10])
#         edens = 2e20 * np.ones((nrho, 1))
#         etemp = 1e4  * np.ones((nrho, 1))
#         idens = 1e20 * np.ones((nrho, 2))
#         itemp = 1e4  * np.ones((nrho, 1))

#         edens[rho>1]   = 1
#         idens[rho>1,:] = 1

#         pls = {
#             "nrho" : nrho, "nion" : 2, "rho" : rho,
#             "anum" : np.array([2, 3]), "znum" : np.array([1, 1]),
#             "mass" : np.array([2.014, 3.016]), "charge" : np.array([1, 1]),
#             "edensity" : edens, "etemperature" : etemp,
#             "idensity" : idens, "itemperature" : itemp,
#             "vtor": 0*rho,
#         }

#         for tag in [PhysTest.tag_afsi_thermal, PhysTest.tag_afsi_beamthermal,
#                     PhysTest.tag_afsi_beambeam]:
#             self.ascot.data.create_input("B_GS", **bgs, desc=tag)
#             self.ascot.data.create_input("plasma_1D", **pls, desc=tag)

#     def run_afsi(self):
#         """Run AFSI tests.
#         """
#         if hasattr(self.ascot.data, PhysTest.tag_afsi_thermal):
#             warnings.warn("Results already present: Test AFSI")
#             return

#         r = np.linspace(5.9, 6.1, 2)
#         phi = np.linspace(0, 360, 2)
#         z = np.linspace(-0.1, 0.1, 2)
#         ppar = np.linspace(-1.0e-19, 1.0e-19, 401) * unyt.kg *unyt.m / unyt.s
#         pperp = np.linspace(0, 1.0e-19, 401) * unyt.kg *unyt.m / unyt.s
#         ekin = np.linspace(0,7e6,100)
#         pitch = np.linspace(-1,1,20)
#         self._activateinputs(PhysTest.tag_afsi_thermal)
#         self.ascot.afsi.thermal(
#             "DT_He4n", nmc=10**6, r=r, phi=phi, z=z,
#             #ekin1=ekin, ekin2=ekin, pitch1=pitch, pitch2=pitch,
#             ppar1=ppar, ppar2=ppar, pperp1=pperp, pperp2=pperp,
#             )
#         self.ascot.data.active.set_desc(PhysTest.tag_afsi_thermal)

#         dist = self.ascot.data[PhysTest.tag_afsi_thermal].getdist("prod1")
#         ppa, ppe = np.meshgrid(dist.abscissa("ppar"), dist.abscissa("pperp"),
#                                indexing="ij")
#         vol = np.diff(dist.abscissa_edges("r")[:2]) \
#             * np.diff(dist.abscissa_edges("z")[:2]) \
#             * np.diff(dist.abscissa_edges("phi")[:2]) \
#             * np.diff(dist.abscissa_edges("time")[:2]) \
#             * np.diff(dist.abscissa_edges("charge")[:2])

#         def maxwellian(mass):
#             mass *= unyt.amu
#             T = 1e4*unyt.eV
#             nprt = 1e20*2*np.pi*6*0.2*0.2
#             return nprt * ppe * np.exp( -(ppe**2 + ppa**2) / (2*mass*T) ) \
#                 / ( np.power(mass*T, 3/2.0) * vol * np.sqrt(2*np.pi) )

#         dist._distribution[0,0,0,:,:,0,0] = maxwellian(2.014)

#         self._activateinputs(PhysTest.tag_afsi_beamthermal)
#         self.ascot.afsi.beamthermal(
#             "DT_He4n", dist, nmc=10**6,
#             ppar1=ppar, ppar2=ppar, pperp1=pperp, pperp2=pperp,
#             #ekin1=ekin, ekin2=ekin, pitch1=pitch, pitch2=pitch,
#             )
#         self.ascot.data.active.set_desc(PhysTest.tag_afsi_beamthermal)

#         dist = self.ascot.data[PhysTest.tag_afsi_thermal].getdist("prod1")
#         dist._distribution[0,0,0,:,:,0,0] = maxwellian(2.014)

#         dist1 = self.ascot.data[PhysTest.tag_afsi_thermal].getdist("prod1")
#         dist1._distribution[0,0,0,:,:,0,0] = maxwellian(3.016)

#         self._activateinputs(PhysTest.tag_afsi_beambeam)
#         self.ascot.afsi.beambeam(
#             "DT_He4n", dist, dist1, nmc=10**6,
#             ppar1=ppar, ppar2=ppar, pperp1=pperp, pperp2=pperp,
#             #ekin1=ekin, ekin2=ekin, pitch1=pitch, pitch2=pitch,
#             )
#         self.ascot.data.active.set_desc(PhysTest.tag_afsi_beambeam)

#     def check_afsi(self):
#         """Check AFSI tests.
#         """
#         alphadist = self.ascot.data[PhysTest.tag_afsi_thermal].getdist(
#             "prod1", exi=True, ekin_edges=np.linspace(0, 4.5e6, 100),
#             pitch_edges=100)

#         alphadist.integrate(
#             time=np.s_[:], charge=np.s_[:], r=np.s_[:], z=np.s_[:],
#             phi=np.s_[:])
#         ekindist = alphadist.integrate(copy=True, pitch=np.s_[:])
#         xidist   = alphadist.integrate(copy=True, ekin=np.s_[:])

#         ekin  = alphadist.abscissa("ekin")
#         pitch = alphadist.abscissa("pitch")

#         emean = (3.537e6 + 1e4 * 2) * unyt.eV
#         gauss = np.exp(-(ekin - emean)**2 / (4*1e4*emean*4.001/(4.001+1.008))) \
#             / unyt.eV

#         pdens = 1e20*unyt.m**(-3)*1e20*unyt.m**(-3)*1.1e-16*unyt.cm**3
#         vol = 6*unyt.m*np.pi*2*0.2*unyt.m*0.2*unyt.m

#         fig = plt.figure()
#         ax1 = fig.add_subplot(2,1,1)
#         ax2 = fig.add_subplot(2,1,2)

#         ekindist.plot(axes=ax1)
#         ax1.plot(ekin, pdens * vol * gauss / np.trapz(gauss, ekin))
#         xidist.plot(axes=ax2)
#         ax2.plot(pitch, 0.5 * pdens * vol * np.ones(pitch.shape))

#         eana  = pdens * vol * gauss / np.trapz(gauss,ekin)
#         xiana = 0.5 * pdens * vol * np.ones(pitch.shape)
#         ekin1 = ekindist.distribution()
#         xi1   = xidist.distribution()

#         alphadist = self.ascot.data[PhysTest.tag_afsi_beamthermal].getdist(
#             "prod1", exi=True, ekin_edges=np.linspace(0, 4.5e6, 100),
#             pitch_edges=100)

#         alphadist.integrate(
#             time=np.s_[:], charge=np.s_[:], r=np.s_[:], z=np.s_[:],
#             phi=np.s_[:])
#         ekindist = alphadist.integrate(copy=True, pitch=np.s_[:])
#         xidist   = alphadist.integrate(copy=True, ekin=np.s_[:])

#         ekindist.plot(axes=ax1)
#         xidist.plot(axes=ax2)
#         ekin2 = ekindist.distribution()
#         xi2   = xidist.distribution()

#         alphadist = self.ascot.data[PhysTest.tag_afsi_beambeam].getdist(
#             "prod1", exi=True, ekin_edges=np.linspace(0, 4.5e6, 100),
#             pitch_edges=100)

#         alphadist.integrate(
#             time=np.s_[:], charge=np.s_[:], r=np.s_[:], z=np.s_[:],
#             phi=np.s_[:])
#         ekindist = alphadist.integrate(copy=True, pitch=np.s_[:])
#         xidist   = alphadist.integrate(copy=True, ekin=np.s_[:])
#         ekin3 = ekindist.distribution()
#         xi3   = xidist.distribution()

#         ekindist.plot(axes=ax1)
#         xidist.plot(axes=ax2)

#         passed = True
#         print("Test AFSI:")
#         err1 = np.square(np.subtract(eana, ekin1)).mean()
#         err2 = np.square(np.subtract(xiana, xi1)).mean()
#         if err1.v > 1e23 or err2.v > 2e35:
#             passed = False
#             print("Thermal %e %e (FAILED)" % (err1.v, err2.v))
#         else:
#             print("Thermal %e %e" % (err1.v, err2.v))

#         err1 = np.square(np.subtract(eana, ekin2)).mean()
#         err2 = np.square(np.subtract(xiana, xi2)).mean()
#         if err1.v > 1e23 or err2.v > 2e35:
#             passed = False
#             print("Beam-Thermal %e %e (FAILED)" % (err1.v, err2.v))
#         else:
#             print("Beam-Thermal %e %e" % (err1.v, err2.v))

#         err1 = np.square(np.subtract(eana, ekin3)).mean()
#         err2 = np.square(np.subtract(xiana, xi3)).mean()
#         if err1.v > 1e23 or err2.v > 2e35:
#             passed = False
#             print("Beam-Beam %e %e (FAILED)" % (err1.v, err2.v))
#         else:
#             print("Beam-Beam %e %e" % (err1.v, err2.v))

#         return passed