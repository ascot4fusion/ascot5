# def check_biosaw(self):
#         """Check and run BioSaw test results.
#         """
#         # Define coil geometry (Circular coil centered at x=r0,y=0,z=0)
#         nseg = 200
#         r0   = 2.0
#         rad  = 1.0
#         coilxyz = np.zeros((nseg,3))
#         coilxyz[:,0] = r0  + rad * np.cos(np.linspace(0, 2*np.pi, nseg))
#         coilxyz[:,1] = 0.0 + rad * np.sin(np.linspace(0, 2*np.pi, nseg))
#         coilxyz[:,2] = 0.0

#         # Coordinates where field is evaluated
#         npnt = 100
#         r   = np.array([r0])
#         phi = np.array([0.0])
#         z   = np.linspace(-10.0, 10.0, npnt)

#         # Calculate analytical and numerical fields
#         Icoil = 1.0
#         br, bphi, bz = self.ascot.biosaw.calculate(
#             r, phi, z, coilxyz, current=Icoil)
#         br   = np.squeeze(br)
#         bphi = np.squeeze(bphi)
#         bz   = np.squeeze(bz)

#         br0   = np.zeros((npnt,))
#         bphi0 = np.zeros((npnt,))
#         bz0   = 0.5 * unyt.mu_0.v * rad**2 * Icoil \
#             / np.power(rad**2 + z**2, 3.0/2)

#         fig = plt.figure()
#         gs = GridSpec(4, 1, figure=fig, hspace=0.9, top=0.95)
#         ax1 = fig.add_subplot(gs[0,0])
#         ax2 = fig.add_subplot(gs[1,0])

#         ax1.plot(bz, color="C0")
#         ax1.plot(bz0, color="black", ls="--")
#         ax2.plot(br**2  + bphi**2, color="C0")
#         ax2.plot(br0**2 + bphi0**2, color="black", ls="--")

#         passed = True
#         print("Test BioSaw:")
#         if any(br**2 + bphi**2 > 1e-40) or any(np.abs(bz - bz0) > 1e-10):
#             print("- Field on axis of a cylindrical loop (FAILED)")
#             passed = False
#         else:
#             print("- Field on axis of a cylindrical loop")

#         # Tokamak field test
#         coilxyz = np.zeros((nseg,3))
#         coilxyz[:,0] = r0  + rad * np.cos(np.linspace(0, 2*np.pi, nseg))
#         coilxyz[:,1] = 0.0
#         coilxyz[:,2] = 0.0 + rad * np.sin(np.linspace(0, 2*np.pi, nseg))

#         npnt = 100
#         r   = np.linspace(r0-rad*0.8, r0+rad*0.8, npnt)
#         phi = np.linspace(0, 2*np.pi, npnt+1)[:-1]
#         z   = np.array([0.0])

#         Icoil = -1.0
#         _, bphi, _ = self.ascot.biosaw.calculate(
#             r, phi, z, coilxyz, current=Icoil, revolve=1)

#         brad  = np.squeeze(bphi[:,0,0])
#         bphi  = np.squeeze(bphi[int(npnt/2),:,0])
#         brad0 = -unyt.mu_0.v * npnt * Icoil / ( 2 * np.pi * r )
#         bphi0 = bphi[0] * np.ones((npnt,))

#         if any(np.abs(bphi - bphi0) > 1e-10) or \
#            any(np.abs(brad - brad0) > 1e-7):
#             print("- Tokamak field (FAILED)")
#             passed = False
#         else:
#             print("- Tokamak field")

#         ax3 = fig.add_subplot(gs[2,0])
#         ax4 = fig.add_subplot(gs[3,0])
#         ax3.plot(r, brad, color="C0")
#         ax3.plot(r, brad0, color="black", ls="--")
#         ax4.plot(phi, bphi, color="C0")
#         ax4.plot(phi, bphi0, color="black", ls="--")

#         ax1.set_xlabel("$z$ [m]")
#         ax1.set_ylabel("$B_z$ [T]")
#         ax2.set_xlabel("$z$ [m]")
#         ax2.set_ylabel("$|B_{xy}|$ [T]")

#         ax3.set_xlabel("$r$ [m]")
#         ax3.set_ylabel("$B_phi$ [m]")
#         ax4.set_xlabel("$phi$ [rad]")
#         ax4.set_ylabel("$B_phi$ [m]")
#         return passed