"""Unit tests for distribution moments.
"""
import numpy as np
import unittest
import subprocess
import unyt

from a5py import Ascot, AscotInitException, AscotIOException

class TestAscot(unittest.TestCase):
    """Class for testing distribution moments.
    """

    @classmethod
    def setUpClass(cls):
        """Create and run a test case.

        Assumes ascot5_main is located in the same folder.
        """
        super(TestAscot, cls).setUpClass()
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("options tutorial")
        a5.data.create_input("bfield analytical iter circular")
        a5.data.create_input("wall rectangular")
        a5.data.create_input("plasma flat")

        from a5py.ascot5io.marker import Marker
        mrk = Marker.generate("gc", n=1, species="alpha")
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.8
        mrk["r"][:]      = 7.8
        a5.data.create_input("gc", **mrk)
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]))

        a5.data.create_input("N0_3D")
        a5.data.create_input("Boozer")
        a5.data.create_input("MHD_STAT")
        a5.data.create_input("asigma_loc")

        name = a5.data.options.active.new(
            ENDCOND_MAX_MILEAGE=0.5e-4, DIST_MIN_CHARGE=1.5,
            DIST_MAX_CHARGE=2.5, DIST_NBIN_PPE=10, DIST_NBIN_PPA=10,
            DIST_NBIN_PHI=1, DIST_NBIN_R=50, DIST_NBIN_Z=100, DIST_NBIN_TIME=1,
            DIST_NBIN_THETA=50, DIST_MIN_R=4, DIST_MAX_R=8, DIST_MIN_Z=-2,
            DIST_MAX_Z=2, DIST_NBIN_RHO=10,
            ENABLE_DIST_6D=0, ENABLE_DIST_RHO6D=1, ENABLE_DIST_COM=0,
            ORBITWRITE_NPOINT=2000)
        a5.data.options[name].activate()

        subprocess.run(["./ascot5_main", "--in=unittest.h5"])

    @classmethod
    def tearDownClass(cls):
        super(TestAscot, cls).tearDownClass()
        subprocess.run(["rm", "-f", "unittest.h5"])

    def test_moments(self):
        a5 = Ascot("unittest.h5")
        ordinate = "electronpowerdep"

        a5.input_init(bfield=True, plasma=True)
        dist = a5.data.active.getdist("5d", exi=True)
        mom = a5.data.active.getdist_moments(dist, ordinate)

        rhodist = a5.data.active.getdist("rho5d", exi=True)
        rhomom = a5.data.active.getdist_moments(rhodist, ordinate)
        #print(np.sum(rhodist.distribution()))
        #print(np.sum(rhomom.ordinate(ordinate)))
        #print(dist.abscissa("pitch"))
        #print(dist.abscissa("ekin").to("eV"))

        rho, r, z, weight, time, charge, energy, vphi = a5.data.active.getorbit(
            "rho", "r", "z", "weight", "mileage", "charge", "ekin", "vphi")
        a5.input_free()

        if ordinate == "density":
            weight *= np.diff(time, prepend=0)
        elif ordinate == "chargedensity":
            weight *= np.diff(time, prepend=0) * charge
        elif ordinate == "energydensity":
            weight *= np.diff(time, prepend=0) * energy
        elif ordinate == "toroidalcurrent":
            weight *= np.diff(time, prepend=0) * charge * vphi
        elif ordinate == "electronpowerdep":
            # This weight is just a placeholder
            weight *= np.diff(time, prepend=0) * charge * vphi
        prt = np.histogram2d(
            r, z, [dist.abscissa_edges("r"), dist.abscissa_edges("z")],
            weights=weight)[0]
        rhoprt = np.histogram(
            rho, rhodist.abscissa_edges("rho"), weights=weight)[0]

        re = dist.abscissa_edges("r")
        ze = dist.abscissa_edges("z")
        v1, v2 = np.meshgrid(re[1:]**2 - re[:-1]**2, ze[1:] - ze[:-1])
        vol = np.pi * v1 * v2

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(1,3,1)
        pm = ax1.pcolormesh(
            dist.abscissa_edges("r").v, dist.abscissa_edges("z").v,
            prt.T / vol, shading="flat")
        plt.colorbar(pm)
        ax1.set_aspect("equal", adjustable="box")
        ax2 = fig.add_subplot(1,3,2)
        a5.data.active.plotdist_moments(mom, ordinate, axes=ax2)

        ax3 = fig.add_subplot(1,3,3)
        pm = ax3.plot(rhodist.abscissa("rho").v,
                      rhoprt / np.sum(np.sum(rhomom.volume, axis=1), axis=1))
        a5.data.active.plotdist_moments(rhomom, ordinate, axes=ax3)

        plt.show()

if __name__ == '__main__':
    unittest.main()
