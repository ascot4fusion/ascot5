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
        mrk["pitch"][:]  = 0.7
        mrk["r"][:]      = 7.8
        a5.data.create_input("gc", **mrk)
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]))

        a5.data.create_input("N0_3D")
        a5.data.create_input("Boozer")
        a5.data.create_input("MHD_STAT")
        a5.data.create_input("asigma_loc")

        name = a5.data.options.active.new(
            ENDCOND_MAX_MILEAGE=1.5e-4, DIST_NBIN_TIME=1, DIST_MIN_CHARGE=1.5,
            DIST_MAX_CHARGE=2.5, DIST_NBIN_CHARGE=1, DIST_MIN_R=4, DIST_MAX_R=8,
            DIST_MIN_Z=-2, DIST_MAX_Z=2,
            ENABLE_DIST_6D=0, ENABLE_DIST_RHO6D=1, ENABLE_DIST_COM=0,
            DIST_NBIN_PHI=1, DIST_NBIN_R=50, DIST_NBIN_Z=100,
            DIST_NBIN_PPE=50, DIST_NBIN_PPA=50,
            DIST_NBIN_THETA=2,  DIST_NBIN_RHO=10,
            ORBITWRITE_NPOINT=6000)
        a5.data.options[name].activate()

        subprocess.run(["./ascot5_main", "--in=unittest.h5"])

    @classmethod
    def tearDownClass(cls):
        super(TestAscot, cls).tearDownClass()
        subprocess.run(["rm", "-f", "unittest.h5"])

    def test_moments(self):
        a5 = Ascot("unittest.h5")
        ordinates = ["density", "chargedensity", "energydensity", "pressure",
                     "toroidalcurrent", "parallelcurrent", "powerdep",
                     "electronpowerdep", "ionpowerdep"]

        a5.input_init(bfield=True, plasma=True)
        dist    = a5.data.active.getdist("5d")
        mom     = a5.data.active.getdist_moments(dist, *ordinates)
        rhodist = a5.data.active.getdist("rho5d")
        rhomom  = a5.data.active.getdist_moments(rhodist, *ordinates)

        re = dist.abscissa_edges("r")
        ze = dist.abscissa_edges("z")
        v1, v2 = np.meshgrid(re[1:]**2 - re[:-1]**2, ze[1:] - ze[:-1])
        vol = np.pi * v1 * v2

        rho, r, z, phi, psi, weight, time,  \
        charge, energy, vphi, vnorm, vpar,  \
        mass, p, pitch, bnorm, bphi, ppar = \
            a5.data.active.getorbit("rho", "r", "z", "phi", "psi", "weight",
                                    "mileage", "charge", "ekin", "vphi",
                                    "vnorm", "vpar", "mass", "pnorm", "pitch",
                                    "bnorm", "bphi", "ppar")
        ei     = a5.data.active.getstate("ekin", state="ini")
        ef, tf = a5.data.active.getstate("ekin", "mileage", state="end")
        dt = np.diff(time, prepend=0)

        k, nu = a5.input_eval_collcoefs(
            mass[0], charge[0], r, phi, z, time, vnorm, "k", "nu",
            grid=False)
        dEtot_d = p * dt * np.sum(k, axis=0)
        dEele_d = p * dt * k[0,:]
        dEion_d = p * dt * np.sum(k[1:,:], axis=0)
        dPsi    = np.diff(psi)
        dPsi    = np.append(dPsi, 0)
        nu      = np.sum(nu, axis=0)
        k       = np.sum(k, axis=0)
        dppar   = mass * k * pitch - p * pitch * nu
        Pphi    = ppar * r * (bphi/bnorm) + charge * psi
        dPphi   = np.diff(Pphi)
        dPphi   = np.append(dPphi, 0)

        a5.input_free()

        mrkdist = {}
        mrkdist["density"]          = weight * dt
        mrkdist["chargedensity"]    = weight * charge * dt
        mrkdist["energydensity"]    = weight * ( energy * dt ).to("J")
        mrkdist["pressure"]         = weight * (mass * vnorm**2 * dt).to("J")/3
        mrkdist["toroidalcurrent"]  = weight * ( charge * vphi * dt ).to("A*m")
        mrkdist["parallelcurrent"]  = weight * ( charge * vpar * dt ).to("A*m")
        mrkdist["powerdep"]         = weight * dEtot_d.to("W")
        mrkdist["electronpowerdep"] = weight * dEele_d.to("W")
        mrkdist["ionpowerdep"]      = weight * dEion_d.to("W")
        mrkdist["jxbtorque"]        = weight * -charge * dPsi
        mrkdist["colltorque"]       = weight * (r*dppar*(bphi/bnorm)*dt).to("J")
        mrkdist["canmomtorque"]     = weight * -charge *dPphi

        print(weight[0]*tf)
        print(((ef-ei)*weight[0]/unyt.s).to("W"))
        for o in ordinates:
            a1 = np.sum(rhomom.ordinate(o) * rhomom.volume)
            a2 = np.sum(mom.ordinate(o) * mom.volume)
            a3 = np.sum(mrkdist[o])
            print(o, a1, a2, a3)
        """
        for w in mrkdist:
            prt = np.histogram2d(
                r, z, [dist.abscissa_edges("r"), dist.abscissa_edges("z")],
                weights=w)[0]
            rhoprt = np.histogram(
                rho, rhodist.abscissa_edges("rho"), weights=w)[0]

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(1,3,1)
        plt.title("Test")
        pm = ax1.pcolormesh(
            dist.abscissa_edges("r").v, dist.abscissa_edges("z").v,
            prt.T / vol, shading="flat")
        plt.colorbar(pm)
        ax1.set_aspect("equal", adjustable="box")
        ax2 = fig.add_subplot(1,3,2)
        plt.title("Ascot <3")
        a5.data.active.plotdist_moments(mom, ordinate, axes=ax2)

        ax3 = fig.add_subplot(1,3,3)
        pm = ax3.plot(rhodist.abscissa("rho").v,
                      rhoprt / np.sum(np.sum(rhomom.volume, axis=1), axis=1))
        a5.data.active.plotdist_moments(rhomom, ordinate, axes=ax3)
        plt.show()
        """

if __name__ == '__main__':
    unittest.main()
