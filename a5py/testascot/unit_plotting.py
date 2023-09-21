"""Test plotting routines.

The plotting routines obtain the data using routines that are tested elsewhere,
so here we can focus on making sure that the plots look correct. These tests
only fail if there was an exception but otherwise the actual "testing" is to
look at the figures.
"""
import numpy as np
import unittest
import subprocess
import matplotlib.pyplot as plt
import unyt

from a5py import Ascot
import a5py.routines.plotting as a5plt
from a5py.ascot5io.marker import Marker

class TestPlotting(unittest.TestCase):

    def tearDown(self):
        try:
            a5 = Ascot("unittest.h5")
            a5.simulation_free()
        except:
            pass
        subprocess.run(["rm", "-f", "unittest.h5"])

    def initandrunsim(self, simtime):
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("options tutorial")
        a5.data.create_input("bfield analytical iter circular")
        a5.data.create_input("wall rectangular")
        a5.data.create_input("plasma flat")
        a5.data.create_input("E_TC", exyz=np.array([0,0,0]))
        a5.data.create_input("N0_3D")
        a5.data.create_input("Boozer")
        a5.data.create_input("MHD_STAT")
        a5.data.create_input("asigma_loc")

        a5.simulation_initinputs()

        mrk = Marker.generate("gc", n=100, species="alpha")
        mrk["weight"][:] = 10**(10 +10*np.random.rand(100,))
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(100,)
        mrk["r"][:]      = np.linspace(6.2, 8.5, 100)
        a5.simulation_initmarkers(**mrk)

        opt = a5.data.options.active.read()
        opt.update(ENDCOND_MAX_MILEAGE=simtime)
        a5.simulation_initoptions(**opt)
        return a5.simulation_run()

    def test_input_plotrz(self):
        a5 = Ascot("unittest.h5", create=True)
        a5.data.create_input("bfield analytical iter circular", splines=True,
                             axisymmetric=True)
        a5.data.create_input("plasma flat")
        a5.input_init(bfield=True, plasma=True)

        fig = plt.figure(figsize=(30,40))

        r = np.linspace( 3.0, 9.0,  50) * unyt.m
        z = np.linspace(-2.0, 2.0, 100) * unyt.m

        ax = fig.add_subplot(2,4,1)
        a5.input_plotrz(r, z, "bnorm", clim=[None, None], cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,5)
        a5.input_plotrz(r, z, "bnorm", clim=[0, 10], cmap="Greys", axes=ax)

        ax = fig.add_subplot(2,4,2)
        a5.input_plotrz(r, z, "br", clim=[None, None], cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,6)
        a5.input_plotrz(r, z, "br", clim=[-1, 1], cmap="PuOr", axes=ax)

        ax = fig.add_subplot(2,4,3)
        a5.input_plotrz(r, z, "log ne", clim=[None, None], cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,7)
        a5.input_plotrz(r, z, "log ne", clim=[1e10, 1e22],
                        cmap="Greys", axes=ax)

        ax = fig.add_subplot(2,4,4)
        a5.input_plotrz(r, z, "log divb", clim=[None, None],
                        cmap=None, axes=ax)

        ax = fig.add_subplot(2,4,8)
        a5.input_plotrz(r, z, "log divb",clim=[-1e-14, 1e-14],
                        cmap="PuOr",axes=ax)

        a5.input_free()
        plt.show(block=False)

    def test_plotstate_scatter(self):
        vr = self.initandrunsim(1e-3)
        fig = plt.figure(figsize=(30,40))

        ax = fig.add_subplot(2,3,1)
        vr.plotstate_scatter("end rho", "log end ekin",
                             axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,2)
        vr.plotstate_scatter("ini rho", "log diff ekin",
                             axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,3)
        vr.plotstate_scatter("log end x", "log end y", c="log diff phi",
                             axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,4, projection="3d")
        vr.plotstate_scatter("end x", "end y", "end z",
                             axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,5, projection="3d")
        vr.plotstate_scatter("end x", "end y", "end z", c="end rho",
                             axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,6, projection="3d")
        vr.plotstate_scatter("log end x", "log end y", "log end z",
                             c="log diff phi", axesequal=True, axes=ax)

        plt.show(block=False)

    def test_plotstate_histogram(self):
        vr = self.initandrunsim(1e-4)
        fig = plt.figure(figsize=(30,40))

        ax = fig.add_subplot(2,3,1)
        vr.plotstate_histogram("end phimod", axes=ax)

        ax = fig.add_subplot(2,3,2)
        vr.plotstate_histogram("log end phi", xbins=10, endcond="WALL", axes=ax)

        ax = fig.add_subplot(2,3,3)
        vr.plotstate_histogram("end ekin", xbins=np.linspace(0,4e6,100),
                               logscale=True, weight=True, axes=ax)

        ax = fig.add_subplot(2,3,4)
        vr.plotstate_histogram("end phimod", "end thetamod", axes=ax)

        ax = fig.add_subplot(2,3,5)
        vr.plotstate_histogram("log end phi", "log end theta", endcond="WALL",
                               axes=ax)

        ax = fig.add_subplot(2,3,6)
        vr.plotstate_histogram("end phimod", "end thetamod", logscale=True,
                               weight=True, axes=ax)

        plt.show(block=False)

    def test_plotorbit_trajectory(self):
        vr = self.initandrunsim(1e-4)
        fig = plt.figure(figsize=(30,40))

        ax = fig.add_subplot(2,3,1)
        vr.plotorbit_trajectory("r", "z",
                                axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,2)
        vr.plotorbit_trajectory("mileage", "reldiff ekin", c="rho",
                                axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,3)
        vr.plotorbit_trajectory("log mileage", "log reldiff ekin",
                                c="log rho",
                                axesequal=False, axes=ax)

        ax = fig.add_subplot(2,3,4, projection="3d")
        vr.plotorbit_trajectory("x", "y", "z",
                                axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,5, projection="3d")
        vr.plotorbit_trajectory("x", "y", "z", c="diff ekin",
                                axesequal=True, axes=ax)

        ax = fig.add_subplot(2,3,6, projection="3d")
        vr.plotorbit_trajectory("log x", "log y", "log z", c="log diff ekin",
                                axesequal=True, axes=ax)

        plt.show(block=False)

if __name__ == '__main__':
    unittest.main()
