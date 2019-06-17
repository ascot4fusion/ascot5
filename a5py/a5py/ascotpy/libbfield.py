"""
Methods to evaluate quantities from magnetic field data.

File: libbfield.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot


import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class LibBfield(LibAscot):

    quantities = ["rho", "psi", "br", "bphi", "bz", "brdr", "brdphi", "brdz",
                  "bphidr", "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz",
                  "divergence", "axis", "bnorm",
                  "ripplecritlarmor"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["rho", "psi"]:
            out = self.eval_bfield(R, phi, z, t, evalrho=True)[quantity]

        if quantity in ["br", "bphi", "bz", "brdr", "brdphi", "brdz", "bphidr",
                        "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz"]:
            out = self.eval_bfield(R, phi, z, t, evalb=True)[quantity]

        if quantity == "divergence":
            out = self.eval_bfield(R, phi, z, t, evalb=True)
            out = out["br"]/R + out["brdr"] + out["bphidphi"]/R + out["bzdz"]
        if quantity == "axis":
            out = self.eval_bfield(R, phi, z, t, evalaxis=True)
        if quantity == "bnorm":
            out = self.eval_bfield(R, phi, z, t, evalb=True)
            out = np.sqrt( out["br"]*out["br"] + out["bphi"]*out["bphi"]
                           + out["bz"]*out["bz"] )
        if quantity == "ripplecritlarmor":
            pass

        assert out is not None, "Unknown quantity"

        return out


    def evaluateripple(self, R, z, t, nphi):
        nin = R.size
        phigrid = np.linspace(0, 2*np.pi, nphi+1)[:-1]
        R = np.meshgrid(R, phigrid, indexing="ij")[0]
        z = np.meshgrid(z, phigrid, indexing="ij")[0]
        t, phigrid = np.meshgrid(t, phigrid, indexing="ij")
        out = np.abs(self.eval_bfield(R.ravel(), phigrid.ravel(), z.ravel(),
                                      t.ravel(), evalb=True)["bphi"])
        out = out.reshape(nin, nphi)
        bmax = np.nanmax(out, axis=1)
        bmin = np.nanmin(out, axis=1)
        return (bmax - bmin) / (bmax + bmin)


    def evaluateripplewell(self, R, z, t, nphi):
        nin = R.size
        phigrid = np.linspace(0, 2*np.pi, nphi+1)[:-1]
        R = np.meshgrid(R, phigrid, indexing="ij")[0]
        z = np.meshgrid(z, phigrid, indexing="ij")[0]
        t, phigrid = np.meshgrid(t, phigrid, indexing="ij")
        out = self.eval_bfield(R.ravel(), phigrid.ravel(), z.ravel(),
                               t.ravel(), evalb=True)

        bnorm = np.sqrt( out["br"]*out["br"] + out["bphi"]*out["bphi"]
                         + out["bz"]*out["bz"] )
        bhat_r   = out["br"]   / bnorm
        bhat_phi = out["bphi"] / bnorm
        bhat_z   = out["bz"]   / bnorm
        dbdl_r   = (bhat_r * out["brdr"]   + bhat_phi * out["bphidr"]   + bhat_z * out["bzdr"]) * bhat_r
        dbdl_phi = (bhat_r * out["brdphi"] + bhat_phi * out["bphidphi"] + bhat_z * out["bzdphi"]) * bhat_phi / R.ravel()
        dbdl_z   = (bhat_r * out["brdz"]   + bhat_phi * out["bphidz"]   + bhat_z * out["bzdz"]) * bhat_z

        dbdl = np.sqrt(dbdl_r * dbdl_r + dbdl_phi * dbdl_phi + dbdl_z * dbdl_z).reshape(nin, nphi)

        print(np.nanmean(dbdl, axis=1).shape)
        phigrid = np.linspace(0, 2*np.pi, nphi+1)[:-1]
        avg = np.meshgrid(np.nanmean(dbdl, axis=1).ravel(), phigrid, indexing="ij")[0]
        return np.nanmax(np.abs(dbdl - avg), axis=1) / np.nanmean(dbdl, axis=1)


    def plotseparatrix(self, R, phi, z, t, axes):
        out = self.evaluate(R, phi, z, t, "rho", grid=True)

        mesh = axes.contour(R, z, np.transpose(out[:,0,:,0]), [1],
                            colors='black',zorder=1)


    def plotripple(self, Rgrid, zgrid, time, nphi, axes=None):
        R, z, t = np.meshgrid(Rgrid, zgrid, time, indexing="ij")
        rip = self.evaluateripple(R, z, t, nphi).reshape(Rgrid.size, zgrid.size)

        newfig = axes is None
        if newfig:
            plt.figure()
            axes = plt.gca()

        CS = axes.contour(Rgrid, zgrid, 100 * rip.transpose(),
                          [0.01, 0.05, 0.1, 0.5, 1, 5])
        axes.clabel(CS)

        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(Rgrid[0], Rgrid[-1])
        axes.set_ylim(zgrid[0], zgrid[-1])

        if newfig:
            plt.show(block=False)


    def plotripplewell(self, Rgrid, zgrid, time, nphi, axes=None):
        R, z, t = np.meshgrid(Rgrid, zgrid, time, indexing="ij")
        rip = self.evaluateripplewell(R, z, t, nphi).reshape(Rgrid.size, zgrid.size)

        newfig = axes is None
        if newfig:
            plt.figure()
            axes = plt.gca()

        CS = axes.contour(Rgrid, zgrid, 100 * rip.transpose(),
                          [0.1, 1, 10, 100, 200])
        axes.clabel(CS)

        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(Rgrid[0], Rgrid[-1])
        axes.set_ylim(zgrid[0], zgrid[-1])

        if newfig:
            plt.show(block=False)
