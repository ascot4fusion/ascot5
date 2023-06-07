"""
Methods to evaluate quantities from magnetic field data.

File: libbfield.py
"""
import numpy as np
from scipy.constants import physical_constants as const

import importlib.util as util

class LibBfield():

    quantities = ["rho", "psi", "br", "bphi", "bz", "brdr", "brdphi", "brdz",
                  "bphidr", "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz",
                  "divergence", "axis", "bnorm",
                  "ripplecritlarmor", "jnorm", "jr", "jphi", "jz"]


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

        bbar = (bhat_r * out["bphidr"] +  bhat_z * out["bphidz"]).reshape(nin, nphi)
        btil = (bhat_phi * out["bphidphi"] / R.ravel()).reshape(nin, nphi)

        return np.nanmean(np.abs(bbar), axis=1) \
            / np.nanmax(np.abs(btil), axis=1)


    


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


    def plotripplewell(self, Rgrid, zgrid, time, nphi, axes=None,
                       clevel=[-1,0,1], clabel=True, **kwargs):
        R, z, t = np.meshgrid(Rgrid, zgrid, time, indexing="ij")
        rip = self.evaluateripplewell(R, z, t, nphi).reshape(Rgrid.size, zgrid.size)

        newfig = axes is None
        if newfig:
            plt.figure()
            axes = plt.gca()

        CS = axes.contour(Rgrid, zgrid, np.log10(rip.transpose()),
                          clevel, **kwargs)
        if clabel:
            axes.clabel(CS)

        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(Rgrid[0], Rgrid[-1])
        axes.set_ylim(zgrid[0], zgrid[-1])

        if newfig:
            plt.show(block=False)
