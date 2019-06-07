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
                  "divergence", "axis", "bnorm"]

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

        assert out is not None, "Unknown quantity"

        return out


    def plotseparatrix(self, R, phi, z, t, axes):
        out = self.evaluate(R, phi, z, t, "rho", grid=True)

        mesh = axes.contour(R, z, np.transpose(out[:,0,:,0]), [1], colors='black',zorder=1)
