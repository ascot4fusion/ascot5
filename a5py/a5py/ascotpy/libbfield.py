"""
Methods to evaluate quantities from magnetic field data.

File: libbfield.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot

class LibBfield(LibAscot):

    quantities = ["rho", "psi", "br", "bphi", "bz", "brdr", "brdphi", "brdz",
                  "bphidr", "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz",
                  "divergence"]

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

        assert out is not None, "Unknown quantity"

        return out
