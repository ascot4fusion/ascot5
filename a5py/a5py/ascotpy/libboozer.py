"""
Methods to evaluate quantities from boozer data.

File: libbbozer.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot


import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class LibBoozer(LibAscot):

    quantities = ["psi (bzr)", "theta", "zeta",
                  "dpsidr (bzr)", "dpsidphi (bzr)", "dpsidz (bzr)",
                  "dthetadr", "dthetadphi", "dthetadz",
                  "dzetadr", "dzetadphi", "dzetadz", "qprof", "jacobian",
                  "jacobianb2"]

    def evaluate(self, R, phi, z, t, quantity, br=None, bphi=None, bz=None):

        out = None
        if quantity == "psi (bzr)":
            out = self.eval_boozer(R, phi, z, t)["psi"]
        elif quantity == "dpsidr (bzr)":
            out = self.eval_boozer(R, phi, z, t)["dpsidr"]
        elif quantity == "dpsidphi (bzr)":
            out = self.eval_boozer(R, phi, z, t)["dpsidphi"]
        elif quantity == "dpsidz (bzr)":
            out = self.eval_boozer(R, phi, z, t)["dpsidz"]
        elif quantity in ["theta", "zeta",
                          "dthetadr", "dthetadphi", "dthetadz", "dzetadr",
                          "dzetadphi", "dzetadz"]:
             out = self.eval_boozer(R, phi, z, t)[quantity]
        elif quantity in ["qprof", "jacobian", "jacobianb2"]:
            out = self.eval_boozer(R, phi, z, t, evalfun=True)[quantity]

        assert out is not None, "Unknown quantity"

        return out
