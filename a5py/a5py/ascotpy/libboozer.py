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

    quantities = ["psi", "theta", "zeta", "dpsidr", "dpsidphi", "dpsidz", 
"dthetadr", "dthetadphi", "dthetadz", "dzetadr", "dzetadphi", "dzetadz"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None

        if quantity in ["psi", "theta", "zeta", "dpsidr", "dpsidphi", "dpsidz", 
"dthetadr", "dthetadphi", "dthetadz", "dzetadr", "dzetadphi", "dzetadz"]:
             out = self.eval_boozer(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
