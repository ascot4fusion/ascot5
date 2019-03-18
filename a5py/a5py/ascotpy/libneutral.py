"""
Methods to evaluate quantities from neutral data.

File: libneutral.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot


import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class LibNeutral(LibAscot):

    quantities = ["density"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["density"]:
            out = self.eval_neutral(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
