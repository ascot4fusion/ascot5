"""
Methods to evaluate quantities from electric field data.

File: libefield.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot


import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class LibEfield(LibAscot):

    quantities = ["er", "ephi", "ez"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["er", "ephi", "ez"]:
            out = self.eval_efield(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
