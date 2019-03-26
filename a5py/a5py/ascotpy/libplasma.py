"""
Methods to evaluate quantities from plasma data.

File: libplasma.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot


import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class LibPlasma(LibAscot):

    quantities = ["ne", "te"]

    def get_plasmaquantities(self):
        quantities = ["ne", "te"]

        spec = self.get_plasmaspecies()
        for i in range(1,spec["nspecies"]):
            quantities.append("ni" + str(i))
            quantities.append("ti" + str(i))

        return quantities


    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in self.get_plasmaquantities():
            out = self.eval_plasma(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
