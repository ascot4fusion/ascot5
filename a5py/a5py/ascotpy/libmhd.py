"""
Methods to evaluate quantities from MHD data.

File: libmhd.py
"""
import numpy as np

from a5py.ascotpy.libascot import LibAscot


import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

class LibMhd(LibAscot):

    quantities = ["mhd_br", "mhd_bphi", "mhd_bz","mhd_er", "mhd_ephi", "mhd_ez"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["mhd_br", "mhd_bphi", "mhd_bz", "mhd_er", "mhd_ephi",
"mhd_ez"]:
            out = self.eval_mhd(R, phi, z, t)[quantity]

        assert out is not None, "Unknown quantity"

        return out
