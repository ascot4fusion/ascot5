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

    quantities = ["mhd_br", "mhd_bphi", "mhd_bz","mhd_er", "mhd_ephi",
                  "mhd_ez","mhd_phi", "db/b"]

    def evaluate(self, R, phi, z, t, quantity):

        out = None
        if quantity in ["mhd_br", "mhd_bphi", "mhd_bz", "mhd_er", "mhd_ephi",
                        "mhd_ez", "mhd_phi"]:
            out = self.eval_mhd(R, phi, z, t)[quantity]

        if quantity in ["db/b"]:
            bpert = self.eval_mhd(R, phi, z, t)
            b = self.eval_bfield(R, phi, z, t, evalb=True)

            bpert = np.sqrt(bpert["mhd_br"]**2 + bpert["mhd_bphi"]**2
                            +  bpert["mhd_bz"]**2)
            b = np.sqrt(b["br"]**2 + b["bphi"]**2 + b["bz"]**2)

            out = bpert/b

        assert out is not None, "Unknown quantity"

        return out
