"""Fusion reaction data from AFSI.
"""
import numpy as np
import h5py
import unyt

class DataContainer():
    pass


class Reaction(DataContainer):
    """Class contining the stored fusion reaction data in AFSI run.
    """

    def read(self):
        """Read raw reaction data to a dictionary.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]

        return out

    def get(self):
        """Return reaction information.

        Returns
        -------
        m1 : float
            Mass of the first reactant.
        m2 : float
            Mass of the second reactant.
        mprod1 : float
            Mass of the first product.
        mprod2 : float
            Mass of the second product.
        q : float
            Energy released.
        """
        with self as f:
            m1     = f["m1"][:][0] * unyt.amu
            m2     = f["m2"][:][0] * unyt.amu
            mprod1 = f["mprod1"][:][0] * unyt.amu
            mprod2 = f["mprod2"][:][0] * unyt.amu
            q      = f["q"][:][0] * unyt.eV

        return m1, m2, mprod1, mprod2, q
