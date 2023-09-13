"""List of commonly encountered marker species and their properties.
"""
import numpy as np
import unyt

def species(name, charge=None):
    """Assumes the particle is fully ionized. Charge state can be given
    explicitly with charge parameter.
    """

    data = None
    if name in ["e", "electron"]:
        data = {"anum":0, "znum":0, "charge":-1*unyt.e, "mass":0.0005486*unyt.amu}

    if name in ["n", "neutron"]:
        data = {"anum":1, "znum":0, "charge":0*unyt.e,  "mass":1.009*unyt.amu}

    if name in ["p", "proton", "H", "H1"]:
        data = {"anum":1, "znum":1, "charge":1*unyt.e,  "mass":1.007*unyt.amu}

    if name in ["deuterium", "D", "H2"]:
        data = {"anum":2, "znum":1, "charge":1*unyt.e,  "mass":2.014*unyt.amu}

    if name in ["tritium", "T", "H3"]:
        data = {"anum":3, "znum":1, "charge":1*unyt.e,  "mass":3.016*unyt.amu}

    if name in ["He4", "alpha"]:
        data = {"anum":4, "znum":2, "charge":2*unyt.e,  "mass":4.003*unyt.amu}

    if charge is not None:
        data["charge"] = charge*unyt.e

    return data
