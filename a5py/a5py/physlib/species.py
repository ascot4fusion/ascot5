"""List of commonly encountered marker species and their properties.
"""
import numpy as np
from unyt import amu

def species(name, charge=None):
    """

    Assumes the particle is fully ionized. Charge state can be given explicitly
    with charge parameter.
    """

    data = None
    if name in ["e", "electron"]:
        data = {"anum":0, "znum":0, "charge":-1, "mass":0.0005486*amu}

    if name in ["n", "neutron"]:
        data = {"anum":1, "znum":0, "charge":0,  "mass":1.009*amu}

    if name in ["p", "proton", "H", "H1"]:
        data = {"anum":1, "znum":1, "charge":1,  "mass":1.007*amu}

    if name in ["deuterium", "D", "H2"]:
        data = {"anum":2, "znum":1, "charge":1,  "mass":2.014*amu}

    if name in ["tritium", "T", "H3"]:
        data = {"anum":3, "znum":1, "charge":1,  "mass":3.016*amu}

    if name in ["He4", "alpha"]:
        data = {"anum":4, "znum":2, "charge":2,  "mass":4.003*amu}

    if charge is not None:
        data["charge"] = charge

    return data
