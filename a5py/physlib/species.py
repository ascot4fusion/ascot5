"""List of commonly encountered marker species and their properties.
"""
import numpy as np
from unyt import kg
import scipy.constants as constants

def species(name, charge=None):
    """

    Assumes the particle is fully ionized. Charge state can be given explicitly
    with charge parameter.
    """

    data = None
    if name in ["e", "electron"]:
        data = {"anum":0, "znum":0, "charge":-1, "mass":constants.electron_mass*kg}

    if name in ["n", "neutron"]:
        data = {"anum":1, "znum":0, "charge":0,  "mass":constants.neutron_mass*kg}

    if name in ["p", "proton", "H", "H1"]:
        data = {"anum":1, "znum":1, "charge":1,  "mass":constants.proton_mass*kg}

    if name in ["deuterium", "D", "H2"]:
        data = {"anum":2, "znum":1, "charge":1,  "mass":constants.physical_constants['deuteron mass'][0]*kg}

    if name in ["tritium", "T", "H3"]:
        data = {"anum":3, "znum":1, "charge":1,  "mass":constants.physical_constants['triton mass'][0]*kg}

    if name in ["helion","He3"]:
        data = {"anum":3, "znum":2, "charge":2,  "mass":constants.physical_constants['helion mass'][0]*kg}

    if name in ["He4", "alpha"]:
        data = {"anum":4, "znum":2, "charge":2,  "mass":constants.physical_constants['alpha particle mass'][0]*kg}

    if charge is not None:
        data["charge"] = charge

    return data


def test_species():
    for s in ['e','n','p','D','T','He3','alpha']:
        print(s,species(s))

def autodetect_species(anum,znum,charge=None):

    if anum == 0 and znum == 0:
        return species('e',charge)
    if anum == 1:
        if znum == 0:
            return species('n',charge)
        if znum == 1:
            return species('p',charge)
    if anum == 2 and znum == 1:
        return species('D', charge)
    if anum == 3:
        if znum == 1:
            return species('T',charge)
        if znum == 2:
            return species('He3',charge)
    if anum == 4 and znum == 2:
        return species('alpha',charge)

    raise(ValueError("Unknown species anum={} znum={}".format(anum,znum)))
