"""List of commonly encountered marker species and their properties.
"""
import numpy as np
import unyt

def species(name, charge=None):
    """Retrieve anum, znum, mass, and charge based on the name of the species.

    Parameters
    ----------
    name : str
        Name of the species.
    charge : int, float
        Charge state or otherwise the ion is assumed to be fully ionized.

    Returns
    -------
    data : dict
        Contains "anum", "znum", "charge", and "mass".
    """
    valid_options = []
    def checkadd(options):
        valid_options.extend(options)
        return options

    data = None
    if name in checkadd(["e", "electron"]):
        data = {"anum":0, "znum":0, "charge":-1*unyt.e,
        "mass":0.0005486*unyt.amu}
    elif name in checkadd(["n", "neutron"]):
        data = {"anum":1, "znum":0, "charge":0*unyt.e, "mass":1.009*unyt.amu}
    elif name in checkadd(["p", "proton", "H", "H1"]):
        data = {"anum":1, "znum":1, "charge":1*unyt.e, "mass":1.007*unyt.amu}
    elif name in checkadd(["deuterium", "D", "H2"]):
        data = {"anum":2, "znum":1, "charge":1*unyt.e, "mass":2.014*unyt.amu}
    elif name in checkadd(["tritium", "T", "H3"]):
        data = {"anum":3, "znum":1, "charge":1*unyt.e, "mass":3.016*unyt.amu}
    elif name in checkadd(["He3"]):
        data = {"anum":3, "znum":2, "charge":2*unyt.e, "mass":3.016*unyt.amu}
    elif name in checkadd(["He4", "alpha"]):
        data = {"anum":4, "znum":2, "charge":2*unyt.e, "mass":4.003*unyt.amu}
    elif name in checkadd(["Be9"]):
        data = {"anum":9, "znum":4, "charge":4*unyt.e, "mass":9.012*unyt.amu}
    elif name in checkadd(["C12"]):
        data = {"anum":12, "znum":6, "charge":6*unyt.e, "mass":12.011*unyt.amu}
    elif name in checkadd(["Ne20"]):
        data = {"anum":20, "znum":10, "charge":10*unyt.e,
                "mass":19.992*unyt.amu}
    elif name in checkadd(["Ar40"]):
        data = {"anum":40, "znum":18, "charge":18*unyt.e,
                "mass":39.962*unyt.amu}
    elif name in checkadd(["Xe132"]):
        data = {"anum":132, "znum":54, "charge":54*unyt.e,
                "mass":131.904*unyt.amu}
    elif name in checkadd(["W184"]):
        data = {"anum":74, "znum":184, "charge":74*unyt.e,
                "mass":183.950*unyt.amu}
    else:
        species = ", ".join(valid_options)
        raise ValueError(
            "Unknown species %s. Known species are: %s" % (name, species))

    if charge is not None:
        data["charge"] = charge*unyt.e

    return data
