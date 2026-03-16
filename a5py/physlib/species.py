"""List of commonly encountered marker species and their properties.
"""
import unyt


speciesdict = {
    "e"     : (  0,  0, -1*unyt.e,   0.0005486*unyt.amu),
    "n"     : (  1,  0,  0*unyt.e,   1.009*unyt.amu),
    "H"     : (  1,  1,  1*unyt.e,   1.007*unyt.amu),
    "D"     : (  2,  1,  1*unyt.e,   2.014*unyt.amu),
    "T"     : (  3,  1,  1*unyt.e,   3.016*unyt.amu),
    "He3"   : (  3,  2,  2*unyt.e,   3.016*unyt.amu),
    "He4"   : (  4,  2,  2*unyt.e,   4.003*unyt.amu),
    "Li6"   : (  6,  3,  3*unyt.e,   6.015*unyt.amu),
    "Li7"   : (  7,  3,  3*unyt.e,   7.016*unyt.amu),
    "Be9"   : (  9,  4,  4*unyt.e,   9.012*unyt.amu),
    "B10"   : ( 10,  5,  5*unyt.e,  10.012*unyt.amu),
    "B11"   : ( 11,  5,  5*unyt.e,  11.009*unyt.amu),
    "C12"   : ( 12,  6,  6*unyt.e,  12.011*unyt.amu),
    "C13"   : ( 13,  6,  6*unyt.e,  13.003*unyt.amu),
    "N14"   : ( 14,  7,  7*unyt.e,  14.003*unyt.amu),
    "N15"   : ( 15,  7,  7*unyt.e,  15.000*unyt.amu),
    "O16"   : ( 16,  8,  8*unyt.e,  15.994*unyt.amu),
    "O17"   : ( 17,  8,  8*unyt.e,  16.999*unyt.amu),
    "O18"   : ( 18,  8,  8*unyt.e,  17.999*unyt.amu),
    "Ne20"  : ( 20, 10, 10*unyt.e,  19.992*unyt.amu),
    "Ne22"  : ( 22, 10, 10*unyt.e,  21.991*unyt.amu),
    "Ar36"  : ( 36, 18, 18*unyt.e,  35.967*unyt.amu),
    "Ar40"  : ( 40, 18, 18*unyt.e,  39.962*unyt.amu),
    "Ni59"  : ( 59, 28, 28*unyt.e,  58.934*unyt.amu),
    "Sn118" : (118, 50, 50*unyt.e, 117.901*unyt.amu),
    "Xe132" : (132, 54, 54*unyt.e, 131.904*unyt.amu),
    "W183"  : (183, 74, 74*unyt.e, 182.950*unyt.amu),
    "W184"  : (184, 74, 74*unyt.e, 183.950*unyt.amu),
}

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
        data = speciesdict["e"]
    elif name in checkadd(["n", "neutron"]):
        data = speciesdict["n"]
    elif name in checkadd(["H", "p", "proton", "H1"]):
        data = speciesdict["H"]
    elif name in checkadd(["D", "deuterium", "H2"]):
        data = speciesdict["D"]
    elif name in checkadd(["T", "tritium", "H3"]):
        data = speciesdict["T"]
    elif name in checkadd(["He3", "helion"]):
        data = speciesdict["He3"]
    elif name in checkadd(["He4", "alpha"]):
        data = speciesdict["He4"]
    elif name in checkadd(["Be9"]):
        data = speciesdict["Be9"]
    elif name in checkadd(["C12"]):
        data = speciesdict["C12"]
    elif name in checkadd(["Ne20"]):
        data = speciesdict["Ne20"]
    elif name in checkadd(["Ar40"]):
        data = speciesdict["Ar40"]
    elif name in checkadd(["Ni59"]):
        data = speciesdict["Ni59"]
    elif name in checkadd(["Xe132"]):
        data = speciesdict["Xe132"]
    elif name in checkadd(["W184"]):
        data = speciesdict["W184"]
    else:
        spec = ", ".join(valid_options)
        raise ValueError(
            f"Unknown species {name}. Known species are: {spec}"
        )

    data = {"anum":data[0], "znum":data[1], "charge":data[2], "mass":data[3]}
    if charge is not None:
        data["charge"] = charge

    return data

def autodetect(anum, znum, charge=None):
    """Return species based on given anum and znum.

    Parameters
    ----------
    anum : int
        Atomic mass number.
    znum : int
        Charge number.
    charge : int, optional
        Charge state of the returned species or fully ionized if None.

    Returns
    -------
    data : dict
        Contains "anum", "znum", "charge", "mass", and "name".
    """
    for key,v in speciesdict.items():
        if anum == v[0] and znum == v[1]:
            if charge is None:
                charge = v[2]
            return {"anum":v[0], "znum":v[1], "charge":charge, "mass":v[3], "name": key}

    raise ValueError(f"Unknown species anum={anum} znum={znum}")
