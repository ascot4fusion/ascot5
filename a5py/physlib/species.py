"""List of commonly encountered marker species and their properties.
"""
from typing import NamedTuple

from unyt import e, amu


class Species(NamedTuple):
    anum: int
    znum: int
    charge: float
    mass: float
    """Named tuple containing species properties.

    Attributes
    ----------
    anum : int
        Atomic mass number.
    znum : int
        Charge number.
    charge : unyt.quantity
        Charge state of the species.
    mass : unyt.quantity
        Mass of the species.
    """


KNOWN_SPECIES = {
    "e"     : Species(  0,  0, -1*e,   0.0005486*amu),
    "n"     : Species(  1,  0,  0*e,   1.009*amu),
    "H1"     : Species(  1,  1,  1*e,   1.007*amu),
    "H2"     : Species(  2,  1,  1*e,   2.014*amu),
    "H3"     : Species(  3,  1,  1*e,   3.016*amu),
    "He3"   : Species(  3,  2,  2*e,   3.016*amu),
    "He4"   : Species(  4,  2,  2*e,   4.003*amu),
    "Be9"   : Species(  9,  4,  4*e,   9.012*amu),
    "C12"   : Species( 12,  6,  6*e,  12.011*amu),
    "Ne20"  : Species( 20, 10, 10*e,  19.992*amu),
    "Ar40"  : Species( 40, 18, 18*e,  39.962*amu),
    "Ni59"  : Species( 59, 28, 28*e,  58.934*amu),
    "Xe132" : Species(132, 54, 54*e, 131.904*amu),
    "W183"  : Species(183, 74, 74*e, 182.950*amu),
    "W184"  : Species(184, 74, 74*e, 183.950*amu),
}
"""Names of the recognized species and their properties."""


def species2properties(name, charge=None):
    """Retrieve nuclear properties (atomic mass number, charge number, mass,
    and charge) based on the name of the particle species.

    For electrons/positrons (anum=0, znum=0) and for neutrons (anum=1, znum=0).

    Parameters
    ----------
    name : str
        Name of the species.
    charge : int, optional
        Charge state.

        If not provided, the ion is assumed to be fully ionized.

    Returns
    -------
    species : Species
        Named tuple containing species properties.
    """
    valid_options = []
    def checkadd(options):
        valid_options.extend(options)
        return options

    species = None
    if name in checkadd(["e", "electron"]):
        species = KNOWN_SPECIES["e"]
    elif name in checkadd(["n", "neutron"]):
        species = KNOWN_SPECIES["n"]
    elif name in checkadd(["H1", "H", "p", "proton"]):
        species = KNOWN_SPECIES["H1"]
    elif name in checkadd(["H2", "D", "deuterium"]):
        species = KNOWN_SPECIES["H2"]
    elif name in checkadd(["H3", "T", "tritium"]):
        species = KNOWN_SPECIES["H3"]
    elif name in checkadd(["He3", "helion"]):
        species = KNOWN_SPECIES["He3"]
    elif name in checkadd(["He4", "alpha"]):
        species = KNOWN_SPECIES["He4"]
    elif name in checkadd(["Be9"]):
        species = KNOWN_SPECIES["Be9"]
    elif name in checkadd(["C12"]):
        species = KNOWN_SPECIES["C12"]
    elif name in checkadd(["Ne20"]):
        species = KNOWN_SPECIES["Ne20"]
    elif name in checkadd(["Ar40"]):
        species = KNOWN_SPECIES["Ar40"]
    elif name in checkadd(["Ni59"]):
        species = KNOWN_SPECIES["Ni59"]
    elif name in checkadd(["Xe132"]):
        species = KNOWN_SPECIES["Xe132"]
    elif name in checkadd(["W184"]):
        species = KNOWN_SPECIES["W184"]
    else:
        all_species_names = ", ".join(valid_options)
        raise ValueError(
            f"Unknown species {name}. Known species are: {all_species_names}"
        )

    if charge is not None:
        if charge != int(charge):
            raise ValueError("Charge must be an integer.")
        species = Species(species.anum, species.znum, charge, species.mass)
    return species


def properties2species(anum ,znum):
    """Identify the species and return its name based on the atomic mass number
    and charge number.

    Parameters
    ----------
    anum : int
        Atomic mass number.
    znum : int
        Charge number.

    Returns
    -------
    name : str
        Name of the species.
    """
    for name, properties in KNOWN_SPECIES.items():
        if anum == properties.anum and znum == properties.znum:
            return name
    raise ValueError(
        f"Unknown species anum={anum} znum={znum}. See "
        "a5py.physlib.KNOWN_SPECIES for a list of known species."
        )


def findmass(anum, znum):
    """Identify the species and return its mass based on the atomic mass number
    and charge number.

    Parameters
    ----------
    anum : int
        Atomic mass number.
    znum : int
        Charge number.

    Returns
    -------
    mass : float
        Mass of the species.
    """
    for properties in KNOWN_SPECIES.values():
        if anum == properties.anum and znum == properties.znum:
            return properties.mass
    raise ValueError(
        f"Unknown species anum={anum} znum={znum}. See "
        "a5py.physlib.KNOWN_SPECIES for a list of known species."
        )
