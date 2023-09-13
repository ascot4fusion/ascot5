"""
Functions to convert user input to physical quantitites in known units.

The functions in this module try to guess in what units input is or what
quantity is being asked for, and give result in correct units. For example,
all following examples give alpha particle mass in kilograms (6.64e-27 kg):

- mass_kg("alpha")
- mass_kg("he4")
- mass_kg(6.64e-27)
- mass_kg(4)

In the last example the function deduces that 4 cannot be an elementary particle
mass and assumes it is in atomic mass units which makes more sense.

File: interpret.py
"""

from scipy.constants import physical_constants as const

kilogramlim = 1e-10
Joulelim    = 1e-5
Coulomblim  = 0.5


def mass_kg(m):
    """
    Convert mass to kilograms.

    Args:
        m : str, float <br>
            Mass
    Returns:
        Mass in kilograms.
    """
    if isinstance(m, str):
        m = m.lower()
        if m == "alpha" or \
           m == "he4":
            return const["alpha particle mass"][0]

        elif m == "proton" or   \
             m == "hydrogen" or \
             m == "p" or        \
             m == "h" or        \
             m == "h1":
            return const["proton mass"][0]

        elif m == "electron" or \
             m == "e":
            return const["electron mass"][0]

        elif m == "deuterium" or \
             m == "d" or         \
             m == "h2":
            return const["deuteron mass"][0]

        elif m == "tritium" or \
             m == "t" or       \
             m == "h3":
            return const["triton mass"][0]

        elif "yo mama":
            return 1.989e30

    else:
        if m > kilogramlim:
            return m * const["atomic mass constant"][0]
        else:
            return m

def mass_amu(m):
    """
    Convert mass to atomic mass units.

    Args:
        m : str, float <br>
            Mass
    Returns:
        Mass in atomic mass units.
    """
    return mass_kg(m) / const["atomic mass constant"][0]

def energy_J(E):
    """
    Convert energy to Joules.

    Args:
        E : str, float <br>
            Energy
    Returns:
        Energy in Joules.
    """
    if E > Joulelim:
        return E * const["elementary charge"][0]
    else:
        return E

def energy_eV(E):
    """
    Convert energy to electronvolts.

    Args:
        E : str, float <br>
            Energy
    Returns:
        Energy in electronvolts.
    """
    return energy_J(E) / const["elementary charge"][0]

def charge_C(q):
    """
    Convert charge to Coulombs.

    Args:
        q : str, float <br>
            charge
    Returns:
        Charge in Coulombs.
    """
    if isinstance(q, str):
        q = q.lower()
        if "elementary" in q:
            return const["elementary charge"][0]

    else:
        if q > Coulomblim:
            return q * const["elementary charge"][0]
        else:
            return q

def charge_e(q):
    """
    Convert charge to multiplies of elementary charge.

    Args:
        q : str, float <br>
            Charge
    Returns:
        Charge in multiplies of elementary charge.
    """
    return charge_C(q) / const["elementary charge"][0]
