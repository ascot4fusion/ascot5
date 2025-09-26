"""ASCOT5 package for high-performance orbit-following simulations in fusion
plasma physics and engineering.
"""
import unyt

from a5py.composite import Ascot
from a5py.exceptions import AscotDataException


def set_up_units():
    """Set up unit system ASCOT5 so likes to use"""
    unyt.define_unit("markers", 1*unyt.Dimensionless)
    unyt.define_unit("particles", 1*unyt.Dimensionless)
    unyt.define_unit("e", -unyt.electron_charge)
    u = unyt.UnitSystem("ascot", "m", "atomic_mass_unit", "s", angle_unit="deg")
    u["energy"] = "eV"
    u["charge"] = "e"
    u["magnetic_field"] = "T"


set_up_units()


__all__ = [
    "Ascot",
    "AscotDataException",
    ]
