"""Package for processing ASCOT5 data and generating inputs.
"""
import unyt

# Define the unit system ascot uses and add our own units and constants
try:
    unyt.define_unit("markers", 1*unyt.Dimensionless)
    unyt.define_unit("particles", 1*unyt.Dimensionless)
    unyt.define_unit("e", -unyt.electron_charge)
    u=unyt.UnitSystem("ascot", "m", "atomic_mass_unit", "s", angle_unit="deg")
    u["energy"] = "eV"
    u["charge"] = "e"
    u["magnetic_field"] = "T"
    unyt.mn = 1.675e-27*unyt.kg   # Neutron mass
    unyt.mD = 3.344e-27*unyt.kg   # Deuterium mass
    unyt.mT = 5.008e-27*unyt.kg   # Tritium mass
    unyt.mHe3 = 5.008e-27*unyt.kg # Helium-3 mass
    unyt.mHe4 = 6.646e-27*unyt.kg # Helium-4 mass
except RuntimeError:
    # We get exception when trying to define unit that is already defined.
    # This can be ignored.
    pass


from a5py.composite import Ascot
from a5py.exceptions import AscotIOException

__all__ = [
    "Ascot",
    "AscotIOException",
    ]
