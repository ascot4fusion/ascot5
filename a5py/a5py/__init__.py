# Define the unit system ascot uses and add our own unit types
import unyt

unyt.define_unit("markers", 1*unyt.Dimensionless)
unyt.define_unit("e", -unyt.electron_charge)

u=unyt.UnitSystem("ascot", "m", "atomic_mass_unit", "s", angle_unit="deg")
u["energy"] = "eV"
u["charge"] = "e"
u["magnetic_field"] = "T"
