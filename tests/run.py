import ctypes
import unyt
import numpy as np
from a5py import Ascot
from a5py.templates import PremadeMagneticField
from a5py.engine.simulate import setup_run, setup_markers, simulate

from a5py.data.bfield.analytical import Struct

print(f"real {ctypes.sizeof(Struct)}\n")

a5 = Ascot()
template = PremadeMagneticField(a5, field="iter-baseline")
template.create_input()
#a5.data.create_bfieldcartesian(
#    bxyz=np.array([1., 1., 1.]), jacobian=np.full((3,3), 0.0), axisrz=(6.2, 0.0), rhoval=1.0,
#    )

opt = a5.data.create_options(
    enable_orbit_following=True,
    activate_simulation_time_limits=True,
    simulation_mode=4,
    activate_real_time_limit=True,
    max_real_time=1.0,
    activate=True,
)

mrk = a5.data.create_fieldlinemarker(
    r=6.*unyt.m,
    z=0.*unyt.m,
    phi=0.0,
)

#a5.data.bfield.active.stage()
run = setup_run(opt, options=a5.data.options.active, bfield=a5.data.bfield.active)

setup_markers(mrk, a5.data.bfield.active)
run._setup(mrk._cdata)

simulate(run)

print(run._diagnostics["endstate"].r)
#a5.init(0)
