import numpy as np
import unyt
#import matplotlib.pyplot as plt
from a5py import Ascot

filename='ascot_vrun.h5'

import os


# Check if the file exists before attempting to delete it
if os.path.exists(filename):
    os.remove(filename)
    print(f"The file {filename} has been deleted.")
else:
    print(f"The file {filename} does not exist.")



a5 = Ascot(filename, create=True)
a5.data.create_input("bfield analytical iter circular")
a5.data.create_input("plasma flat", density=1e21)
a5.data.create_input("wall_2D")
a5.data.create_input("E_TC")
a5.data.create_input("N0_1D")
a5.data.create_input("Boozer")
a5.data.create_input("MHD_STAT")
a5.data.create_input("asigma_loc")

from a5py.ascot5io.marker import Marker
nmrk = 1000
mrk = Marker.generate("gc", n=nmrk, species="alpha")
mrk["energy"][:] = 3.5e6
mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(nmrk,)
mrk["r"][:]      = 4.5 + 3 * np.random.rand(nmrk,)
a5.data.create_input("gc", **mrk)

from a5py.ascot5io.options import Opt
opt = Opt.get_default()
opt.update({
    "SIM_MODE":2, "ENABLE_ADAPTIVE":1,
    "ENDCOND_ENERGYLIM":1, "ENDCOND_MIN_ENERGY":2.0e3, "ENDCOND_MIN_THERMAL":2.0,
    "ENABLE_ORBIT_FOLLOWING":1, "ENABLE_COULOMB_COLLISIONS":1,
})

print("Inputs created")



opt.update({
    # Distribution output
    "ENABLE_DIST_5D":1, "ENABLE_DIST_RHO5D":1,
    # (R,z) abscissae for the 5D distribution
    "DIST_MIN_R":4.3,  "DIST_MAX_R":8.3, "DIST_NBIN_R":50,
    "DIST_MIN_Z":-2.0, "DIST_MAX_Z":2.0, "DIST_NBIN_Z":50,
    # (rho, theta) abscissae for the rho5D distribution. Most of the time a single
    # theta slot is sufficient but please verify it in your case.
    "DIST_MIN_RHO"  :0, "DIST_MAX_RHO"  :1.0, "DIST_NBIN_RHO"  :100,
    "DIST_MIN_THETA":0, "DIST_MAX_THETA":360, "DIST_NBIN_THETA":1,
    # Single phi slot since this is not a stellarator.
    # These values are shared between other distributions
    "DIST_MIN_PHI":0,        "DIST_MAX_PHI":360,     "DIST_NBIN_PHI":1,
    # The momentum abscissae are shared by 5D distributions
    "DIST_MIN_PPA":-1.3e-19, "DIST_MAX_PPA":1.3e-19, "DIST_NBIN_PPA":100,
    "DIST_MIN_PPE":0,        "DIST_MAX_PPE":1.3e-19, "DIST_NBIN_PPE":50,
    # One time slot, the span doesn't matter as long as it covers the whole simulation time
    "DIST_MIN_TIME":0,       "DIST_MAX_TIME":1.0,    "DIST_NBIN_TIME":1,
    # One charge slot exactly at q=2 since we are simulating alphas
    "DIST_MIN_CHARGE":1,     "DIST_MAX_CHARGE":3,    "DIST_NBIN_CHARGE":1,
})
a5.data.create_input("opt", **opt)

a5 = Ascot(filename)

# You may use any options and markers but in this example we take them from
# the input file.
mrk = a5.data.marker.active.read()
opt = a5.data.options.active.read()

a5.simulation_initinputs()
a5.simulation_initmarkers(**mrk)
a5.simulation_initoptions(**opt)

print('running simulation')
# Running simulations returns a virtual run that acts the same way as normal runs
vrun = a5.simulation_run()


# Get us a distribution object for testing
import a5py.ascot5io.imas
dist = a5py.ascot5io.imas.distributions()


dist.write(user='akaslos', tokamak='test', version='3', shot=1, run=1, runobject=vrun, metadata={} )



print(''' Try the following:

vrun.getdist_list()
d=vrun.getdist('5d')
vrun.getdist_moments?
dist.write(user='akaslos', tokamak='test', version='3', shot=1, run=1, runobject=vrun, metadata={} )
''')

