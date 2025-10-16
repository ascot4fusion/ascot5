import numpy as np
import unyt
import matplotlib.pyplot as plt
from a5py import Ascot

inp = {}
dryrun=True

a5 = Ascot( create=False)
inp['bfield']  = a5.data.create_input("bfield analytical iter circular", dryrun=dryrun)
inp['plasma']  = a5.data.create_input("plasma flat", density=1e21,       dryrun=dryrun)
inp['wall']    = a5.data.create_input("wall_2D",                         dryrun=dryrun)
inp['efield']  = a5.data.create_input("E_TC",                            dryrun=dryrun)
inp['neutral'] = a5.data.create_input("N0_1D",                           dryrun=dryrun)
inp['boozer']  = a5.data.create_input("Boozer",                          dryrun=dryrun)
inp['mhd']     = a5.data.create_input("MHD_STAT",                        dryrun=dryrun)
inp['asigma']  = a5.data.create_input("asigma_loc",                      dryrun=dryrun)

from a5py.ascot5io.marker import Marker
nmrk = 1000
mrk = Marker.generate("gc", n=nmrk, species="alpha")
mrk["energy"][:] = 3.5e6
mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(nmrk,)
mrk["r"][:]      = 4.5 + 3 * np.random.rand(nmrk,)
#a5.data.create_input("gc", **mrk)

from a5py.ascot5io.options import Opt
opt = Opt.get_default()
opt.update({
    "SIM_MODE":2, "ENABLE_ADAPTIVE":1,
    "ENDCOND_ENERGYLIM":1, "ENDCOND_MIN_ENERGY":2.0e3, "ENDCOND_MIN_THERMAL":2.0,
    "ENABLE_ORBIT_FOLLOWING":1, "ENABLE_COULOMB_COLLISIONS":1,
})

print("Inputs created:")

#print(inp)

print('Updating options')

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
# a5.data.create_input("opt", **opt)


print("Setting up distributions in inputs done")

print("init inputs")

a5.simulation_initinputs(**inp)

print("init markers")

a5.simulation_initmarkers(**mrk)


print("init options")
a5.simulation_initoptions(**opt)

print("Starting simulations")

vrun = a5.simulation_run()

print("Simulation completed")

vrun.getdist_list()
d=vrun.getdist('5d')

pickle_file = 'dist5d.pickle'

print("pickling distribution to '"+pickle_file+"'")

import pickle

with open(pickle_file, 'wb') as f:
   pickle.dump(d,f,pickle.HIGHEST_PROTOCOL)

print('done')


# to access the pickle
#import pickle
#with open('dist5d.pickle', 'rb') as f:
#   d=pickle.load(f)
