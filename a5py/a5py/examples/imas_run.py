# coding: utf-8

# For Simppa at least, following recipe worked to get the below working:
# 0. compile ascotlib.so and add the path to LD_LIBRARY_PATH
# 1. in ascot5.yml, set the python version: python=3.6, and create the environment
# 2. module load imasenv
# 3. conda activate ascot5
# the below works.

import ctypes
from   a5py import Ascot
import a5py.ascot5io.imas
import a5py.ascot5io.wall
from   a5py.ascotpy import ascot2py
from   a5py.ascot5io.marker import Marker
from a5py.templates.optiontemplates import OptionTemplates
import numpy as np

# Read from IMAS

print('Reading 2D wall')
# 2D wall
w2d=a5py.ascot5io.imas.wall_2d()
#wdict2=w2d.read("g2jvarje","test","3",92436,272)
wdict2=w2d.read("akaslos","test","3",92436,272)
a5py.ascot5io.wall_2D.write_hdf5('from_imas.h5',**wdict2)

print('Reading 3D wall')
# 3D wall
w3d=a5py.ascot5io.imas.wall_3d()
wdict3=w3d.read("akaslos","mywall","3",201,101)
#print(wdict3['x1x2x3'].shape, wdict3['y1y2y3'].shape, wdict3['z1z2z3'].shape, )
a5py.ascot5io.wall_3D.write_hdf5('from_imas.h5',**wdict3)

print('Reading Bfield (stellarator)')
# 3D Bfield in Stellarator format
bsts=a5py.ascot5io.imas.B_STS()
bdict=bsts.read("akaslos","ggdtest","3",32,3)
a5py.ascot5io.B_STS.write_hdf5('from_imas.h5',**bdict)


print('Initializing ascot5')


a5 = Ascot('helloworld.h5')
#M=a5py.ascotpy.ascot5_main.ascot5_main(input_filename=b'helloworld.h5',output_filename=b'helloworld_out.h5')



#print("Inializing.")
#M.init()

# don't read wall input
#what_to_read = ctypes.c_int32()
#what_to_read.value =    ascotpy2.hdf5_input_options | \
#                        ascotpy2.hdf5_input_efield  | ascotpy2.hdf5_input_plasma | \
#                        ascotpy2.hdf5_input_neutral | ascotpy2.hdf5_input_marker |\
#                        ascotpy2.hdf5_input_boozer  | ascotpy2.hdf5_input_mhd

#  ascotpy2.hdf5_input_wall |
# ascotpy2.hdf5_input_bfield | \

#M.read_input(what_to_read=what_to_read)
print('Initializing inputs')

#a5.simulation_initinputs(bfield=True, efield=True, plasma=True, neutral=True,
#                        wall=True, boozer=True, mhd=True, switch=True)

bfield=bdict # True
wall=wdict3  # wdict2, True
a5.simulation_initinputs(bfield=bdict, efield=True, plasma=True, neutral=True,
                         wall=wdict3, boozer=True, mhd=True, switch=True)


print('generating markers')
nmarkers = 10
mrk = Marker.generate("gc", n=nmarkers, species="alpha")
mrk["energy"][:] = 3.5e6
mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(nmarkers,)
mrk["r"][:]      = np.linspace(3.5, 5.4, nmarkers)
a5.data.create_input("gc", **mrk)

#print('Reading markers')
#mrk = a5.data.marker.active.read()
#print(mrk.keys())
print('providing markers')
a5.simulation_initmarkers(**mrk)


print('generating options')
#a5.data.create_input("options tutorial")
opt = OptionTemplates().options_tutorial()

#print('Reading options')
#opt = a5.data.options.active.read()
print('providing options')
a5.simulation_initoptions(**opt[1])



# Reduce the simulation time a bit:
max_simtime=1.0e-4
print("Reducing max_simtime to {}s.".format(max_simtime))
#M.sim.endcond_max_simtime = max_simtime
a5._sim.endcond_lim_simtime=max_simtime


#print("Injecting 2D wall from imas")
#M.inject_wall_2d( wdict['r'], wdict['z']  )

#print("Injecting 3D wall from imas")
#a5.provide_wall_3d( wdict3['x1x2x3'], wdict3['y1y2y3'], wdict3['z1z2z3']  )

#print("Injecting stellarator B-field from imas")
#a5.provide_BSTS(bdict)


print('starting simulation')
vrun_output = a5.simulation_run()





#print("offloading")
#M.offload()

#print("simulating")
#M.run_simulation()

#print("gathering output")
#M.gather_output()

#print("printing marker summary")
#M.print_marker_summary()

#print("writing output h5")
#M.write_output_h5()

print("freeing memory")
#M.free_c()
a5.simulation_free(diagnostics=True)

#print("finalizing MPI")
#M.finalize()

print('Finished test')


#print('Stopping test')
#quit()
