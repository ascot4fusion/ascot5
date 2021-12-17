# coding: utf-8

# For Simppa at least, following recipe worked to get the below working:
# 0. compile ascotlib.so and add the path to LD_LIBRARY_PATH
# 1. in ascot5.yml, set the python version: python=3.6, and create the environment
# 2. module load imasenv
# 3. conda activate ascot5
# the below works.

import ctypes
import a5py.ascotpy.ascot5_main
import a5py.ascot5io.imas
import a5py.ascot5io.wall_2D
from   a5py.ascotpy import ascotpy2


# Read from IMAS

# 2D wall
w2d=a5py.ascot5io.imas.wall_2d()
wdict=w2d.read("g2jvarje","test","3",92436,272)
a5py.ascot5io.wall_2D.write_hdf5('from_imas.h5',**wdict)




M=a5py.ascotpy.ascot5_main.ascot5_main(input_filename=b'helloworld.h5',output_filename=b'helloworld_out.h5')
M.init()

# don't read wall input
what_to_read = ctypes.c_int32()
what_to_read.value =    ascotpy2.hdf5_input_options | ascotpy2.hdf5_input_bfield | \
                        ascotpy2.hdf5_input_efield  | ascotpy2.hdf5_input_plasma | \
                        ascotpy2.hdf5_input_neutral | ascotpy2.hdf5_input_marker |\
                        ascotpy2.hdf5_input_boozer  | ascotpy2.hdf5_input_mhd

#  ascotpy2.hdf5_input_wall |

M.read_input(what_to_read=what_to_read)

n_prts=3

# Reduce the simulation time a bit:
max_simtime=1.0e-4
print("Reducing max_simtime to {}s.".format(max_simtime))
M.sim.endcond_max_simtime = max_simtime

M.inject_wall_2d( wdict['r'], wdict['z']  )


M.offload()

M.run_simulation()

M.gather_output()

M.print_marker_summary()

M.write_output_h5()

M.free_c()

M.finalize()
