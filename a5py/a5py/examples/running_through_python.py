# coding: utf-8

import ctypes
from   a5py.ascotpy import ascotpy2
import a5py.ascotpy.ascot5_main
import numpy as np


def create_wall_2d():

    r=np.array([4.96, 5.82, 6.68, 7.54,
                8.4, 8.4, 8.4, 8.4, 8.4, 8.4,
                7.54, 6.68, 5.82, 4.96,
                4.1, 4.1, 4.1, 4.1, 4.1, 4.1 ] )

    z=np.array([ -3.9, -3.9,  -3.9,  -3.9,
                 -3.9, -2.34, -0.78, 0.78, 2.34, 3.9,
                 3.9,   3.9,   3.9,  3.9,
                 3.9,   2.34,  0.78, -0.78, -2.34, -3.9])

    return r,z

def create_particles(n_prts=3):


    prt_array = ascotpy2.libascot_allocate_input_particles(n_prts)
    #prt_array = (ascotpy2.struct_c__SA_input_particle * n_prts )() #This cannot be freed on the C-side
    
    prt_type = {v: k for k, v in ascotpy2.input_particle_type__enumvalues.items()}['input_particle_type_p']

    # p_phi is a function of i    
    for i in range(n_prts):
        #print('creating marker {}'.format(i) )
        prt_array[i].type = prt_type
       
        prt_array[i].c__SA_input_particle_0.p.id     = i+1

        prt_array[i].c__SA_input_particle_0.p.weight = 1.0

       
        prt_array[i].c__SA_input_particle_0.p.mileage  = 0.0
        prt_array[i].c__SA_input_particle_0.p.time     = 0.0
        
        prt_array[i].c__SA_input_particle_0.p.anum   = 4
        prt_array[i].c__SA_input_particle_0.p.mass   = 6.642155684e-27
        prt_array[i].c__SA_input_particle_0.p.znum   = 2
        prt_array[i].c__SA_input_particle_0.p.charge = 3.20435313e-19
        
        prt_array[i].c__SA_input_particle_0.p.r      =  4.5
        prt_array[i].c__SA_input_particle_0.p.phi    =  0.0
        prt_array[i].c__SA_input_particle_0.p.z      =  0.0

        prt_array[i].c__SA_input_particle_0.p.p_r    =  7.0e-20
        prt_array[i].c__SA_input_particle_0.p.p_phi  =  (1.0+i) * 1.0e-20
        prt_array[i].c__SA_input_particle_0.p.p_z    =  0.0
    
    #prt_point = ctypes.cast(prt_array, ctypes.POINTER(ascotpy2.struct_c__SA_input_particle))
        
    return (prt_array,n_prts)

M=a5py.ascotpy.ascot5_main.ascot5_main(input_filename=b'helloworld.h5',output_filename=b'helloworld_out.h5')
M.init()

# don't read marker input
what_to_read = ctypes.c_int32()
what_to_read.value =    ascotpy2.hdf5_input_options | ascotpy2.hdf5_input_bfield | \
                        ascotpy2.hdf5_input_efield  | ascotpy2.hdf5_input_plasma | \
                        ascotpy2.hdf5_input_neutral |\
                        ascotpy2.hdf5_input_boozer  | ascotpy2.hdf5_input_mhd

#  ascotpy2.hdf5_input_wall |

M.read_input(what_to_read=what_to_read)

n_prts=3

# Reduce the simulation time a bit:
max_simtime=1.0e-4
print("Reducing max_simtime to {}s.".format(max_simtime))
M.sim.endcond_max_simtime = max_simtime

# Generate a few markers of our own.
(prt_array,n_tot) = create_particles(n_prts=n_prts)
M.n_tot = n_tot
M.prt.contents = prt_array.contents
#M.prt.contents = prt_array

#print("id=",M.prt[n_prts-1].c__SA_input_particle_0.p.id)


R,z = create_wall_2d()
M.inject_wall_2d( R, z  )


M.offload()

M.run_simulation()

M.gather_output()

print('marker summary:')
M.print_marker_summary()

M.write_output_h5()

M.free_c()

M.finalize()
