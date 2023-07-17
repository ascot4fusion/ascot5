# coding: utf-8

import ctypes
from   a5py.ascotpy import ascotpy2
import a5py.ascotpy.ascot5_main
import numpy as np

def create_particles(R, nprt):

    prt_array = ascotpy2.libascot_allocate_input_particles(nprt)
    particletypes = ascotpy2.input_particle_type__enumvalues.items()
    prt_type = {v: k for k, v in particletypes}['input_particle_type_p']
    for i in range(nprt):
        prt_array[i].type = prt_type

        prt_array[i].c__SA_input_particle_0.p.id     = i+1

        prt_array[i].c__SA_input_particle_0.p.weight = 1.0

        prt_array[i].c__SA_input_particle_0.p.mileage  = 0.0
        prt_array[i].c__SA_input_particle_0.p.time     = 0.0

        prt_array[i].c__SA_input_particle_0.p.anum   = 4
        prt_array[i].c__SA_input_particle_0.p.mass   = 6.642155684e-27
        prt_array[i].c__SA_input_particle_0.p.znum   = 2
        prt_array[i].c__SA_input_particle_0.p.charge = 3.20435313e-19

        prt_array[i].c__SA_input_particle_0.p.r      =  R
        prt_array[i].c__SA_input_particle_0.p.phi    =  0.0
        prt_array[i].c__SA_input_particle_0.p.z      =  0.0

        prt_array[i].c__SA_input_particle_0.p.p_r    =  7.0e-20
        prt_array[i].c__SA_input_particle_0.p.p_phi  =  (1.0) * 1.0e-20
        prt_array[i].c__SA_input_particle_0.p.p_z    =  0.0

    return (prt_array, nprt)

def initandrun(fn, rcoord, init=True):

    M = a5py.ascotpy.ascot5_main.ascot5_main(input_filename=fn, \
                                             output_filename=fn)
    M.init(init)

    # Skip marker input
    inputs = ctypes.c_int32()
    inputs = ascotpy2.hdf5_input_options | ascotpy2.hdf5_input_bfield |\
             ascotpy2.hdf5_input_efield  | ascotpy2.hdf5_input_plasma |\
             ascotpy2.hdf5_input_neutral | ascotpy2.hdf5_input_wall   |\
             ascotpy2.hdf5_input_boozer  | ascotpy2.hdf5_input_mhd

    M.read_input(what_to_read=inputs)

    nprt = 1

    # Set simulation time
    max_simtime = 1.0e-4
    M.sim.endcond_max_simtime = max_simtime

    (prt_array,n_tot) = create_particles(rcoord, nprt)
    M.n_tot = n_tot
    M.prt.contents = prt_array.contents

    M.offload()

    M.run_simulation()

    M.gather_output()

    return M



#M.write_output_h5()
def clean(M):
    M.free_c()
    #M.finalize()
