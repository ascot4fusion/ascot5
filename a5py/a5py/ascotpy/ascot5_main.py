'''
Created on Nov 16, 2021

@author: sjjamsa
'''

import ctypes

# generate ascotpy2 in the source folder with
# % clang2py -l ./libascot.so -o ../python/a5py/a5py/ascotpy/ascotpy2.py  particle.h hdf5_interface.h ascot5.h mpi_interface.h --clang-args="-I/usr/include/hdf5/serial"
from a5py.ascotpy import ascotpy2



class ascot5_main(object):
    '''
    classdocs
    '''

    sim = ascotpy2.struct_c__SA_sim_offload_data()

    
    mpi_rank = ctypes.c_int32()
    mpi_size = ctypes.c_int32()
    mpi_root = ctypes.c_int32()

    # Pointers to hold temporary arrays
    B_offload_array       = ctypes.POINTER(ctypes.c_double)()
    E_offload_array       = ctypes.POINTER(ctypes.c_double)()
    plasma_offload_array  = ctypes.POINTER(ctypes.c_double)()
    neutral_offload_array = ctypes.POINTER(ctypes.c_double)()
    wall_offload_array    = ctypes.POINTER(ctypes.c_double)()
    boozer_offload_array  = ctypes.POINTER(ctypes.c_double)()
    mhd_offload_array     = ctypes.POINTER(ctypes.c_double)()
    prt                   = ctypes.POINTER(ascotpy2.struct_c__SA_input_particle)()



    def __init__(self, input_filename=b'input.h5'):
        '''
        Constructor
        '''
        
        self.sim.hdf5_in=input_filename
        


    def init(self):

        # No argumets given at the moment.
        argc=ctypes.c_int32()
        argc.value=0
        argv=ctypes.POINTER(ctypes.c_char)()
        argv.contents=None
        
        ascotpy2.mpi_interface_init(argc,
                                    ctypes.byref(argv),
                                    ctypes.byref(self.sim), 
                                    ctypes.byref(self.mpi_rank),
                                    ctypes.byref(self.mpi_size),
                                    ctypes.byref(self.mpi_root)  ) 
        
    def read_input(self):
        
         
        what_to_read = ctypes.c_int32()
        what_to_read.value =    ascotpy2.hdf5_input_options | ascotpy2.hdf5_input_bfield | \
                                ascotpy2.hdf5_input_efield  | ascotpy2.hdf5_input_plasma | \
                                ascotpy2.hdf5_input_neutral | ascotpy2.hdf5_input_wall   | \
                                ascotpy2.hdf5_input_marker  | ascotpy2.hdf5_input_boozer | \
                                ascotpy2.hdf5_input_mhd
                                
        ascotpy2.hdf5_interface_read_input(
            ctypes.byref(self.sim),
            what_to_read,
            ctypes.byref(self.B_offload_array),
            ctypes.byref(self.E_offload_array),
            ctypes.byref(self.plasma_offload_array),
            ctypes.byref(self.neutral_offload_array),
            ctypes.byref(self.wall_offload_array),
            ctypes.byref(self.boozer_offload_array),
            ctypes.byref(self.mhd_offload_array),
            ctypes.byref(self.prt)
            )
