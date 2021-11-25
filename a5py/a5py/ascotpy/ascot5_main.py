'''
Created on Nov 16, 2021

@author: sjjamsa
'''

import ctypes

# generate ascotpy2 in the source folder with
# % clang2py -l ./libascot.so -o ../python/a5py/a5py/ascotpy/ascotpy2.py  particle.h hdf5_interface.h ascot5.h mpi_interface.h simulate.h ascot5_main.h offload.h diag.h libascot_mem.h --clang-args="-I/usr/include/hdf5/serial"
from a5py.ascotpy import ascotpy2



class ascot5_main(object):
    '''
    This class plays the role of the C-function ascot5_main.c main() function.
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
    offload_array         = ctypes.POINTER(ctypes.c_double)() # Some sort of joint array for all (?)
    
    diag_offload_array    = ctypes.POINTER(ctypes.c_double)() 

    offload_package       = ascotpy2.struct_c__SA_offload_package()

    # number of particles read
    n_tot   = ctypes.c_int32()
    
    #number of particles to simulate (?)
    nprts   = ctypes.c_int32()

    #number of particles gathered with MPI (?)
    n_gathered   = ctypes.c_int32()
    
    
    # QID for the new results.
    qid_str = (ctypes.c_char * 12)()
    qid     = ctypes.POINTER(ctypes.c_char)(qid_str)
    
    # pointer to the current particle
    prt          = ctypes.POINTER(ascotpy2.struct_c__SA_input_particle)()
    ps           = ctypes.POINTER(ascotpy2.struct_c__SA_particle_state)()
    ps_gathered  = ctypes.POINTER(ascotpy2.struct_c__SA_particle_state)()

    def __init__(self, input_filename=b'input.h5', output_filename=b'output.h5'):
        '''
        Constructor
        '''
        
        self.sim.hdf5_in  =  input_filename
        self.sim.hdf5_out = output_filename
        


    def init(self):

        # No argumets given at the moment.
        argc=ctypes.c_int32()
        argc.value=0
        argv=ctypes.POINTER(ctypes.c_char)() # This points to null.

        
        ascotpy2.mpi_interface_init(argc,
                                    ctypes.byref(argv),
                                    ctypes.byref(self.sim), 
                                    ctypes.byref(self.mpi_rank),
                                    ctypes.byref(self.mpi_size),
                                    ctypes.byref(self.mpi_root)  ) 

        ascotpy2.hdf5_generate_qid(self.qid)

        
    def read_input(self,what_to_read=None):
        
        if what_to_read is None:         
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
            ctypes.byref(self.prt),
            ctypes.byref(self.n_tot)
            )

    def offload(self):        
        '''
        int offload(
        sim_offload_data *sim,
        real** B_offload_array,
        real** E_offload_array,
        real** plasma_offload_array,
        real** neutral_offload_array,
        real** wall_offload_array,
        real** boozer_offload_array,
        real** mhd_offload_array,
        int n_tot,
        int mpi_rank,
        int mpi_size,
        int mpi_root,
        char *qid,
        int *nprts,
        input_particle **p,
        int* n_gathered,
        real **offload_array,
        offload_package *offload_data,
        particle_state** ps,
        real** diag_offload_array_mic0,
        real** diag_offload_array_mic1,
        real** diag_offload_array_host
);
        '''
        '''
        offload.argtypes = [
        ctypes.POINTER(struct_c__SA_sim_offload_data), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, 
        ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_int32), 
        ctypes.POINTER(ctypes.POINTER(struct_c__SA_input_particle)), ctypes.POINTER(ctypes.c_int32), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(struct_c__SA_offload_package), 
        ctypes.POINTER(ctypes.POINTER(struct_c__SA_particle_state)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
        '''
        

        retval = ascotpy2.offload(
            ctypes.byref(self.sim),
            ctypes.byref(self.B_offload_array),
            ctypes.byref(self.E_offload_array),
            ctypes.byref(self.plasma_offload_array),
            ctypes.byref(self.neutral_offload_array),
            ctypes.byref(self.wall_offload_array),
            ctypes.byref(self.boozer_offload_array),
            ctypes.byref(self.mhd_offload_array),
            self.n_tot,  self.mpi_rank, self.mpi_size, self.mpi_root,
            self.qid, ctypes.byref(self.nprts),
            ctypes.byref(self.prt), ctypes.byref(self.n_gathered),
            ctypes.byref(self.offload_array),
            ctypes.byref(self.offload_package),
            ctypes.byref(self.ps),
            ctypes.byref(self.diag_offload_array)
            )

        if retval != 0:
            print('offload() returned',retval)

            '''
            int cleanup( sim_offload_data sim,    particle_state* ps,     particle_state* ps_gathered,
            real** diag_offload_array,
            real* offload_array,
            offload_package *offload_data
            ){
            '''
            
            ascotpy2.cleanup(self.sim, 
                             self.ps, self.ps_gathered,
                             ctypes.byref(self.diag_offload_array),
                             self.offload_array,
                             self.offload_package)
                             
    def run_simulation(self):
         
        retval = ascotpy2.run(self.nprts,
                              self.mpi_rank,
                              self.ps, 
                              self.offload_array,
                              self.diag_offload_array,
                              ctypes.byref(self.sim), 
                              self.offload_package)   
        
        return retval
        
    def gather_output(self):
        
        ascotpy2.gather_output(self.ps, ctypes.byref(self.ps_gathered),
                               ctypes.byref(self.n_gathered),
                               self.n_tot, self.mpi_rank, self.mpi_size, self.mpi_root,
                               self.sim, self.diag_offload_array)
        
    def print_marker_summary(self):
        
        ascotpy2.marker_summary(self.ps_gathered, self.n_gathered)
    
    def write_output_h5(self):
        
        '''
        int write_output(sim_offload_data sim, int mpi_rank, int mpi_root,
        particle_state *ps_gathered, int n_gathered,
        real* diag_offload_array);
        '''
        
        retval = ascotpy2.write_output(
            self.sim, 
            self.mpi_rank, self.mpi_root, 
            self.ps_gathered, self.n_gathered,
            self.diag_offload_array)
        
        if retval != 0:
            print('offload() returned',retval)
            
             
    def free_c(self):
        
        ascotpy2.offload_free_offload( 
            ctypes.byref(self.offload_package),
            ctypes.byref(self.offload_array) ) 
        
        ascotpy2.diag_free_offload(
            ctypes.byref(self.sim.diag_offload_data),
            ctypes.byref(self.diag_offload_array)
            )
        
        ascotpy2.free_ps(self.ps)
        ascotpy2.free_ps(self.ps_gathered)
        
    
    def finalize(self):
        ascotpy2.mpi_interface_finalize()
        
