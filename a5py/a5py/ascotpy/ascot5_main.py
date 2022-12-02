'''
Created on Nov 16, 2021

@author: sjjamsa
'''

import ctypes
from math import pi

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
        
    def inject_wall_2d(self,R,z):

        if( len(R) != len(z) ):
            raise ValueError("R and z must be of equal length.")

        nelements = len(R)

        # Create temporary variable for the wall
        R_c = (ctypes.c_double * nelements)()
        z_c = (ctypes.c_double * nelements)()

        for i in range(nelements):
            R_c[i] = R[i]
            z_c[i] = z[i]

        ascotpy2.hdf5_wall_2d_to_offload(
            ctypes.byref(self.sim.wall_offload_data.w2d),
            ctypes.byref(self.wall_offload_array),
            nelements,
            R_c, z_c
            )

        self.sim.wall_offload_data.type = ascotpy2.wall_type_2D
        
        ascotpy2.wall_init_offload(
            ctypes.byref(self.sim.wall_offload_data),
            self.wall_offload_array
            )

    def inject_wall_3d(self,x1x2x3,y1y2y3,z1z2z3):

        nelements = int(x1x2x3.shape[0])

        # Create temporary variable for the wall
        x1x2x3_c = (ctypes.c_double * (3*nelements) )()
        y1y2y3_c = (ctypes.c_double * (3*nelements) )()
        z1z2z3_c = (ctypes.c_double * (3*nelements) )()

        # Let's hope the ordering is correct...
        x1x2x3_c[:] = x1x2x3.flatten()[:]
        y1y2y3_c[:] = y1y2y3.flatten()[:]
        z1z2z3_c[:] = z1z2z3.flatten()[:]

        
        ascotpy2.hdf5_wall_3d_to_offload(
            ctypes.byref(self.sim.wall_offload_data.w3d),
            ctypes.byref(self.wall_offload_array),
            nelements,
            x1x2x3_c, y1y2y3_c, z1z2z3_c,
            )

        self.sim.wall_offload_data.type = ascotpy2.wall_type_3D

        ascotpy2.wall_init_offload(
            ctypes.byref(self.sim.wall_offload_data),
            self.wall_offload_array
            )

        
    
    def inject_BSTS(self,bsts):
        # bsts is the dictionary that comes from reading the hdf5



        # 1. First fill in the meta-data struct
        #--------------------------------------
        # Mimic the C-function    hdf5_bfield_read_STS()

        #phimin/max deg2rad
        
        #B_STS_offload_data
        #sts = struct_c__SA_B_STS_offload_data()
        BSTS = self.sim.B_offload_data.BSTS

        BSTS.psigrid_n_r     = bsts['psi_nr'][0]
        BSTS.psigrid_n_z     = bsts['psi_nz'][0]
        BSTS.psigrid_n_phi   = bsts['psi_nphi'][0]
        BSTS.psigrid_r_min   = bsts['psi_rmin'][0]
        BSTS.psigrid_r_max   = bsts['psi_rmax'][0]
        BSTS.psigrid_z_min   = bsts['psi_zmin'][0]
        BSTS.psigrid_z_max   = bsts['psi_zmax'][0]
        BSTS.psigrid_phi_min = bsts['psi_phimin'][0] * pi * 2.0 / 360.0
        BSTS.psigrid_phi_max = bsts['psi_phimax'][0] * pi * 2.0 / 360.0

        BSTS.Bgrid_n_r       = bsts['b_nr'][0]
        BSTS.Bgrid_n_z       = bsts['b_nz'][0]
        BSTS.Bgrid_n_phi     = bsts['b_nphi'][0]
        BSTS.Bgrid_r_min     = bsts['b_rmin'][0]
        BSTS.Bgrid_r_max     = bsts['b_rmax'][0]
        BSTS.Bgrid_z_min     = bsts['b_zmin'][0]
        BSTS.Bgrid_z_max     = bsts['b_zmax'][0]
        BSTS.Bgrid_phi_min   = bsts['b_phimin'][0] * pi * 2.0 / 360.0
        BSTS.Bgrid_phi_max   = bsts['b_phimax'][0] * pi * 2.0 / 360.0

        BSTS.psi0            = bsts['psi0'][0]
        BSTS.psi1            = bsts['psi1'][0]


        BSTS.n_axis          = bsts['axis_nphi'][0]
        BSTS.axis_min        = bsts['axis_phimin'][0] * pi * 2.0 / 360.0
        BSTS.axis_max        = bsts['axis_phimax'][0] * pi * 2.0 / 360.0
        # BSTS.axis_grid       = bsts['']   # Not really used
        
        # 2. Get the right sized offload array
        #--------------------------------------
        # Does this deed to be allocated by C code or can we do it with a python array?

        '''
        /* Allocate offload_array storing psi and the three components of B */
        int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z
           * offload_data->psigrid_n_phi;
        int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
           * offload_data->Bgrid_n_phi;
        int axis_size = offload_data->n_axis;

        *offload_array = (real*) malloc((psi_size + 3 * B_size + 2 * axis_size)
                                    * sizeof(real));
        offload_data->offload_array_length = psi_size + 3 * B_size + 2 * axis_size;
        '''

        offload_size = 0

        npsi = bsts['psi_nr'][0] * bsts['psi_nz'][0] * bsts[ 'psi_nphi'][0]
        nB   = bsts[  'b_nr'][0] * bsts[  'b_nz'][0] * bsts[   'b_nphi'][0]
        naxis=                                         bsts['axis_nphi'][0]
        
        # psi_size
        offload_size +=  1 * npsi
        # B_size
        offload_size +=  3 * nB
        # axis_size
        offload_size +=  2 * naxis
        
        BSTS.offload_array_length = offload_size

        # Python side allocation
        #-----------------------
        # B_offload_array = (ctypes.c_double * (offload_size) )()
        #self.B_offload_array.contents = B_offload_array


        # C side allocation
        #------------------
        self.B_offload_array =  ascotpy2.libascot_allocate_reals(offload_size)

        # Cast the pointer into an array
        B_offload_array = ctypes.cast( self.B_offload_array, ctypes.POINTER(ctypes.c_double*offload_size) )[0]
        
        
        # 3. copy the large arrays to the offload array
        #----------------------------------------------

        offset = 0
        order = 'F'

        B_offload_array[(offset):(offset+nB  )] = bsts['br'].flatten(order=order)
        offset += nB
        
        B_offload_array[(offset):(offset+nB  )] = bsts['bphi'].flatten(order=order)
        offset += nB

        B_offload_array[(offset):(offset+nB  )] = bsts['bz'].flatten(order=order)
        offset += nB

        B_offload_array[(offset):(offset+npsi)] = bsts['psi'].flatten(order=order)
        offset += npsi

        B_offload_array[(offset):(offset+naxis)] = bsts['axisr']
        offset += naxis

        B_offload_array[(offset):(offset+naxis)] = bsts['axisz']
        offset += naxis





        # 4. Set the correct data type
        #------------------------------

        self.sim.B_offload_data.type = ascotpy2.B_field_type_STS

        # 5. Do the init offload
        #-----------------------

        # B_field_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_field_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]

        print('B_field_init_offload')
        ascotpy2.B_field_init_offload(
            ctypes.byref(self.sim.B_offload_data),
            ctypes.byref(self.B_offload_array)
            )
        print('                    ..complete')

        
    
