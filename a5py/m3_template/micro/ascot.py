import imas

from mpi4py import MPI

from a5py._nonfunctional.distsource_fun import *


def check_ids(ids, rank, size):
    print('IDS [ ', ids.metadata.name, ' : ', rank, '/', size, '] time:' , ids.time)
    pass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                         MOCKS ASCOT MAIN METHOD
#   This is simple, sequential version of the code. To benefit from the parallel runs the data
#    should be
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def ascot_main(equilibrium_3d_in, equilibrium_in, distribution_sources_in, core_profiles_in, wall2d_in, wall3d_in, distributions_out, code_parameters):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # print('=======================================')
    # print('START OF PHYSICS CODE')

    distsource_run(equilibrium_3d_in, equilibrium_in, distribution_sources_in, core_profiles_in, wall2d_in, wall3d_in,distributions_out, code_parameters)

    
    print('rank, size', rank, size)
    #check_ids(equilibrium_in, rank, size)
    #check_ids(distribution_sources_in, rank, size)
    #check_ids(wall2d_in, rank, size)
    #check_ids(wall3d_in, rank, size)

    # print('PROCESS [', rank, ':', size, ']')
    # print('=======================================')
    print(code_parameters)

    # MANDATORY FLAG (UNIFORM TIME HERE)
    distributions_out.ids_properties.homogeneous_time = 1

    distributions_out.code.name   = 'EXAMPLE: code_restart'
    distributions_out.code.version   = '1.0'

    distributions_out.time.resize(1)
    distributions_out.time[0] = 5

    # FINAL DISPLAY
    # print('END OF PHYSICS CODE')
    # print('=======================================')






