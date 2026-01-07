import imaspy
import libmuscle
import ymmsl

import ascot

from imas import imasdef

from mpi4py import MPI


def read_file(file_path: str):
    with open(file_path, 'r') as file:
        file_content = file.read()

    return file_content

def main():

    # ----  IN/OUT IDSes  ----

    ids_factory = imaspy.IDSFactory()
    equilibrium_3d_in = ids_factory.equilibrium()
    equilibrium_in = ids_factory.equilibrium()
    core_profiles_in = ids_factory.core_profiles()
    distribution_sources_in = ids_factory.distribution_sources()
    wall2d_in = ids_factory.wall()
    wall3d_in = ids_factory.wall()
    distributions_out = ids_factory.distributions()

    # Get parent communicator
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print('PROCESS [', rank, ':', size, ']')



    # ##################################################################################
    #                      ROOT NODE ACTIONS
    #          (only ROOT communicates with MUSCLE3)
    # ##################################################################################
    if rank == 0:

        print('ASCOT RUNS @ 0')

        # # # #  - - - M3 INITIALISATION - - - # # # #
        ports = { ymmsl.Operator.F_INIT: ["equilibrium_3d_in","equilibrium_in", "core_profiles_in","distribution_sources_in", "wall2d_in", "wall3d_in"],
                  ymmsl.Operator.O_F: ["distributions_out"]}

        m3_instance = libmuscle.Instance(ports)

        # # # #  - - - MUSCLE3 LOOP  - - - # # # #

        reuse_flag = m3_instance.reuse_instance()
        comm.bcast(reuse_flag, root=0)

        code_parameters_file = m3_instance.get_setting("code_parameters_file", "str")
        comm.bcast(code_parameters_file, root=0)

        code_parameters = read_file(code_parameters_file)

        while reuse_flag:
            print('M# LOOP  @ 0')
            # F_INIT
            # # # #  - - - RECEIVING INPUT IDSes - - - # # # #

            # equilibrium_3d_in
            m3_message = m3_instance.receive("equilibrium_3d_in")
            ids_bytes = m3_message.data
            t_cur = m3_message.timestamp

            # Broadcast data to all other processes
            comm.bcast(ids_bytes, root=0)

            equilibrium_3d_in.deserialize(ids_bytes)

            # equilibrium_in
            m3_message = m3_instance.receive("equilibrium_in")
            ids_bytes = m3_message.data
            t_cur = m3_message.timestamp

            # Broadcast data to all other processes
            comm.bcast(ids_bytes, root=0)

            equilibrium_in.deserialize(ids_bytes)

            # core_profiles_in
            m3_message = m3_instance.receive("core_profiles_in")
            ids_bytes = m3_message.data
            t_cur = m3_message.timestamp

            # Broadcast data to all other processes
            comm.bcast(ids_bytes, root=0)

            core_profiles_in.deserialize(ids_bytes)

            # distribution_sources_in
            m3_message = m3_instance.receive("distribution_sources_in")

            ids_bytes = m3_message.data
            t_cur = m3_message.timestamp

            # Broadcast data to all other processes
            comm.bcast(ids_bytes, root=0)

            distribution_sources_in.deserialize(ids_bytes)

            # wall2d_in
            m3_message = m3_instance.receive("wall2d_in")

            ids_bytes = m3_message.data
            t_cur = m3_message.timestamp

            # Broadcast data to all other processes
            comm.bcast(ids_bytes, root=0)

            wall2d_in.deserialize(ids_bytes)

            # wall3d_in
            m3_message = m3_instance.receive("wall3d_in")

            ids_bytes = m3_message.data
            t_cur = m3_message.timestamp

            # Broadcast data to all other processes
            comm.bcast(ids_bytes, root=0)

            wall3d_in.deserialize(ids_bytes)

            # # # #  - - - ASCOT CODE CALL - - - # # # #
            print('ASCOT RUNS @ 0')
            ascot.ascot_main(equilibrium_3d_in,equilibrium_in, core_profiles_in, distribution_sources_in, wall2d_in, wall3d_in, distributions_out, code_parameters)
            print('ASCOT ENDS @ 0')
            # O_F
            # # # #  - - - SENDING OUTPUT IDSes - - - # # # #
            print('SENDING OUT')
            ids_bytes = distributions_out.serialize(imasdef.DEFAULT_SERIALIZER_PROTOCOL)

            m3_message = libmuscle.Message(t_cur, None, ids_bytes)

            print('SENDING OUT - M3')
            m3_instance.send("distributions_out", m3_message)
            print('SENDING OUT - ENDS')

            reuse_flag = m3_instance.reuse_instance()
            comm.bcast(reuse_flag, root=0)
            # > > > > THE END OF MUSCLE3 LOOP

    else:
        # ##################################################################################
        #                     'NON-ROOT' NODE ACTIONS
        #               (all processes, beside ROOT, are MUSCLE3 agnostic)
        # ##################################################################################
        # Non-rank-0 processes receive broadcasted serialized IDSes
        print(f'ASCOT RUNS @ {rank}')
        reuse_flag = comm.bcast(None, root=0)
        code_parameters_file = comm.bcast(None, root=0)

        code_parameters = read_file(code_parameters_file)

        while reuse_flag:

            ids_bytes = comm.bcast(None, root=0)
            equilibrium_3d_in.deserialize(ids_bytes)

            ids_bytes = comm.bcast(None, root=0)
            equilibrium_in.deserialize(ids_bytes)

            ids_bytes = comm.bcast(None, root=0)
            core_profiles_in.deserialize(ids_bytes)

            ids_bytes = comm.bcast(None, root=0)
            distribution_sources_in.deserialize(ids_bytes)

            ids_bytes = comm.bcast(None, root=0)
            wall2d_in.deserialize(ids_bytes)

            ids_bytes = comm.bcast(None, root=0)
            wall3d_in.deserialize(ids_bytes)

            # # # #  - - - ASCOT CODE CALL - - - # # # #
            ascot.ascot_main(equilibrium_3d_in, equilibrium_in, core_profiles_in, distribution_sources_in, wall2d_in, wall3d_in, distributions_out, code_parameters)

            reuse_flag = comm.bcast(None, root=0)
            # Do not send results back to rank 0


if __name__ == "__main__":
    main()



