import imaspy
import libmuscle
from ymmsl import Operator

import tools

# ----  IN/OUT IDSes  ----
equilibrium_3d, equilibrium, core_profiles, distribution_sources, wall2d, wall3d = tools.get_inputs()
distributions = imaspy.IDSFactory().distributions()

# # # #  - - - M3 INITIALISATION - - - # # # #
ports = {Operator.O_I: ["equilibrium_3d_out","equilibrium_out", "core_profiles_out","distribution_sources_out", "wall2d_out", "wall3d_out"],
         Operator.S: ["distributions_in"]}

m3_instance = libmuscle.Instance(ports)

# # # #  - - - MUSCLE3 LOOP  - - - # # # #
while m3_instance.reuse_instance():

    # F_INIT
    t_max = m3_instance.get_setting('t_max', 'float')
    dt = m3_instance.get_setting('dt', 'float')

    t_cur = 0.0

    while t_cur + dt <= t_max:
        # O_I
        t_next = t_cur + dt

        if t_next + dt > t_max:
            t_next = None

        # # # #  - - - SENDING OUTPUT IDSes - - - # # # #
        tools.send_ids(m3_instance, equilibrium_3d, "equilibrium_3d_out", t_cur, t_next)
        tools.send_ids(m3_instance, equilibrium, "equilibrium_out", t_cur, t_next)
        tools.send_ids(m3_instance, core_profiles, "core_profiles_out", t_cur, t_next)
        tools.send_ids(m3_instance, distribution_sources, "distribution_sources_out", t_cur, t_next)
        tools.send_ids(m3_instance, wall2d, "wall2d_out", t_cur, t_next)
        tools.send_ids(m3_instance, wall3d, "wall3d_out", t_cur, t_next)

        # # # #  - - - RECEIVING INPUT IDSes - - - # # # #

        print('Waiting for distributions')
        # core_profiles_in
        m3_message = m3_instance.receive('distributions_in')
        print('DISTRIBUTIONS received')

        if m3_message.timestamp > t_cur + dt:
            m3_instance.error_shutdown('Received a message from the future!')

        ids_bytes = m3_message.data
        distributions.deserialize(ids_bytes)

        if len(distributions.time) < 1:
            m3_instance.error_shutdown("ERROR: Received empty TIME vector!")
            exit(1)

        print("TIME received ", distributions.time)

        t_cur = t_cur + dt
        # > > > > THE END OF MUSCLE3 LOOP
