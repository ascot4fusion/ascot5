from mpi4py import MPI

import unyt
import numpy as np

from a5py import Ascot, SimulationOptions
from a5py.templates import PremadeMagneticField

if MPI is not None:
    comm = MPI.COMM_WORLD

rank = 0 if comm is None else comm.Get_rank()
root = rank == 0

if root:
    a5 = Ascot()
    template = PremadeMagneticField(a5, field="iter-baseline")
    template.create_input()

    parameters = SimulationOptions(
        enable_orbit_following=True,
        activate_simulation_time_limits=True,
        simulation_mode=4,
        activate_real_time_limit=True,
        max_real_time=1.0,
    )

    a5.data.create_fieldlinemarker(r=6.3*unyt.m, z=0.001*unyt.m)
    run = a5.simulate(params=parameters, comm=MPI.COMM_WORLD, time=3600)
    print(run.getstate("r", state="end"))
else:
    a5 = Ascot()
    a5.simulate(comm=MPI.COMM_WORLD)
