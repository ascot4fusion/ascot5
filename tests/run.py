from mpi4py import MPI

import unyt
import numpy as np

from a5py import Ascot, SimulationOptions
from a5py.templates import PremadeMagneticField

if MPI is not None:
    comm = MPI.COMM_WORLD

rank = 0 if comm is None else comm.rank
isroot = rank == 0

if isroot:
    ascot = Ascot()
    template = PremadeMagneticField(ascot, field="iter-baseline")
    template.create_input()
    ascot.data.create_efieldcartesian(exyz=np.array([0., 0., 0.])*unyt.V/unyt.m)

    parameters = SimulationOptions.from_dict(
        simulation={
            "mode": "gyro-orbit",
            "timestep": 1e-11,
        },
        physics={
            "enable_orbit_following": True,
        },
        endconditions={
            "activate_simulation_time_limits": True,
            "max_mileage": 2e-9,
        },
        orbit={
            "collect": "interval",
            "buffer_size": 5,
            "interval": 1e-11,
        },
        histograms=[
            {
                "dimensions": [
                    ("ekin", 0.*unyt.eV, 1e7*unyt.eV, 5),
                    ("pitch", -1, 1, 3),
                    ],
                "charge_interval": (-1, -1,),
            },
        ]
    )

    ascot.data.create_guidingcentermarker(
        species="electron",
        r=np.array([6.3, 6.5])*unyt.m,
        z=0.*unyt.m,
        ekin=3.5e4*unyt.eV,
        pitch=0.5,
        )
    run = ascot.simulate(params=parameters, comm=comm)
else:
    ascot = Ascot()
    ascot.simulate(comm=comm)
