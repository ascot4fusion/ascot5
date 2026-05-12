from mpi4py import MPI

import unyt
import numpy as np
import pytest
from unittest.mock import patch

from a5py import Ascot, SimulationOptions
from a5py.templates import PremadeMagneticField


def test_solve_setup_inputs():
    pass

@pytest.mark.mpi
def test_solve_mock_execute():
    """Test solve function without calling ``execute``."""
    a5 = Ascot()

    #comm = MPI.COMM_WORLD
    comm = None
    rank = 0 if comm is None else comm.Get_rank()
    root = rank == 0

    if root:
        template = PremadeMagneticField(a5, field="iter-baseline")
        template.create_input()

        parameters = SimulationOptions(
            enable_orbit_following=True,
            activate_simulation_time_limits=True,
            simulation_mode=4,
            activate_real_time_limit=True,
            max_real_time=5.0,
        )

        arr = np.ones(100,)
        a5.data.create_fieldlinemarker(r=arr*6.3*unyt.m, z=arr*0.0*unyt.m)

    # with (patch("a5py.composite.ascot.execute") as mock_execute,
    #      patch("a5py.composite.ascot.finalize") as mock_finalize):
    #     if root:
    #         run = a5.simulate(params=parameters)
    #     else:
    #         pass
    #         #a5.simulate()
    #     #mock_execute.assert_called_once()
    #     #mock_finalize.assert_called_once()

    with patch("a5py.composite.ascot.executee") as mock_execute:
        if root:
            run = a5.simulate(params=parameters, comm=comm)
            print(run.getstate("r", state="end")[0].size)
        else:
            a5.simulate(comm=comm)
        #mock_execute.assert_called_once()
        #mock_finalize.assert_called_once()

    assert False


def test_solve():
    pass
