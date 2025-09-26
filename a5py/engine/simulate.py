"""Simulation routines and tools to prepare for the simulation.
"""
import ctypes

#from mpi4py import MPI
#MPI=None

import numpy as np

from ..libascot import LIBASCOT
from ..data.options import Options, OptionsStruct
from ..data.bfield import Bfield
from ..data.efield import Efield
from ..data.plasma import Plasma
from ..data.neutral import Neutral
from ..data.wall import Wall
from ..data.boozer import Struct as BoozerMap
from ..data.mhd import Mhd
from ..data.asigma import Atomic
from ..data.nbi import NbiStruct

from a5py.data.marker.state import Structure

from .functions import init_fun, PTR_DOUBLE, PTR_INT
from .interpolate import evaluate

from a5py.composite.run import Run


class Rfof(ctypes.Structure):

    _fields_ = [
        ("rfof_input_params", ctypes.c_void_p),
        ("rfglobal", ctypes.c_void_p),
        ]


class SimData(ctypes.Structure):

    _fields_ = [
        ("B_data", Bfield),
        ("E_data", Efield),
        ("plasma_data", Plasma),
        ("neutral_data", Neutral),
        ("wall_data", Wall),
        ("boozer_data", ctypes.POINTER(BoozerMap)),
        ("mhd_data", Mhd),
        ("asigma_data", Atomic),
        ("nbi_data", NbiStruct),
        ("rfof_data", Rfof),
        ("orbit", ctypes.c_void_p),
        ("dist5d", ctypes.c_void_p),
        ("dist6d", ctypes.c_void_p),
        ("dist5drho", ctypes.c_void_p),
        ("dist6drho", ctypes.c_void_p),
        ("distcom", ctypes.c_void_p),
        ("transport_coefficient", ctypes.c_void_p),
        ("random_data", ctypes.c_void_p),
        ("mccc_data", ctypes.c_void_p),
        ("params", ctypes.POINTER(OptionsStruct)),
        ]

init_fun(
    "simulate",
    ctypes.c_int,
    ctypes.POINTER(Structure),
    ctypes.POINTER(SimData),
    )


def setup_run(params, **inputs):
    """Check that we have sufficient inputs to carry out this simulation.
    """
    required = {"bfield"}
    mf = params.simulation.simulation_mode == 4
    if params.physics.enable_orbit_following and not mf:
        required.add("efield")
    if params.physics.enable_icrh:
        required.add("rfof")
    if params.physics.enable_coulomb_collisions or params.endconditions.activate_energy_limits:
        required.add("plasma")
    if params.physics.enable_mhd:
        required.update({"boozer", "mhd"})

    missing_inputs = {req for req in required if not req in inputs.keys()}
    if missing_inputs:
        raise ValueError(f"Missing inputs: {missing_inputs}")

    return Run(inputs=inputs)


def setup_markers(mrk, bfield):
    axisr, axisz = 6.2, 0.0
    for i in range(len(mrk._cdata)):
        mrk._cdata[i].rprt = mrk._cdata[i].r
        mrk._cdata[i].phiprt = mrk._cdata[i].phi
        mrk._cdata[i].zprt = mrk._cdata[i].z
        mrk._cdata[i].p_r = 0.0
        mrk._cdata[i].p_phi = 0.0
        mrk._cdata[i].p_z = 0.0
        mrk._cdata[i].mass = 0.0
        mrk._cdata[i].charge = 0.0
        mrk._cdata[i].anum = 0
        mrk._cdata[i].znum = 0
        mrk._cdata[i].weight = 0.0
        mrk._cdata[i].theta = np.arctan2(mrk._cdata[i].z - axisz,
                                            mrk._cdata[i].r - axisr)
        mrk._cdata[i].endcond = 0
        mrk._cdata[i].walltile = 0
        mrk._cdata[i].cputime = 0.0
        mrk._cdata[i].mileage = 0.0
        mrk._cdata[i].mu = 0.0
        mrk._cdata[i].zeta = 0.0




def init(mpirank, simtype, params: Options, **inputs):
    """Initialize data for the simulation.

    This function does the following:

    - Check what input is needed according to the parameters and simulation
      type. (root)
    - Generate Run Group (root)
    - Initialize marker state (and randomize it) (all)
    - Broadcast input data (root/all)
    - Initialize sim struct (all)
    - Initialize diagnostics (all)
    """

    def bcast_data(data):
        if MPI is None:
            return data
        return MPI.COMM_WORLD.bcast(data, root=0)
    data = None
    if mpirank == 0:
        data = self.data.bfield.active.export()
    data = bcast_data(data)
    if mpirank > 0:
        self.data.create_bfieldanalytical(**data, activate=True)

    mrk = self.data.create_fieldlinemarker(
        ids=np.arange(100)+1,
        r=6 + np.random.rand(100),
        phi=np.zeros((100,)),
        z=np.zeros((100,)),
        direction=np.ones((100,)),
    )
    self.data.bfield.active.stage()
    self.data.options.active.stage()
    

    self._sim = Simulate.sim_data()
    self._sim.B_data.use(self.data.bfield.active)
    self._sim.params = ctypes.pointer(self.data.options.active._struct_)


def simulate(run):
    sim = SimData()
    sim.B_data.use(run.bfield)
    run.options.stage()
    sim.params = ctypes.pointer(run.options._cdata)
    LIBASCOT.simulate(len(run._diagnostics["endstate"]._cdata), run._diagnostics["endstate"]._cdata, ctypes.byref(sim))
    for i in range(len(run._diagnostics["endstate"]._cdata)):
        print(run._diagnostics["endstate"]._cdata[i].z)

class Simulate():
    pass
