"""Simulation routines and tools to prepare for the simulation.
"""
import ctypes

from mpi4py import MPI
MPI=None

from ..libascot import LIBASCOT
from ..data.options import Options
from ..data.bfield import Bfield
from ..data.efield import Efield
from ..data.plasma import Plasma
from ..data.neutral import Neutral
from ..data.wall import Wall
from ..data.boozer import BoozerMap
from ..data.mhd import Mhd
from ..data.asigma import Atomic
from ..data.nbi import NbiStruct

from .functions import init_fun, PTR_DOUBLE, PTR_INT




class Rfof(ctypes.Structure):
    #_pack_ = 1
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
        ("boozer_data", ctypes.POINTER(BoozerMap.Struct)),
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
        ("params", ctypes.POINTER(Options.Struct)),
        ]


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
    if mpirank == 0:
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

        missing_inputs = {req for req in required if not req in input.keys()}
        if missing_inputs:
            raise ValueError(f"Missing inputs: {missing_inputs}")

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

    import numpy as np
    from a5py.data.marker.cstructs import particle_state
    mrk = self.data.create_fieldlinemarker(
        ids=np.arange(100)+1,
        r=6 + np.random.rand(100),
        phi=np.zeros((100,)),
        z=np.zeros((100,)),
        direction=np.ones((100,)),
    )
    self.data.bfield.active.stage()
    self.data.options.active.stage()
    from .interpolate import evaluate
    qnt_name = ["br", "bphi", "bz", "brdr", "brdphi", "brdz",
                "bphidr", "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz"]
    qnt = evaluate(mrk.r, mrk.phi, mrk.z, mrk.time, *qnt_name, bfield=self.data.bfield.active)
    qnt = {n: q for n, q in zip(qnt_name, qnt)}
    axisr, axisz, rho = 6.2, 0.0, 0.0
    for i in range(len(mrk._struct_)):
        mrk._struct_[i].rprt = mrk._struct_[i].r
        mrk._struct_[i].phiprt = mrk._struct_[i].phi
        mrk._struct_[i].zprt = mrk._struct_[i].z
        mrk._struct_[i].p_r = 0.0
        mrk._struct_[i].p_phi = 0.0
        mrk._struct_[i].p_z = 0.0
        mrk._struct_[i].mass = 0.0
        mrk._struct_[i].charge = 0.0
        mrk._struct_[i].anum = 0
        mrk._struct_[i].znum = 0
        mrk._struct_[i].weight = 0.0
        mrk._struct_[i].theta = np.arctan2(mrk._struct_[i].z - axisz,
                                            mrk._struct_[i].r - axisr)
        mrk._struct_[i].endcond = 0
        mrk._struct_[i].walltile = 0
        mrk._struct_[i].cputime = 0.0
        mrk._struct_[i].mileage = 0.0
        mrk._struct_[i].mu = 0.0
        mrk._struct_[i].zeta = 0.0
        mrk._struct_[i].rho = rho
        mrk._struct_[i].B_r = qnt["br"][i]
        mrk._struct_[i].B_phi = qnt["bphi"][i]
        mrk._struct_[i].B_z = qnt["bz"][i]
        mrk._struct_[i].B_r_dr = qnt["brdr"][i]
        mrk._struct_[i].B_r_dphi = qnt["brdphi"][i]
        mrk._struct_[i].B_r_dz = qnt["brdz"][i]
        mrk._struct_[i].B_phi_dr = qnt["bphidr"][i]
        mrk._struct_[i].B_phi_dphi = qnt["bphidphi"][i]
        mrk._struct_[i].B_phi_dz = qnt["bphidz"][i]
        mrk._struct_[i].B_z_dr = qnt["bzdr"][i]
        mrk._struct_[i].B_z_dphi = qnt["bzdphi"][i]
        mrk._struct_[i].B_z_dz = qnt["bzdz"][i]

    self._sim = Simulate.sim_data()
    self._sim.B_data.use(self.data.bfield.active)
    self._sim.params = ctypes.pointer(self.data.options.active._struct_)


def simulate(sim, state):
    init_fun(
        "simulate",
        ctypes.c_int,
        ctypes.POINTER(particle_state),
        ctypes.POINTER(Simulate.sim_data),
        )

    for field_name, field_type in Simulate.sim_data._fields_:
        print(field_name, getattr(Simulate.sim_data, field_name).offset)


    LIBASCOT.simulate(100, state, ctypes.byref(sim))
    for i in range(len(state)):
        print(state[i].z)


class Simulate():

    pass

    
