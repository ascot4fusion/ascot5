"""Simulation routines and tools to prepare for the simulation.
"""
import ctypes

#from mpi4py import MPI
#MPI=None

import numpy as np

from a5py.exceptions import AscotDataException

from ..libascot import LIBASCOT
from ..data.options import Options, OptionsStruct

from a5py.data import (
    AscotData, Bfield, Efield, Plasma, Neutral, Wall, Mhd, Atomic,
    NbiStruct,
    )
from a5py.data.boozer import Struct as BoozerMap
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


def simulate(
        data, time=None, run=None, params=None, comm=None, **priority_inputs,
        ):
    """Let the markers loose.

    The simulation can either be newly initialized or the simulation can be
    continued from a previous run.

    For new simulation, options must be provided. By default the inputs that are
    set as `active` are used but these can be overridden with ``**inputs``:

    .. code-block:: python

        simulate(options=Options(...), bfield=mybfield)

    For continuing a simulation, the run must be provided and no ``options`` or
    ``**inputs``:

    .. code-block:: python

        simulate(run=run)

    In both cases, ``time`` can be supplied to control the real time for how
    long the simulation should run.


    Parameters
    ----------
    run : Run
        The run to simulate.
    """
    root = True if comm is None else comm.Get_rank() == 0

    inputs = {}
    if root:
        inputs = preflight_check(data, params, **priority_inputs)

    inputs, unstage = setup_inputs(inputs, comm)
    run = setup_run(params, inputs, comm)
    execute(run, time)
    finalize(run, unstage, comm)
    # Add to tree

    return run


def preflight_check(data, params, **priority_inputs):
    """Check that inputs are consistent and return only those inputs that are
    needed.

    Parameters
    ----------
    data : :class:`.AscotData`
        Simulation data structure whose `active` inputs will be used
    params : :class:`.Options`
        Simulation parameters.
    priority_inputs : dict
        Inputs to use in place of the ones that are `active` in ``data``.

    Returns
    -------
    inputs : dict[str, :class:`.InputVariant`]
        Inputs that will be used in the simulation.
    """
    preflight_check_parameters(params)

    required = {"bfield", "marker"}
    mf = params.simulation.simulation_mode == 4
    if params.physics.enable_orbit_following and not mf:
        required.add("efield")

    if params.physics.enable_icrh:
        required.add("rfof")

    if (params.physics.enable_coulomb_collisions
            or params.endconditions.activate_energy_limits):
        required.add("plasma")

    if params.physics.enable_mhd:
        required.update({"boozer", "mhd"})

    inputs = {}
    for req in required:
        try:
            inputs[req] = data[req].active
        except AscotDataException:
            if req in priority_inputs:
                inputs[req] = priority_inputs[req]

    missing_inputs = {req for req in required if not req in inputs.keys()}
    if missing_inputs:
        raise ValueError(f"Missing inputs: {missing_inputs}")

    return inputs


def preflight_check_parameters(params):
    """Check that options are internally consistent."""
    return True


def setup_inputs(inputs, comm=None):
    root = True if comm is None else comm.Get_rank() == 0
    if root:
        input_names = list(inputs.keys())

    if comm is not None:
        input_names = comm.bcast(None if not root else input_names, root=0)
        input_names.remove("marker")

        for variant in input_names:
            exported_data = None if not root else inputs[variant].export()

            exported_data = comm.bcast(exported_data, root=0)
            if not root:
                inputs[variant] = AscotData.from_export(variant, exported_data)

    inputs["marker"] = setup_markers(inputs["marker"], inputs["bfield"])
    return inputs, []


def setup_markers(mrk, bfield):
    axisr, axisz = 6.2, 0.0
    if True:#isinstance(mrk, FieldlineMarker):
        for i in range(len(mrk._cdata)):
            mrk._cdata[i].rprt = mrk._cdata[i].r
            mrk._cdata[i].phiprt = mrk._cdata[i].phi
            mrk._cdata[i].zprt = mrk._cdata[i].z
            mrk._cdata[i].pr = 0.0
            mrk._cdata[i].pphi = 0.0
            mrk._cdata[i].pz = 0.0
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
            mrk._cdata[i].ekin = 0.0
            mrk._cdata[i].zeta = 0.0

    return mrk


def setup_run(params, inputs, comm=None):
    """Create run object and setup diagnostics."""
    inputs.update({"options": params})
    run = Run(inputs=inputs)
    run._setup(inputs["marker"]._cdata)
    run.options.stage()
    return run


def execute(run, time):
    sim = SimData()
    sim.B_data.use(run.bfield)
    sim.params = ctypes.pointer(run.options._cdata)
    mrk = run._diagnostics["endstate"]._cdata
    LIBASCOT.simulate(len(mrk), mrk, ctypes.byref(sim))


def finalize(run, unstage, comm):
    """Combine results to root process and free any used resources."""
    pass


class Simulate():
    pass

    def resume():
        pass
