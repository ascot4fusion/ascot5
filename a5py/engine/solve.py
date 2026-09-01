"""Simulation routines and tools to prepare for the simulation.
"""
import ctypes

import numpy as np

from a5py.exceptions import AscotDataException

from ..libascot import LIBASCOT
from ..data.options import SimulationOptions, OptionsStruct
import a5py.data.orbit as orbit
import a5py.data.hist as hist

from a5py.data import (
    AscotData, Bfield, Efield, Plasma, Neutral, Wall, Mhd, Atomic,
    NbiStruct,
    )
from a5py.data.boozer import Struct as BoozerMap
from a5py.data.marker.state import Structure, MarkerState

from .functions import init_fun, PTR_DOUBLE, PTR_INT
from .interpolate import evaluate

from a5py.composite.run import Run



class Rfof(ctypes.Structure):

    _fields_ = [
        ("rfof_input_params", ctypes.c_void_p),
        ("rfglobal", ctypes.c_void_p),
        ]

class Diag(ctypes.Structure):

    _fields_ = [
        ("orbit", ctypes.POINTER(orbit.Struct)),
        ("hist", ctypes.POINTER(ctypes.POINTER(hist.Struct))),
        ("nhist", ctypes.c_size_t),
        ]


class SimData(ctypes.Structure):

    _fields_ = [
        ("mhd", Mhd),
        ("wall", Wall),
        ("rfof", ctypes.POINTER(Rfof)),
        ("bfield", Bfield),
        ("efield", Efield),
        ("plasma", Plasma),
        ("boozer", ctypes.POINTER(BoozerMap)),
        ("atomic", ctypes.POINTER(Atomic)),
        ("neutral", Neutral),
        ("options", ctypes.POINTER(OptionsStruct)),
        ("diagnostics", Diag),
        ("random_data", ctypes.c_void_p),
        ("mccc_data", ctypes.c_void_p),
        ]

init_fun(
    "ascot_solve_distribution",
    ctypes.POINTER(SimData),
    ctypes.c_size_t,
    ctypes.POINTER(Structure),
    )

def _get_marker_idx(comm, nmrk):
    if comm is None:
        return np.arange(nmrk)
    size, rank = comm.Get_size(), comm.Get_rank()
    base = nmrk // size
    remainder  = nmrk % size

    start = rank * base + min(rank, remainder)
    count = base + (1 if rank < remainder else 0)
    end = start + count
    return np.arange(start, end)

def setup(inputs, params, comm=None):
    """Prepare inputs for the simulation.

    Parameters
    ----------
    inputs : dict
        Input data.
    comm : MPI.Comm, *optional*
        MPI communicator.
    """
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

    idx = _get_marker_idx(comm, inputs["marker"].n)
    bjac = np.zeros((12,idx.size))
    (
        bjac[0,:], bjac[1:,:], bjac[2:,:], bjac[3:,:], bjac[4:,:], bjac[5:],
        bjac[6:,:], bjac[7:,:], bjac[8:,:], bjac[9:,:], bjac[10:,:], bjac[11,:],
     ) = evaluate(
        inputs["marker"].r[idx],
        inputs["marker"].phi[idx],
        inputs["marker"].z[idx],
        inputs["marker"].time[idx],
        "br", "bphi", "bz", "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz",
        bfield=inputs["bfield"],
        )
    state = MarkerState.from_markerinput(inputs["marker"], bjac.T, idx)
    run = Run.make_fresh(
        mrk=state,
        inputs=inputs,
        params=params,
        )
    return run, inputs, []


def execute(run, time=None):
    """Initialize the simulation object and run the simulation.

    Parameters
    ----------
    run : Run
        The run to simulate.
    time : float, optional
        Real time before the simulation is paused.
    """
    sim = SimData()
    for inp in ["bfield", "efield", "plasma", "neutral", "boozer", "mhd"]:
        if hasattr(run, inp):
            getattr(sim, inp).use(getattr(run, inp))
    options = OptionsStruct.from_params(run.options)
    sim.options = ctypes.pointer(options)
    if "orbit" in run._diagnostics:
        sim.diagnostics.orbit = ctypes.pointer(run._diagnostics["orbit"]._cdata)

    i = 0
    while True:
        if f"hist_{i}" not in run._diagnostics:
            break
        i += 1
    sim.diagnostics.nhist = i
    ptr_array_type = ctypes.POINTER(hist.Struct) * sim.diagnostics.nhist
    sim.diagnostics.hist = ptr_array_type()

    for i in range(sim.diagnostics.nhist):
        obj = run._diagnostics[f"hist_{i}"]._cdata
        sim.diagnostics.hist[i] = ctypes.pointer(obj)

    mrk = run._diagnostics["endstate"]
    LIBASCOT.ascot_solve_distribution(ctypes.byref(sim), mrk.n, mrk._cdata)


def finalize(run, unstage, comm):
    """Combine results to root process and free any used resources."""
    root = True if comm is None else comm.Get_rank() == 0
    if root:
        pass

    if comm is not None:
        mrk = run._diagnostics["endstate"]
        n_local = mrk.n
        counts = comm.gather(n_local, root=0)

        if root:
            displs = np.cumsum([0] + counts[:-1])
            n_total = sum(counts)

            cdata = (Structure * n_total)()
            base_ptr = ctypes.addressof(cdata)

            # copy own data first
            ctypes.memmove(
                base_ptr + displs[0] * ctypes.sizeof(Structure),
                ctypes.addressof(mrk._cdata),
                n_local * ctypes.sizeof(Structure)
            )

            # receive others
            for i in range(1, comm.Get_size()):
                offset = displs[i] * ctypes.sizeof(Structure)
                comm.Recv(
                    [base_ptr + offset,
                    counts[i] * ctypes.sizeof(Structure),
                    comm.BYTE],
                    source=i
                )
        else:
            comm.Send(
                [ctypes.addressof(mrk._cdata),
                n_local * ctypes.sizeof(Structure),
                comm.BYTE],
                dest=0
            )
