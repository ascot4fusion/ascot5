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
        params = comm.bcast(params, root=0)
        input_names = comm.bcast(None if not root else input_names, root=0)
        #input_names.remove("marker")

        for variant in input_names:
            exported_data = None if not root else inputs[variant].export()
            name = None if not root else inputs[variant].variant

            exported_data = comm.bcast(exported_data, root=0)
            name = comm.bcast(name, root=0)
            if not root:
                inputs[variant] = AscotData.from_export(name, exported_data)

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

    if comm is not None:
        from mpi4py import MPI
        mrk = run._diagnostics["endstate"]
        n_local = mrk.n
        counts = comm.gather(n_local, root=0)
        structure_size = ctypes.sizeof(Structure)

        if root:
            displs = np.cumsum([0] + counts[:-1])
            n_total = sum(counts)

            cdata = (Structure * n_total)()

            # Copy root's own data
            ctypes.memmove(
                ctypes.addressof(cdata),
                ctypes.addressof(mrk._cdata),
                n_local * structure_size
            )

            # Receive data from other ranks
            for i in range(1, comm.Get_size()):
                nbytes = counts[i] * structure_size
                offset = displs[i] * structure_size

                # Create a ctypes byte buffer pointing into cdata
                recv_buf = (ctypes.c_char * nbytes).from_buffer(
                    cdata, offset
                )

                comm.Recv(
                    [recv_buf, nbytes, MPI.BYTE],
                    source=i
                )
            run._diagnostics["endstate"] = MarkerState.from_params(cdata)

        else:
            nbytes = n_local * structure_size
            send_buf = (ctypes.c_char * nbytes).from_buffer(mrk._cdata)

            comm.Send(
                [send_buf, nbytes, MPI.BYTE],
                dest=0
            )

        i = 0
        while True:
            if f"hist_{i}" not in run._diagnostics:
                break
            values = run._diagnostics[f"hist_{i}"].values

            if root:
                values.flags.writeable = True
                comm.Reduce(MPI.IN_PLACE, values, op=MPI.SUM, root=0)
                values.flags.writeable = False
            else:
                comm.Reduce(values, None, op=MPI.SUM, root=0)
            i += 1

        if "orbit" in run._diagnostics:
            orbit = run._diagnostics["orbit"]

            # Number of elements on this process.
            local_count = np.array(orbit.id.size, dtype=np.int64)

            if root:
                counts = np.empty(comm.size, dtype=np.int64)
            else:
                counts = None

            # Gather the sizes first.
            comm.Gather(
                local_count,
                counts,
                root=0
            )

            for attr in ["r", "z", "phi", "p1", "p2", "p3", "mileage",
                         "id", "charge", "poincare", "simmode"]:
                local_array = np.ascontiguousarray(getattr(orbit, attr))
                if root:
                    # Displacement of each process' data in the final array.
                    displacements = np.zeros(comm.size, dtype=np.int64)

                    if comm.size > 1:
                        displacements[1:] = np.cumsum(counts[:-1])

                    # Allocate the final concatenated array.
                    total_size = int(np.sum(counts))

                    result = np.empty(
                        total_size,
                        dtype=local_array.dtype
                    )

                else:
                    result = None
                    displacements = None

                # Determine the MPI datatype from the NumPy dtype.
                mpi_dtype = MPI._typedict[local_array.dtype.char]

                # Gather the actual data.
                comm.Gatherv(
                    local_array,
                    (
                        result,
                        counts,
                        displacements,
                        mpi_dtype
                    ) if root else None,
                    root=0
                )

                if root:
                    setattr(orbit._cdata, attr + "_ref", result)
                    arr = getattr(orbit._cdata, attr + "_ref")
                    arr.flags.writeable = False
                    ctype = np.ctypeslib.as_ctypes_type(getattr(orbit._cdata, attr)._type_)
                    setattr(orbit._cdata, attr, arr.ctypes.data_as(ctypes.POINTER(ctype)))
