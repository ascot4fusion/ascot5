"""Simulation routines and tools to prepare for the simulation.
"""
import ctypes

from mpi4py import MPI
MPI=None

from ..libascot import LIBASCOT
from ..data.bfield import Bfield
from ..data.efield import Efield
from ..data.plasma import Plasma
from ..data.neutral import Neutral
from ..data.wall import Wall
from ..data.boozer import BoozerMap
from ..data.mhd import Mhd
from ..data.asigma import Atomic
from ..data.nbi import Nbi


VARIANT_TYPES = {
    "EfieldCartesian": ("ETC", 0),
    "EfieldRadialPotential": ("E1DS", 1),

    "Plasma1D": ("Plasma_1D", 0),
    "Plasma1DDynamic": ("Plasma_1Dt", 1),

    "Neutral1D": ("N0_1D", 0),
    "Neutral3D": ("N0_3D", 1),

    "MHDStationary": ("MHD_STAT", 0),
    "MHDDynamic": ("MHD_NONSTAT", 1),

    "AtomicData": ("asigma_loc", 0),

    "Wall2D": ("w2d", 0),
    "Wall3D": ("w3d", 1),
}
"""Dictionary that connects the variant name to it's type in C and
the corresponding enum value."""


class Diag(ctypes.Structure):
    pass

class Simulate():

    #particle_input_ml_to_state(&p->p_ml, ps, Bdata)

    def init(self, mpirank):
        """
        - Check what input is needed
        - Generate Run Group
        - Initialize marker state (and randomize it)
        - Broadcast input data
        - Initialize sim struct
        - Initialize diagnostics
        - offload input and output data (?).
        """
        if mpirank == 0:
            required = ["bfield", "efield"]
            options = self.data.options.active
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
        br, bphi, bz = self.input_eval(mrk.r, mrk.phi, mrk.z, mrk.time, "br", "bphi", "bz")
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
            mrk._struct_[i].B_r = br
            mrk._struct_[i].B_phi = bphi
            mrk._struct_[i].B_z = bz
            mrk._struct_[i].B_r_dr = 0.0
            mrk._struct_[i].B_r_dphi = 0.0
            mrk._struct_[i].B_r_dz = 0.0
            mrk._struct_[i].B_phi_dr = 0.0
            mrk._struct_[i].B_phi_dphi = 0.0
            mrk._struct_[i].B_phi_dz = 0.0
            mrk._struct_[i].B_z_dr = 0.0
            mrk._struct_[i].B_z_dphi = 0.0
            mrk._struct_[i].B_z_dz = 0.0
