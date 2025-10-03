"""Defines the atomic reaction input class and the corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from a5py import utils
from a5py.libascot import (
    LIBASCOT, DataStruct, interp1D_data, interp2D_data, interp3D_data, init_fun,
    )
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in asigma_loc.h."""

    _fields_ = [
        ('N_reac', ctypes.c_int32),
        ('z_1', ctypes.POINTER(ctypes.c_int32)),
        ('a_1', ctypes.POINTER(ctypes.c_int32)),
        ('z_2', ctypes.POINTER(ctypes.c_int32)),
        ('a_2', ctypes.POINTER(ctypes.c_int32)),
        ('reac_type', ctypes.POINTER(ctypes.c_int32)),
        ('sigma', ctypes.POINTER(interp1D_data)),
        ('sigmav', ctypes.POINTER(interp2D_data)),
        ('BMSsigmav', ctypes.POINTER(interp3D_data)),
        ]


@Leaf.register
class Asigma_loc(InputVariant):
    """Local atomic data."""

    @property
    def reactions(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._format == Format.HDF5:
            return self._read_hdf5("exyz")
        if self._format == Format.CSTRUCT:
            return self._from_struct_("E", shape=(3,), units="V/m")
        return self._exyz.copy()

    def export(self):
        data = {
            "bxyz":self.bxyz,
        }
        return data

    def stage(self):
        init = LIBASCOT.E_TC_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.exyz,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._exyz
            self._staged = True

    def unstage(self):
        free = LIBASCOT.E_TC_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._exyz = self.exyz
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateAsigmaLocMixin():
    """Mixin class used by `Data` to create Asimga_loc input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_atomicreactions(
            self,
            reactions: utils.ArrayLike | None = None,
            chargenumber: utils.ArrayLike | None = None,
            massnumber: utils.ArrayLike | None = None,
            energygrid: Tuple[utils.ArrayLike] | None = None,
            densitygrid: Tuple[utils.ArrayLike] | None = None,
            temperaturegrid: Tuple[utils.ArrayLike] | None = None,
            crosssection: utils.ArrayLike | None = None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> Asigma_loc:
        r"""Atomic reaction data stored locally.

        Parameters
        ----------
        reactions : array_like (nreaction, 1)
            Type of atomic reaction.
        chargenumber : array_like (nreaction, 2)
            Atomic numbers of test particles and the bulk particles.
        massnumber : array_like (nreaction, 2)
            Atomic mass number of test particle.
        energygrid : array_like (nreaction,)
            Energy grids for each reaction.
        densitygrid : array_like (nreaction,)
            Density grids for each reaction.
        temperaturegrid : array_like (nreaction,)
            Temperature grids for each reaction.
        crosssection : array_like (nreaction,)
            Reaction cross-section data [cm^2].
        """


class Atomic(ctypes.Structure):
    """Wrapper for the atomic data in libascot.so."""

    _fields_ = [
        ('asigma_loc', ctypes.POINTER(Struct)),
        ('type', ctypes.c_uint32),
        ]
