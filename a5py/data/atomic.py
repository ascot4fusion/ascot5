"""Defines the atomic reaction input class and the corresponding factory method."""

import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from a5py import utils
from a5py.physlib import Species
from a5py.libascot import (
    LIBASCOT,
    DataStruct,
    Linear1D,
    Linear2D,
    Linear3D,
    init_fun,
)
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in asigma_loc.h."""

    class AtomicData(DataStruct):

        class LinearUnion(ctypes.Union):

            _fields_ = [
                ("linear1d", Linear1D),
                ("linear2d", Linear2D),
                ("linear3d", Linear3D),
            ]

        _anonymous_ = ("linear",)
        _fields_ = [
            ("logscale", ctypes.c_int32 * 3),
            ("linear", LinearUnion),
        ]

    _fields_ = [
        ("EI", ctypes.POINTER(AtomicData)),
        ("II", ctypes.POINTER(AtomicData)),
        ("CX", ctypes.POINTER(AtomicData)),
        ("BMS", ctypes.POINTER(AtomicData)),
    ]


class _DataTemplate:
    """Base class for atomic data classes."""

    def __init__(
        self,
        beam_species: str | Species,
        beam_charge: int,
        target_species: str | Species,
        target_charge: int,
        sigmav,
        temperature,
        energy=None,
        density=None,
    ):

        def isloggrid(abscissa, name):
            """Check if abscissa is logarithmic or linear grid."""
            try:
                utils.validate_abscissa(abscissa, name, uniform=True)
                return False
            except ValueError:
                try:
                    utils.validate_abscissa(np.log10(abscissa), name, uniform=True)
                    return True
                except ValueError as e:
                    raise e from None

        units = 1
        if hasattr(sigmav, "units"):
            units = sigmav.units
        sigmav = np.asarray(sigmav) * units
        datashape = sigmav.shape
        if energy is not None and density is not None:
            if sigmav.ndim != 3:
                raise ValueError(f"Expected sigmav.ndim == 3, got {sigmav.ndim}.")
        elif energy is not None:
            if sigmav.ndim != 2:
                raise ValueError(f"Expected sigmav.ndim == 2, got {sigmav.ndim}.")
        else:
            if sigmav.ndim != 1:
                raise ValueError(f"Expected sigmav.ndim == 1, got {sigmav.ndim}.")

        with utils.validate_variables() as v:
            temperature = v.validate("temperature", temperature, (datashape[0],), "eV")
            if energy is not None:
                energy = v.validate("energy", energy, (datashape[1],), "eV")
            if density is not None:
                density = v.validate("density", density, (datashape[2],), "m**-3")
            sigmav = v.validate("sigmav", sigmav, datashape, "m**3/s")

        self.beam_species = (
            beam_species
            if isinstance(beam_species, Species)
            else Species.from_string(beam_species)
        )
        self.target_species = (
            target_species
            if isinstance(target_species, Species)
            else Species.from_string(target_species)
        )
        if not isinstance(target_charge, int):
            raise ValueError("`target_charge` must be an integer.")
        if not isinstance(beam_charge, int):
            raise ValueError("`beam_charge` must be an integer.")
        self.beam_charge = beam_charge
        self.target_charge = target_charge

        self.temperature = temperature
        self.temperature_log = isloggrid(temperature, "temperature")
        if energy is not None:
            self.energy = energy
            self.energy_log = isloggrid(energy, "energy")
        if density is not None:
            self.density = density
            self.density_log = isloggrid(density, "density")
        self.sigmav = sigmav


class EIdata(_DataTemplate):
    """Electron impact ionization rate coefficient."""

    def __init__(self, beam_species, beam_charge, temperature, sigmav):
        super().__init__(
            beam_species, beam_charge, "electron", -1, sigmav, temperature, None, None
        )


class IIdata(_DataTemplate):
    """Ion impact ionization rate coefficient."""

    def __init__(
        self,
        beam_species,
        beam_charge,
        target_species,
        target_charge,
        temperature,
        energy,
        sigmav,
    ):
        super().__init__(
            beam_species,
            beam_charge,
            target_species,
            target_charge,
            sigmav,
            temperature,
            energy,
            None,
        )


class CXdata(_DataTemplate):
    """Charge exchange rate coefficient."""

    def __init__(
        self, beam_species, beam_charge, target_species, temperature, energy, sigmav
    ):
        super().__init__(
            beam_species,
            beam_charge,
            target_species,
            0,
            sigmav,
            temperature,
            energy,
            None,
        )


class BMSdata(_DataTemplate):
    """Beam stopping rate coefficient."""

    def __init__(
        self,
        beam_species,
        target_species,
        target_charge,
        temperature,
        energy,
        density,
        sigmav,
    ):
        super().__init__(
            beam_species,
            0,
            target_species,
            target_charge,
            sigmav,
            temperature,
            energy,
            density,
        )


@Leaf.register
class Atomic(InputVariant):
    """Local atomic data."""

    @property
    def EI(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._format == Format.HDF5:
            return self._read_hdf5("exyz")
        if self._format == Format.CSTRUCT:
            return self._from_struct_("E", shape=(3,), units="V/m")
        return self._exyz.copy()

    @property
    def II(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._format == Format.HDF5:
            return self._read_hdf5("exyz")
        if self._format == Format.CSTRUCT:
            return self._from_struct_("E", shape=(3,), units="V/m")
        return self._exyz.copy()

    @property
    def CX(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._format == Format.HDF5:
            return self._read_hdf5("exyz")
        if self._format == Format.CSTRUCT:
            return self._from_struct_("E", shape=(3,), units="V/m")
        return self._exyz.copy()

    @property
    def BMS(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._format == Format.HDF5:
            return self._read_hdf5("exyz")
        if self._format == Format.CSTRUCT:
            return self._from_struct_("E", shape=(3,), units="V/m")
        return self._exyz.copy()

    def export(self):
        data = {
            "bxyz": self.bxyz,
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
class CreateAtomicMixin:
    """Mixin class used by `Data` to create Asimga_loc input."""

    # pylint: disable=protected-access, too-many-arguments
    def create_atomic(
        self,
        CX: list[CXdata],
        BMS: list[BMSdata],
        II: list[IIdata],
        EI: Optional[EIdata] = None,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
    ) -> Atomic:
        r"""Atomic reaction data stored locally.

        Parameters
        ----------
        reactions : array_like (nreaction, 1)
            Type of atomic reaction.
        crosssection : array_like (nreaction,)
            Reaction cross-section data [cm^2].
        """
        leaf = Atomic(note=note)
        leaf.beam = []
        leaf.target = []
        leaf._cdata = Struct()

        leaf._cdata.BMS = (Struct.AtomicData * len(BMS))()
        if not len(BMS):
            leaf._cdata.BMS = None
        for i, reaction in enumerate(BMS):
            ngrid = reaction.density.to("m**-3").v
            Tgrid = reaction.temperature.to("J").v
            egrid = reaction.energy.to("J").v

            s = leaf._cdata.BMS[i]
            s.logscale[0] = reaction.energy_log
            s.logscale[1] = reaction.temperature_log
            s.logscale[2] = reaction.density_log
            s.ngrid = ngrid if not reaction.density_log else np.log10(ngrid)
            s.Tgrid = Tgrid if not reaction.temperature_log else np.log10(Tgrid)
            s.egrid = egrid if not reaction.energy_log else np.log10(egrid)
            s.sigmav = np.ascontiguousarray(np.log(reaction.sigmav.to("m**3/s").v))
            LIBASCOT.Linear3D_init(
                ctypes.byref(s.linear3d),
                s.egrid.size,
                s.Tgrid.size,
                s.ngrid.size,
                0,
                0,
                0,
                s.egrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.Tgrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.ngrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.sigmav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )

        leaf._cdata.CX = (Struct.AtomicData * len(CX))()
        if not len(CX):
            leaf._cdata.CX = None
        for i, reaction in enumerate(CX):
            Tgrid = reaction.temperature.to("J").v
            egrid = reaction.energy.to("J").v

            s = leaf._cdata.CX[i]
            s.logscale[0] = reaction.energy_log
            s.logscale[1] = reaction.temperature_log
            s.Tgrid = Tgrid if not reaction.temperature_log else np.log10(Tgrid)
            s.egrid = egrid if not reaction.energy_log else np.log10(egrid)
            s.sigmav = np.ascontiguousarray(np.log(reaction.sigmav.to("m**3/s").v))
            LIBASCOT.Linear2D_init(
                ctypes.byref(s.linear2d),
                s.egrid.size,
                s.Tgrid.size,
                0,
                0,
                s.egrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.Tgrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.sigmav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )

        leaf._cdata.II = (Struct.AtomicData * len(II))()
        if not len(II):
            leaf._cdata.II = None
        for i, reaction in enumerate(II):
            Tgrid = reaction.temperature.to("J").v
            egrid = reaction.energy.to("J").v

            s = leaf._cdata.II[i]
            s.logscale[0] = reaction.energy_log
            s.logscale[1] = reaction.temperature_log
            s.Tgrid = Tgrid if not reaction.temperature_log else np.log10(Tgrid)
            s.egrid = egrid if not reaction.energy_log else np.log10(egrid)
            s.sigmav = np.ascontiguousarray(np.log(reaction.sigmav.to("m**3/s").v))
            LIBASCOT.Linear2D_init(
                ctypes.byref(s.linear2d),
                s.egrid.size,
                s.Tgrid.size,
                0,
                0,
                s.egrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.Tgrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.sigmav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )

        if EI is None:
            leaf._cdata.EI = None
        else:
            leaf._cdata.EI = (Struct.AtomicData * 1)()
            Tgrid = EI.temperature.to("J").v

            s = leaf._cdata.EI[0]
            s.logscale[0] = EI.temperature_log
            s.Tgrid = Tgrid if not EI.temperature_log else np.log10(Tgrid)
            s.sigmav = np.ascontiguousarray(np.log(EI.sigmav.to("m**3/s").v))
            LIBASCOT.Linear1D_init(
                ctypes.byref(s.linear1d),
                s.Tgrid.size,
                0,
                s.Tgrid[[0, -1]].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s.sigmav.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )

        return leaf
