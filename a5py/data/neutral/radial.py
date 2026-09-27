"""Defines Neutral1D radial neutral density input class and the corresponding
factory method.
"""

import ctypes
from typing import Tuple, List, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.physlib import Species
from a5py.libascot import LIBASCOT, DataStruct, Linear1D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in N0_1D.h."""

    _fields_ = [
        ("nspecies", ctypes.c_int32),
        ("anum", ctypes.POINTER(ctypes.c_int32)),
        ("znum", ctypes.POINTER(ctypes.c_int32)),
        ("n", ctypes.POINTER(Linear1D)),
        ("T", ctypes.POINTER(Linear1D)),
    ]


init_fun(
    "NeutralRadial_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
    ctypes.c_size_t,
    *(3*[ndpointer(ctypes.c_double)]),
)

init_fun("NeutralRadial_free", ctypes.POINTER(Struct))


@Leaf.register
class NeutralRadial(InputVariant):
    """Radial neutral profile."""

    @property
    def nspecies(self) -> int:
        """Number of ion species."""
        if self._cdata is not None:
            return int(self._cdata.readonly_carray("nspecies", ()) - 1)
        assert self._file is not None
        return self._file.read("znum").size

    @property
    def anum(self) -> np.ndarray:
        """Atomic mass number of each ion species."""
        return np.array([s.anum for s in self.species], dtype="i4")

    @property
    def znum(self) -> np.ndarray:
        """Atomic number of each ion species."""
        return np.array([s.znum for s in self.species], dtype="i4")

    @property
    def mass(self) -> unyt.unyt_array:
        """Mass of each ion species."""
        return unyt.unyt_array([s.mass for s in self.species], dtype="f8")

    @property
    def species(self) -> list[Species]:
        """The ion species that make up the plasma."""
        if self._cdata is not None:
            anum = self._cdata.readonly_carray("anum", (self.nspecies,))
            znum = self._cdata.readonly_carray("znum", (self.nspecies,))
        else:
            assert self._file is not None
            anum, znum = self._file.read("anum"), self._file.read("znum")
        return [Species.from_znumanum(z, a) for a, z in zip(anum, znum)]

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in :math:`\rho` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "1", "density")
        assert self._file is not None
        return self._file.read("rhogrid")

    @property
    def temperature(self):
        """Species-wise temperature."""
        if self._staged:
            nspecies = self._from_struct_("n_species", shape=())
            data = self._from_struct_("t0", idx=0)
            for i in range(1, nspecies):
                data = np.stack((data, self._from_struct_("t0", idx=i)), axis=1)
            if nspecies == 1:
                data = np.expand_dims(data, axis=1)
            return (data * unyt.J).to("eV")
        if self._format == Format.HDF5:
            return self._read_hdf5("temperature")
        return self._temperature.copy()

    @property
    def density(self):
        """Species-wise density."""
        if self._staged:
            nspecies = self._from_struct_("n_species", shape=())
            data = self._from_struct_("n0", idx=0)
            for i in range(1, nspecies):
                data = np.stack((data, self._from_struct_("n0", idx=i)), axis=1)
            if nspecies == 1:
                data = np.expand_dims(data, axis=1)
            return data * unyt.m ** (-3)
        if self._format == Format.HDF5:
            return self._read_hdf5("density")
        return self._density.copy()

    def _stage(
        self,
        species: list[Species],
        rhogrid: unyt.unyt_array,
        density: unyt.unyt_array,
        temperature: unyt.unyt_array,
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.NeutralRadial_init(
            ctypes.byref(self._cdata),
            len(species),
            rhogrid.size,
            rhogrid[[0, -1]].v,
            density.v,
            temperature.v,
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("rhogrid", self.rhogrid)
        self._file.write("density", self.density)
        self._file.write("temperature", self.temperature)

    def export(self) -> dict[str, unyt.unyt_array]:
        data = {
            "species": self.species,
            "rhogrid": self.rhogrid,
            "density": self.density,
            "temperature": self.temperature,
        }
        return data

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.NeutralRadial_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create `NeutralRadial` input."""

    # pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_neutralradial(
        self,
        rhogrid: utils.ArrayLike,
        species: List[str] | Tuple[str],
        density: utils.ArrayLike,
        temperature: utils.ArrayLike,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
    ) -> NeutralRadial:
        r"""Create radial neutral density.

        The data is interpolated linearly.

        Parameters
        ----------
        rhogrid : array_like (nrho,1)
            Uniform radial grid in rho in which the data is tabulated.
        species : list[str] or tuple[str] (nspecies,)
            Name(s) of the neutral species.
        density : array_like (nrho,nspecies)
            Density for each species.
        temperature : array_like (nrho,nspecies)
            Temperature of each species.
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        preview : bool, *optional*
            If True, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : ~a5py.data.neutral.Neutral1D
            Input variant created from the given parameters.
        """
        with utils.validate_variables() as v:
            rhogrid = v.validate("rhogrid", rhogrid, (-1,), "1")

        species = [
            s if isinstance(s, Species) else Species.from_string(s) for s in species
        ]

        nrho, nspecies = rhogrid.size, len(species)
        with utils.validate_variables() as v:
            density = v.validate("density", density, (nrho, nspecies), "m**(-3)")
            temperature = v.validate("temperature", temperature, (nrho, nspecies), "eV")

        utils.validate_abscissa(rhogrid, "rhogrid")
        leaf = NeutralRadial(note=note)
        leaf._stage(
            rhogrid=rhogrid,
            density=density,
            temperature=temperature,
            species=species,
        )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf,
            activate=activate,
            save=save,
            category="neutral",
        )
        return leaf
