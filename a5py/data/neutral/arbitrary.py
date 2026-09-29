"""Defines Neutral3D arbitrary neutral density input class and the corresponding
factory method.
"""

import ctypes
from typing import Tuple, List, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.physlib import Species
from a5py.libascot import LIBASCOT, DataStruct, Linear3D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in N0_3D.h."""

    _fields_ = [
        ("nspecies", ctypes.c_int32),
        ("anum", ctypes.POINTER(ctypes.c_int32)),
        ("znum", ctypes.POINTER(ctypes.c_int32)),
        ("n", ctypes.POINTER(Linear3D)),
        ("T", ctypes.POINTER(Linear3D)),
    ]


init_fun(
    "NeutralArbitrary_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    *(5*[ndpointer(ctypes.c_double)]),
)

init_fun("NeutralArbitrary_free", ctypes.POINTER(Struct))


@Leaf.register
class NeutralArbitrary(InputVariant):
    """Arbitrary neutral density in 3D."""

    @property
    def nspecies(self) -> int:
        r"""Number of ion species."""
        if self._cdata is not None:
            return int(self._cdata.readonly_carray("nspecies", ()) - 1)
        assert self._file is not None
        return self._file.read("znum").size

    @property
    def anum(self) -> np.ndarray:
        r"""Atomic mass number of each ion species."""
        return np.array([s.anum for s in self.species], dtype="i4")

    @property
    def znum(self) -> np.ndarray:
        r"""Atomic number of each ion species."""
        return np.array([s.znum for s in self.species], dtype="i4")

    @property
    def mass(self) -> unyt.unyt_array:
        r"""Mass of each ion species."""
        return unyt.unyt_array([s.mass for s in self.species], dtype="f8")

    @property
    def species(self) -> list[Species]:
        r"""The ion species that make up the plasma."""
        if self._cdata is not None:
            anum = self._cdata.readonly_carray("anum", (self.nspecies,))
            znum = self._cdata.readonly_carray("znum", (self.nspecies,))
        else:
            assert self._file is not None
            anum, znum = self._file.read("anum"), self._file.read("znum")
        return [Species.from_znumanum(z, a) for a, z in zip(anum, znum)]

    @property
    def rgrid(self) -> unyt.unyt_array:
        r"""Radial grid in :math:`R` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "m", "density")
        assert self._file is not None
        return self._file.read("rgrid")

    @property
    def zgrid(self) -> unyt.unyt_array:
        r"""Axial grid in :math:`z` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("z", "m", "density")
        assert self._file is not None
        return self._file.read("zgrid")

    @property
    def phigrid(self) -> unyt.unyt_array:
        r"""Toroidal grid in :math:`\phi` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "rad", "density").to("deg")
        assert self._file is not None
        return self._file.read("phigrid")

    @property
    def temperature(self):
        r"""Species-wise temperature."""
        if self._staged:
            nspecies = self._from_struct_("n_species", shape=())
            data = self._from_struct_("t0", idx=0)
            for i in range(1, nspecies):
                data = np.stack((data, self._from_struct_("t0", idx=i)), axis=3)
            if nspecies == 1:
                data = np.expand_dims(data, axis=3)
            return (data * unyt.J).to("eV")
        if self._format == Format.HDF5:
            return self._read_hdf5("temperature")
        return self._temperature.copy()

    @property
    def density(self):
        r"""Species-wise density."""
        if self._staged:
            nspecies = self._from_struct_("n_species", shape=())
            data = self._from_struct_("n0", idx=0)
            for i in range(1, nspecies):
                data = np.stack((data, self._from_struct_("n0", idx=i)), axis=3)
            if nspecies == 1:
                data = np.expand_dims(data, axis=3)
            return data * unyt.m ** (-3)
        if self._format == Format.HDF5:
            return self._read_hdf5("density")
        return self._density.copy()

    def _stage(
        self,
        species: list[Species],
        rgrid: unyt.unyt_array,
        zgrid: unyt.unyt_array,
        phigrid: unyt.unyt_array,
        density: unyt.unyt_array,
        temperature: unyt.unyt_array,
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.NeutralArbitrary_init(
            ctypes.byref(self._cdata),
            len(species),
            rgrid.size,
            phigrid.size,
            zgrid.size,
            rgrid[[0, -1]],
            phigrid[[0, -1]],
            zgrid[[0, -1]],
            density,
            temperature,
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for attr in ["rgrid", "phigrid", "zgrid", "density", "temperature"]:
            self._file.write(attr, getattr(self, attr))

    def export(self) -> dict[str, unyt.unyt_array]:
        data = {
            "rgrid": self.rgrid,
            "phigrid": self.phigrid,
            "zgrid": self.zgrid,
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
        LIBASCOT.NeutralArbitrary_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create Neutral3D input."""

    # pylint: disable=protected-access, too-many-arguments
    def create_neutralarbitrary(
        self,
        rgrid: utils.ArrayLike,
        phigrid: utils.ArrayLike,
        zgrid: utils.ArrayLike,
        species: List[str] | Tuple[str],
        density: utils.ArrayLike,
        temperature: utils.ArrayLike,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
    ) -> NeutralArbitrary:
        r"""Create arbitrary 3D neutral density input.

        The data is interpolated linearly.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in R in which data is tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in phi in which data is tabulated.
        zgrid : array_like (nz,)
            The uniform grid in z in which data is tabulated.
        species : list[str] or tuple[str]
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
        inputdata : ~a5py.data.neutral.Neutral3D
            Input variant created from the given parameters.
        """
        with utils.validate_variables() as v:
            rgrid = v.validate("rgrid", rgrid, (-1,), "m")
            zgrid = v.validate("zgrid", zgrid, (-1,), "m")
            phigrid = v.validate("phigrid", phigrid, (-1,), "deg")

        nr, nz, nphi, nspecies = rgrid.size, zgrid.size, phigrid.size, len(species)
        with utils.validate_variables() as v:
            density = v.validate(
                "density", density, (nr, nphi, nz, nspecies), "m**(-3)"
            )
            temperature = v.validate(
                "temperature", temperature, (nr, nphi, nz, nspecies), "eV"
            )

        utils.validate_abscissa(rgrid, "rgrid")
        utils.validate_abscissa(zgrid, "zgrid")
        utils.validate_abscissa(phigrid, "phigrid", periodic=True)
        leaf = NeutralArbitrary(note=note)
        leaf._stage(
            rgrid=rgrid,
            phigrid=phigrid,
            zgrid=zgrid,
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
