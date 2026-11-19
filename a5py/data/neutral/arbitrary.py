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
        ("density", ctypes.POINTER(Linear3D)),
        ("temperature", ctypes.POINTER(Linear3D)),
    ]


init_fun(
    "NeutralArbitrary_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
    *(2 * [ndpointer(ctypes.c_int32)]),
    *(5 * [ndpointer(ctypes.c_double)]),
)

init_fun("NeutralArbitrary_free", ctypes.POINTER(Struct))


@Leaf.register
class NeutralArbitrary(InputVariant):
    """Arbitrary neutral density in 3D."""

    @property
    def nspecies(self) -> int:
        r"""Number of neutral species."""
        if self._cdata is not None:
            return int(self._cdata.readonly_carray("nspecies", ()))
        assert self._file is not None
        return self._file.read("znum").size

    @property
    def anum(self) -> np.ndarray:
        r"""Atomic mass number of each neutral species."""
        return np.array([s.anum for s in self.species], dtype="i4")

    @property
    def znum(self) -> np.ndarray:
        r"""Atomic number of each neutral species."""
        return np.array([s.znum for s in self.species], dtype="i4")

    @property
    def mass(self) -> unyt.unyt_array:
        r"""Mass of each neutral species."""
        return unyt.unyt_array([s.mass for s in self.species], dtype="f8")

    @property
    def species(self) -> tuple[Species]:
        r"""The ion species that make up the plasma."""
        if self._cdata is not None:
            anum = self._cdata.readonly_carray("anum", (self.nspecies,))
            znum = self._cdata.readonly_carray("znum", (self.nspecies,))
        else:
            assert self._file is not None
            anum, znum = self._file.read("anum"), self._file.read("znum")
        return tuple(Species.from_znumanum(z, a) for a, z in zip(anum, znum))

    @property
    def rgrid(self) -> unyt.unyt_array:
        r"""Radial grid in :math:`R` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "m", "density", 0)
        assert self._file is not None
        return self._file.read("rgrid")

    @property
    def zgrid(self) -> unyt.unyt_array:
        r"""Axial grid in :math:`z` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("z", "m", "density", 0)
        assert self._file is not None
        return self._file.read("zgrid")

    @property
    def phigrid(self) -> unyt.unyt_array:
        r"""Toroidal grid in :math:`\phi` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "rad", "density", 0, convert="deg")
        assert self._file is not None
        return self._file.read("phigrid")

    @property
    def T0(self) -> unyt.unyt_array:
        r"""Temperature of neutral species."""
        if self._cdata is not None:
            nr, nphi, nz = self.rgrid.size, self.phigrid.size-1, self.zgrid.size
            nspecies = self.nspecies
            data = np.zeros((nr, nphi, nz, nspecies), dtype="f8") * unyt.eV
            for i in range(self.nspecies):
                data[:,:,:,i] = self._cdata.readonly_interp("temperature", "eV", i)
            data.flags.writeable = False
            return data
        assert self._file is not None
        return self._file.read("T0")

    @property
    def n0(self) -> unyt.unyt_array:
        r"""Density of neutral species."""
        if self._cdata is not None:
            nr, nphi, nz = self.rgrid.size, self.phigrid.size-1, self.zgrid.size
            nspecies = self.nspecies
            data = np.zeros((nr, nphi, nz, nspecies), dtype="f8") * unyt.m**-3
            for i in range(self.nspecies):
                data[:,:,:,i] = self._cdata.readonly_interp("density", "m**(-3)", i)
            data.flags.writeable = False
            return data
        assert self._file is not None
        return self._file.read("n0")

    def _stage(
        self,
        species: list[Species],
        rgrid: unyt.unyt_array,
        zgrid: unyt.unyt_array,
        phigrid: unyt.unyt_array,
        n0: unyt.unyt_array,
        T0: unyt.unyt_array,
    ) -> None:
        anum = np.array([s.anum for s in species], dtype="i4")
        znum = np.array([s.znum for s in species], dtype="i4")
        self._cdata = Struct()
        if LIBASCOT.NeutralArbitrary_init(
            ctypes.byref(self._cdata),
            len(species),
            rgrid.size,
            phigrid.size-1,
            zgrid.size,
            anum,
            znum,
            rgrid[[0, -1]],
            phigrid[[0, -1]].to("rad"),
            zgrid[[0, -1]],
            np.ascontiguousarray(n0.transpose(3, 0, 1, 2)),
            np.ascontiguousarray(T0.transpose(3, 0, 1, 2)),
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for attr in ["anum", "znum", "rgrid", "phigrid", "zgrid", "n0", "T0"]:
            self._file.write(attr, getattr(self, attr))

    def export(self) -> dict[str, unyt.unyt_array]:
        data = {
            "species": list(self.species),
            "rgrid": self.rgrid.copy(),
            "phigrid": self.phigrid.copy(),
            "zgrid": self.zgrid.copy(),
            "n0": self.n0.copy(),
            "T0": self.T0.copy(),
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
        species: List[str] | Tuple[str],
        rgrid: utils.ArrayLike,
        phigrid: utils.ArrayLike,
        zgrid: utils.ArrayLike,
        n0: utils.ArrayLike,
        T0: utils.ArrayLike,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
    ) -> NeutralArbitrary:
        r"""Create arbitrary 3D neutral density input.

        The data is interpolated linearly.

        Parameters
        ----------
        species : list[str] or tuple[str]
            Neutral species.
        rgrid : array_like (nr,)
            The uniform grid in R in which data is tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in phi in which data is tabulated.
        zgrid : array_like (nz,)
            The uniform grid in z in which data is tabulated.
        n0 : array_like (nr,nphi-1,nz,nspecies)
            Density for each species.
        T0 : array_like (nr,nphi-1,nz,nspecies)
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

        species = [
            s if isinstance(s, Species) else Species.from_string(s) for s in species
        ]
        if len(species) == 0:
            raise ValueError("At least one species must be specified.")

        nr, nz, nphi, nspecies = rgrid.size, zgrid.size, phigrid.size, len(species)
        with utils.validate_variables() as v:
            n0 = v.validate(
                "n0", n0, (nr, nphi-1, nz, nspecies), "m**(-3)"
            )
            T0 = v.validate(
                "T0", T0, (nr, nphi-1, nz, nspecies), "eV"
            )

        utils.validate_abscissa(rgrid, "rgrid")
        utils.validate_abscissa(zgrid, "zgrid")
        utils.validate_abscissa(phigrid, "phigrid", periodic=True)
        leaf = NeutralArbitrary(note=note)
        leaf._stage(
            rgrid=rgrid,
            phigrid=phigrid,
            zgrid=zgrid,
            n0=n0,
            T0=T0,
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
