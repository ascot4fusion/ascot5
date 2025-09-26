"""Defines Plasma1D radial plasma input class and the corresponding factory
method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.physlib import Species
from a5py.libascot import LIBASCOT, DataStruct, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in plasma_1D.h."""

    _fields_ = [
        ("n_rho", ctypes.c_int32),
        ("n_species", ctypes.c_int32),
        ("mass", ctypes.POINTER(ctypes.c_double)),
        ("charge", ctypes.POINTER(ctypes.c_double)),
        ("anum", ctypes.POINTER(ctypes.c_int32)),
        ("znum", ctypes.POINTER(ctypes.c_int32)),
        ("rho", ctypes.POINTER(ctypes.c_double)),
        ("temp", ctypes.POINTER(ctypes.c_double)),
        ("dens", ctypes.POINTER(ctypes.c_double)),
        ("vtor", ctypes.POINTER(ctypes.c_double)),
        ]


init_fun(
    "plasma_1D_init",
    ctypes.POINTER(Struct),
    *(2*[ctypes.c_int32]),
    ndpointer(ctypes.c_double),
    *(2*[ndpointer(ctypes.c_int32)]),
    *(7*[ndpointer(ctypes.c_double)]),
    restype=ctypes.c_int32,
    )

init_fun("plasma_1D_free", ctypes.POINTER(Struct))


@Leaf.register
class Plasma1D(InputVariant):
    """Radial plasma profile."""

    _cdata: Optional[Struct]

    @property
    def nrho(self) -> int:
        """Number of radial grid points."""
        if self._cdata is not None:
            return int(self._cdata.readonly_carray("n_rho", ()))
        assert self._file is not None
        return self._file.read("rhogrid").size

    @property
    def nion(self) -> int:
        """Number of ion species."""
        if self._cdata is not None:
            return int(self._cdata.readonly_carray("n_species", ()) - 1)
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
            anum = self._cdata.readonly_carray("anum", (self.nion,))
            znum = self._cdata.readonly_carray("znum", (self.nion,))
        else:
            assert self._file is not None
            anum, znum = self._file.read("anum"), self._file.read("znum")
        return [Species.from_znumanum(z, a) for a, z in zip(anum, znum)]

    @property
    def rhogrid(self) -> unyt.unyt_array:
        r"""Radial grid in :math:`\rho` which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("rho", (self.nrho,), "1")
        assert self._file is not None
        return self._file.read("rhogrid")

    @property
    def ni(self) -> unyt.unyt_array:
        """Density for each ion species."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray(
                "dens", (self.nion, self.nrho), "m**(-3)",
                )
            return data.T[:,1:]
        assert self._file is not None
        return self._file.read("ni")

    @property
    def Ti(self) -> unyt.unyt_array:
        """Ion temperature."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray("temp", (self.nrho, 2), "J")
            return data[:,0].to("eV")
        assert self._file is not None
        return self._file.read("Ti")

    @property
    def ne(self) -> unyt.unyt_array:
        """Electron density."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray(
                "dens", (self.nion, self.nrho), "m**(-3)",
                )
            return data.T[:,0]
        assert self._file is not None
        return self._file.read("ne")

    @property
    def Te(self) -> unyt.unyt_array:
        """Electron temperature."""
        if self._cdata is not None:
            nrho = self.rhogrid.size
            data = self._cdata.readonly_carray("temp", (nrho, 2), "J")
            return data[:,1].to("eV")
        assert self._file is not None
        return self._file.read("Te")

    @property
    def charge(self) -> unyt.unyt_array:
        """Ion charge states."""
        if self._cdata is not None:
            return self._cdata.readonly_carray(
                "charge", (self.nion,), "C"
                ).to("e")
        assert self._file is not None
        return self._file.read("charge")

    @property
    def rotation(self) -> unyt.unyt_array:
        """Toroidal rotation of the plasma."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("vtor", (self.nrho,), "rad/s")
        assert self._file is not None
        return self._file.read("rotation")

    #pylint: disable=too-many-arguments
    def _stage(
            self, species: list[Species],
            rhogrid: unyt.unyt_array,
            ni: unyt.unyt_array,
            Ti: unyt.unyt_array,
            ne: unyt.unyt_array,
            Te: unyt.unyt_array,
            charge: unyt.unyt_array,
            rotation: unyt.unyt_array,
            ) -> None:
        anum = np.array([s.anum for s in species], dtype="i4")
        znum = np.array([s.znum for s in species], dtype="i4")
        mass = unyt.unyt_array([s.mass for s in species], dtype="f8")
        self._cdata = Struct()
        if LIBASCOT.plasma_1D_init(
            ctypes.byref(self._cdata), rhogrid.size, len(species), rhogrid.v,
            anum, znum, mass.to("kg").v, charge.to("C").v.astype("f8"),
            Te.to("J").v, Ti.to("J").v, ni.v, ne.v, rotation.v,
            ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "rhogrid", "ni", "Ti", "ne", "Te", "charge", "rotation",
            ]:
            self._file.write(field, getattr(self, field))

        self._file.write("anum", self.anum)
        self._file.write("znum", self.znum)

    def export(self) -> dict[str, unyt.unyt_array | list[Species]]:
        fields = [
            "rhogrid", "ni", "Ti", "ne", "Te", "charge", "rotation", "species",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            species=self.species, rhogrid=self.rhogrid, ni=self.ni, Ti=self.Ti,
            ne=self.ne, Te=self.Te, charge=self.charge, rotation=self.rotation,
            )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.plasma_1D_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create Plasma1D input."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_plasma1d(
            self,
            species: list[str] | list[Species],
            rhogrid: unyt.unyt_array,
            ni: unyt.unyt_array,
            Ti: unyt.unyt_array,
            ne: Optional[unyt.unyt_array]=None,
            Te: Optional[unyt.unyt_array]=None,
            charge: Optional[unyt.unyt_array]=None,
            rotation: Optional[unyt.unyt_array]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> Plasma1D:
        r"""Create plasma profiles that have only radial dependency.

        This is the most usual plasma input. The profiles should extend beyond
        separatrix since they are not automatically extrapolated and marker will
        be terminated if it enters a region where plasma data is not defined.

        The data is interpolated linearly.

        Parameters
        ----------
        species : list[str] *or* list[:class:`.Species`] (nion,)
            The ion species that make up the plasma.
        rhogrid : array_like (nrho,)
            Radial grid in :math:`\rho` in which the data is tabulated.

            This grid doesn't have to be uniform.
        ni : array_like (nrho,nion)
            Density for each ion species.
        Ti : array_like (nrho,)
            Ion temperature.
        ne : array_like (nrho,), *optional*
            Electron density.

            By default, the electron density is determined from the ion charge
            density so that the plasma is quasi-neutral.
        Te : array_like (nrho,), *optional*
            Electron temperature.

            Same as ion temperature by default.
        charge : array_like (nion,), *optional*
            Ion charge states.

            Ions are fully ionized by default.
        rotation : array_like (nrho,), *optional*
            Toroidal rotation of the plasma.

            Zero by default.
        note : str, *optional*
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, *optional*
            Set this input as active on creation.
        preview : bool, *optional*
            If True, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : :class:`.Plasma1D`
            Input variant created from the given parameters.
        """
        species = [s if isinstance(s, Species) else Species.from_string(s)
                   for s in species]
        nion = len(species)
        znum = np.array([s.znum for s in species])

        with utils.validate_variables() as v:
            rhogrid = v.validate("rhogrid", rhogrid, (-1,), "1")

        nrho = rhogrid.size
        with utils.validate_variables() as v:
            ni = v.validate("ni", ni, (nrho, nion), "m**(-3)")
            Ti = v.validate("Ti", Ti, (nrho,), "eV")
            charge = v.validate("charge", charge, (nion,), "e", default=znum)
            rotation = v.validate(
                "rotation", rotation, (nrho,), "rad/s", default=np.full(nrho, 0))

        if charge is None:
            charge_density = np.matmul(ni, znum)
        else:
            charge_density = np.matmul(ni, charge) / unyt.e

        with utils.validate_variables() as v:
            ne = v.validate("ne", ne, (nrho,), "m**(-3)", default=charge_density)
            Te = v.validate("Te", Te, (nrho,), "eV", default=Ti.v)

        leaf = Plasma1D(note=note)
        leaf._stage(
            species=species, rhogrid=rhogrid, ni=ni, Ti=Ti, ne=ne, Te=Te,
            charge=charge, rotation=rotation,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="plasma",
            )
        return leaf
