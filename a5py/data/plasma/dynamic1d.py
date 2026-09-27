"""Defines Plasma1DDynamic time-dependent radial plasma input class and
the corresponding factory method.
"""

import ctypes
from typing import Tuple, List, Optional

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
    """Python wrapper for the struct in plasma_1Dt.h."""

    _fields_ = [
        ("nrho", ctypes.c_int32),
        ("ntime", ctypes.c_int32),
        ("nspecies", ctypes.c_int32),
        ("anum", ctypes.POINTER(ctypes.c_int32)),
        ("znum", ctypes.POINTER(ctypes.c_int32)),
        ("mass", ctypes.POINTER(ctypes.c_double)),
        ("charge", ctypes.POINTER(ctypes.c_double)),
        ("rho", ctypes.POINTER(ctypes.c_double)),
        ("time", ctypes.POINTER(ctypes.c_double)),
        ("temp", ctypes.POINTER(ctypes.c_double)),
        ("dens", ctypes.POINTER(ctypes.c_double)),
        ("vtor", ctypes.POINTER(ctypes.c_double)),
    ]


init_fun(
    "PlasmaDynamic1D_init",
    ctypes.POINTER(Struct),
    *(3 * [ctypes.c_size_t]),
    *(2 * [ndpointer(ctypes.c_int32)]),
    *(9 * [ndpointer(ctypes.c_double)]),
)

init_fun("PlasmaDynamic1D_free", ctypes.POINTER(Struct))


@Leaf.register
class PlasmaDynamic1D(InputVariant):
    """Time-dependent radial plasma profile."""

    @property
    def rhogrid(self) -> unyt.unyt_array:
        r"""Radial grid in :math:`\rho` which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("rho", (self.nrho,), "1")
        assert self._file is not None
        return self._file.read("rhogrid")

    @property
    def rhogrid(self) -> unyt.unyt_array:
        r"""Time grid in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("rho", (self.ntime,), "s")
        assert self._file is not None
        return self._file.read("timegrid")

    @property
    def ni(self) -> unyt.unyt_array:
        """Density for each ion species."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray(
                "density",
                (self.nion + 1, self.nrho, self.ntime),
                "m**(-3)",
            )
            return data.T[:, 1:]
        assert self._file is not None
        return self._file.read("ni")

    @property
    def Ti(self) -> unyt.unyt_array:
        """Ion temperature."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray(
                "temperature", (self.nrho, self.ntime, 2), "J"
            )
            return data[:, 0].to("eV")
        assert self._file is not None
        return self._file.read("Ti")

    @property
    def ne(self) -> unyt.unyt_array:
        """Electron density."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray(
                "density",
                (self.nion + 1, self.nrho, self.ntime),
                "m**(-3)",
            )
            return data.T[:, 0]
        assert self._file is not None
        return self._file.read("ne")

    @property
    def Te(self) -> unyt.unyt_array:
        """Electron temperature."""
        if self._cdata is not None:
            data = self._cdata.readonly_carray(
                "temperature", (self.nrho, self.ntime, 2), "J"
            )
            return data[:, 1].to("eV")
        assert self._file is not None
        return self._file.read("Te")

    @property
    def charge(self) -> unyt.unyt_array:
        """Ion charge states."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("charge", (self.nion,), "C").to("e")
        assert self._file is not None
        return self._file.read("charge")

    @property
    def rotation(self) -> unyt.unyt_array:
        """Toroidal rotation of the plasma."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("vtor", (self.nrho, self.ntime), "rad/s")
        assert self._file is not None
        return self._file.read("rotation")

    # pylint: disable=too-many-arguments
    def _stage(
        self,
        species: list[Species],
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
        if LIBASCOT.PlasmaDynamic1D_init(
            ctypes.byref(self._cdata),
            rhogrid.size,
            len(species),
            anum,
            znum,
            mass.to("kg").v,
            charge.to("C").v.astype("f8"),
            rhogrid.v,
            Te.to("J").v,
            Ti.to("J").v,
            ne.v,
            ni.v,
            rotation.v,
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "rhogrid",
            "timegrid",
            "ni",
            "Ti",
            "ne",
            "Te",
            "charge",
            "rotation",
        ]:
            self._file.write(field, getattr(self, field))

        self._file.write("anum", self.anum)
        self._file.write("znum", self.znum)

    def export(self) -> dict[str, unyt.unyt_array | list[Species]]:
        fields = [
            "rhogrid",
            "timegrid",
            "ni",
            "Ti",
            "ne",
            "Te",
            "charge",
            "rotation",
            "species",
        ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            species=self.species,
            rhogrid=self.rhogrid,
            timegrid=self.timegrid,
            ni=self.ni,
            Ti=self.Ti,
            ne=self.ne,
            Te=self.Te,
            charge=self.charge,
            rotation=self.rotation,
        )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.PlasmaDynamic1D_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    # pylint: disable=protected-access, too-many-arguments
    def create_plasmadynamic1d(
        self,
        species: List[str] | Tuple[str],
        rhogrid: unyt.unyt_array,
        timegrid: unyt.unyt_array,
        ni: unyt.unyt_array,
        Ti: unyt.unyt_array,
        ne: Optional[unyt.unyt_array] = None,
        Te: Optional[unyt.unyt_array] = None,
        charge: Optional[unyt.unyt_array] = None,
        rotation: Optional[unyt.unyt_array] = None,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
    ) -> PlasmaDynamic1D:
        r"""Create radial plasma profiles that evolve with time.

        This is the dynamic version of :class:`~a5py.data.plasma.Plasma1D`. The
        data is interpolated linearly both in space and time.

        Parameters
        ----------
        species : list[str] or tuple[str] (nspecies,)
            Name(s) of the ion species.
        rhogrid : array_like (nrho,)
            Radial grid in rho in which the data is tabulated.

            This grid doesn't have to be uniform.
        timegrid : array_like (nrho,)
            Time grid in which the data is tabulated.

            This grid doesn't have to be uniform.
        ni : array_like (nrho,nspecies)
            Density for each ion species.
        Ti : array_like (nrho,)
            Ion temperature.
        ne : array_like (nrho,), optional
            Electron density.

            By default, the electron density is determined from the ion charge
            density so that the plasma is quasi-neutral.
        Te : array_like (nrho,), optional
            Electron temperature.

            Same as ion temperature by default.
        charge : array_like (nspecies,), optional
            Ion charge states.

            Ions are fully ionized by default.
        rotation : array_like (nrho,1), optional
            Toroidal rotation of the plasma.

            Zero by default.
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        inputdata : ~a5py.data.plasma.Plasma1DDynamic
            Freshly minted input data object.
        """
        species = [
            s if isinstance(s, Species) else Species.from_string(s) for s in species
        ]
        nion = len(species)
        znum = np.array([s.znum for s in species])

        with utils.validate_variables() as v:
            rhogrid = v.validate("rhogrid", rhogrid, (-1,), "m")
            timegrid = v.validate("timegrid", timegrid, (-1,), "m")

        nrho, ntime = rhogrid.size, timegrid.size
        ni = utils.scalar2array(ni, (nrho, ntime, nion))
        Ti = utils.scalar2array(Ti, (nrho, ntime))
        with utils.validate_variables() as v:
            ni = v.validate("ni", ni, (nrho, ntime, nion), "m**(-3)")
            Ti = v.validate("Ti", Ti, (nrho, ntime), "eV")
            charge = v.validate("charge", charge, (nion,), "e", default=znum)
            rotation = v.validate(
                "rotation",
                rotation,
                (nrho, ntime),
                "rad/s",
                default=np.full(nrho, ntime, 0),
            )

        if charge is None:
            charge_density = np.matmul(ni, znum)
        else:
            charge_density = np.matmul(ni, charge) / unyt.e

        ne = utils.scalar2array(ne, (nrho, ntime))
        Te = utils.scalar2array(Te, (nrho, ntime))
        with utils.validate_variables() as v:
            ne = v.validate("ne", ne, (nrho, ntime), "m**(-3)", default=charge_density)
            Te = v.validate("Te", Te, (nrho, ntime), "eV", default=Ti.v)

        utils.validate_abscissa(rhogrid, "rhogrid", uniform=False)
        utils.validate_abscissa(timegrid, "timegrid", uniform=False)
        leaf = PlasmaDynamic1D(note=note)
        leaf._stage(
            species=species,
            rhogrid=rhogrid,
            timegrid=timegrid,
            ni=ni,
            Ti=Ti,
            ne=ne,
            Te=Te,
            charge=charge,
            rotation=rotation,
        )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf,
            activate=activate,
            save=save,
            category="plasma",
        )
        return leaf
