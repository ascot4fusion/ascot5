"""Defines Plasma1DDynamic time-dependent radial plasma input class and
the corresponding factory method.
"""
import ctypes
from typing import Tuple, List, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, interp2D_data, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin



# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in plasma_1Dt.h."""

    _fields_ = [
        ('nrho', ctypes.c_int32),
        ('ntime', ctypes.c_int32),
        ('nspecies', ctypes.c_int32),
        ('anum', ctypes.POINTER(ctypes.c_int32)),
        ('znum', ctypes.POINTER(ctypes.c_int32)),
        ('mass', ctypes.POINTER(ctypes.c_double)),
        ('charge', ctypes.POINTER(ctypes.c_double)),
        ('rho', ctypes.POINTER(ctypes.c_double)),
        ('time', ctypes.POINTER(ctypes.c_double)),
        ('temp', ctypes.POINTER(ctypes.c_double)),
        ('dens', ctypes.POINTER(ctypes.c_double)),
        ('vtor', ctypes.POINTER(ctypes.c_double)),
        ]


@Leaf.register
class PlasmaDynamic1D(InputVariant):
    """Time-dependent radial plasma profile."""

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Plasma1DDynamic",
            struct=Plasma1DDynamic.Struct(),
            )
        self._charge: unyt.unyt_array
        self._species: tuple[str]
        self._rhogrid: unyt.unyt_array
        self._timegrid: unyt.unyt_array
        self._rotation: unyt.unyt_array
        self._electrondensity: unyt.unyt_array
        self._iondensity: unyt.unyt_array
        self._electrontemperature: unyt.unyt_array
        self._iontemperature: unyt.unyt_array

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in which the data is tabulated."""
        if self._staged:
            nrho = self._from_struct_("n_rho", shape=())
            return self._from_struct_("rho", shape=(nrho,), units="1")
        if self._format == Format.HDF5:
            return self._read_hdf5("rhogrid")
        return self._rhogrid.copy()

    @property
    def timegrid(self) -> unyt.unyt_array:
        """Time grid in which the data is tabulated."""
        if self._staged:
            ntime = self._from_struct_("n_time", shape=())
            return self._from_struct_("time", shape=(ntime,), units="s")
        if self._format == Format.HDF5:
            return self._read_hdf5("timegrid")
        return self._timegrid.copy()

    @property
    def rotation(self) -> unyt.unyt_array:
        """Toroidal rotation of the plasma."""
        if self._staged:
            nrho, ntime = self.rhogrid.size, self.timegrid.size
            return self._from_struct_("vtor", shape=(nrho,ntime), units="rad/s")
        if self._format == Format.HDF5:
            for key in ["vtor", "rotation"]:
                try:
                    return self._read_hdf5(key)
                except KeyError:
                    continue
            raise KeyError(
                "Unable to synchronously open object (object 'rotation' "
                "doesn'texist)")
        return self._rotation.copy()

    @property
    def charge(self) -> unyt.unyt_array:
        """Ion charge states."""
        if self._staged:
            nion = len(self.species)
            return self._from_struct_("charge", shape=(nion,), units="C").to("e")
        if self._format == Format.HDF5:
            return self._read_hdf5("charge")
        return self._charge.copy()

    @property
    def species(self):
        """Names of the ion species."""
        if self._staged:
            nion = self._from_struct_("n_species", shape=()) - 1
            anum = self._from_struct_("anum", shape=(nion,))
            znum = self._from_struct_("znum", shape=(nion,))
            return [physlib.properties2species(anum[i], znum[i])
                    for i in range(nion)]
        if self._format == Format.HDF5:
            anum, znum = self._read_hdf5("anum"), self._read_hdf5("znum")
            return [physlib.properties2species(anum[i], znum[i])
                    for i in range(anum.size)]
        return self._species.copy()

    @property
    def electrondensity(self) -> unyt.unyt_array:
        """Radial grid in which the data is tabulated."""
        if self._staged:
            nrho, ntime, nspecies = (
                self.rhogrid.size, self.timegrid.size, len(self.species)
                )
            data = self._from_struct_(
                "dens", shape=(nspecies,ntime,nrho), units="m**(-3)"
                )
            return data.T[:,:,0]
        if self._format == Format.HDF5:
            return self._read_hdf5("electrondensity")
        return self._electrondensity.copy()

    @property
    def electrontemperature(self) -> unyt.unyt_array:
        """Radial grid in which the data is tabulated."""
        if self._staged:
            nrho, ntime = self.rhogrid.size, self.timegrid.size
            data = self._from_struct_(
                "temp", shape=(nrho,ntime,2), units="J"
                ).to("eV")
            return data[:,:,0]
        if self._format == Format.HDF5:
            return self._read_hdf5("electrontemperature")
        return self._electrontemperature.copy()

    @property
    def iondensity(self) -> unyt.unyt_array:
        """Radial grid in which the data is tabulated."""
        if self._staged:
            nrho, ntime, nspecies = (
                self.rhogrid.size, self.timegrid.size, len(self.species)
                )
            data = self._from_struct_(
                "dens", shape=(nspecies,ntime,nrho), units="m**(-3)"
                )
            return data.T[:,:,1:]
        if self._format == Format.HDF5:
            return self._read_hdf5("iondensity")
        return self._iondensity.copy()

    @property
    def iontemperature(self) -> unyt.unyt_array:
        """Radial grid in which the data is tabulated."""
        if self._staged:
            nrho = self.rhogrid.size
            nrho, ntime = self.rhogrid.size, self.timegrid.size
            data = self._from_struct_(
                "temp", shape=(nrho,ntime,2), units="J"
                ).to("eV")
            return data[:,:,1]
        if self._format == Format.HDF5:
            return self._read_hdf5("iontemperature")
        return self._iontemperature.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        data = self.export()
        data["anum"], data["znum"] = [], []
        for species in data["species"]:
            s = physlib.species2properties(species)
            data["anum"].append(s.anum)
            data["znum"].append(s.znum)
        del data["species"]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "rhogrid":self.rhogrid,
            "timegrid":self.timegrid,
            "species":self.species,
            "iondensity":self.iondensity,
            "iontemperature":self.iontemperature,
            "electrondensity":self.electrondensity,
            "electrontemperature":self.electrontemperature,
            "charge":self.charge,
            "rotation":self.rotation,
        }
        return data

    def stage(self):
        init = LIBASCOT.plasma_1Dt_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_int32,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_int32),
            ndpointer(ctypes.c_int32),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            ns = len(self.species)
            anum, znum, mass = (
                np.zeros((ns,), dtype="int32"),
                np.zeros((ns,), dtype="int32"),
                np.zeros((ns,), dtype="f8") * unyt.amu,
                )
            for i, s in enumerate(self.species):
                species = physlib.species2properties(s)
                anum[i], znum[i], mass[i] = (
                    species.anum, species.znum, species.mass
                    )
            if init(
                ctypes.byref(self._struct_),
                self.rhogrid.size,
                self.timegrid.size,
                ns,
                self.rhogrid.v,
                self.timegrid.v,
                anum,
                znum,
                mass.to("kg").v,
                self.charge.to("C").v.astype("float64"),
                self.electrontemperature.to("J").v,
                self.iontemperature.to("J").v,
                self.electrondensity.v,
                self.iondensity.v,
                self.rotation.v,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._rotation
                del self._electrondensity
                del self._iondensity
                del self._electrontemperature
                del self._iontemperature
            self._staged = True

    def unstage(self):
        free = LIBASCOT.plasma_1Dt_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._rotation = self.rotation
                self._electrondensity = self.electrondensity
                self._iondensity = self.iondensity
                self._electrontemperature = self.electrontemperature
                self._iontemperature = self.iontemperature
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments
    def create_plasmadynamic1d(
            self,
            species: List[str] | Tuple[str] | None = None,
            rhogrid: unyt.unyt_array | None = None,
            timegrid: unyt.unyt_array | None = None,
            iondensity: unyt.unyt_array | None = None,
            iontemperature: unyt.unyt_array | None = None,
            electrondensity: Optional[unyt.unyt_array] = None,
            electrontemperature: Optional[unyt.unyt_array] = None,
            charge: Optional[unyt.unyt_array] = None,
            rotation: Optional[unyt.unyt_array] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
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
        iondensity : array_like (nrho,nspecies)
            Density for each ion species.
        iontemperature : array_like (nrho,)
            Ion temperature.
        electrondensity : array_like (nrho,), optional
            Electron density.

            By default, the electron density is determined from the ion charge
            density so that the plasma is quasi-neutral.
        electrontemperature : array_like (nrho,), optional
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
        parameters = _variants.parse_parameters(
            species, rhogrid, timegrid, iondensity, iontemperature,
            electrondensity, electrontemperature, charge, rotation,
        )
        if parameters["species"] is not None:
            nion = parameters["species"].size
            for s in parameters["species"]:
                try:
                    physlib.species2properties(s).znum
                except KeyError as e:
                    raise e from None
        else:
            nion = 2
        default_rhogrid, default_timegrid = (
            np.linspace(0., 1., 3), np.linspace(0., 1., 5)
            )
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        ntime = (default_timegrid.size if parameters["timegrid"] is None
                else parameters["timegrid"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rhogrid", "timegrid", "iondensity", "iontemperature",
                   "species"],
            units=["1", "s", "m**(-3)", "eV", ""],
            shape=[(nrho,), (ntime,), (nrho,ntime,nion), (nrho,ntime), (nion,)],
            dtype=["f8", "f8", "f8", "f8", "s"],
            default=[
                default_rhogrid, default_timegrid, np.ones((nrho,ntime,nion)),
                np.ones((nrho,ntime)), np.array(["H1", "H2"]),
                ],
        )
        znum = []
        for s in parameters["species"]:
            znum.append(physlib.species2properties(s).znum)
        znum = unyt.unyt_array(znum)
        if parameters["charge"] is None:
            charge_density = np.matmul(parameters["iondensity"], znum) / unyt.e
        else:
            charge_density = np.matmul(
                parameters["iondensity"], parameters["charge"]
                ) / unyt.e

        _variants.validate_optional_parameters(
            parameters,
            ["electrondensity", "electrontemperature", "charge", "rotation"],
            ["m**(-3)", "eV", "e", "rad/s"],
            [(nrho,ntime), (nrho,ntime), (nion,), (nrho,ntime)],
            ["f8", "f8", "i4", "f8"],
            [charge_density.v, parameters["iontemperature"].v, znum,
             np.zeros((nrho,ntime)),],
        )
        meta = _variants.new_metadata("Plasma1DDynamic", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
