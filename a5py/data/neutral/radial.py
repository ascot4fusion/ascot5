"""Defines Neutral1D radial neutral density input class and the corresponding
factory method.
"""
import ctypes
from typing import Tuple, List, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import _variants, InputVariant, Format, TreeCreateClassMixin
from ..cstructs import linint1D_data
from ... import utils, physlib
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class Neutral1D(InputVariant):
    """Radial neutral profile."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in N0_1D.h."""
        _pack_ = 1
        _fields_ = [
            ('n_species', ctypes.c_int32),
            ('PADDING_0', ctypes.c_ubyte * 4),
            ('anum', ctypes.POINTER(ctypes.c_int32)),
            ('znum', ctypes.POINTER(ctypes.c_int32)),
            ('n0', ctypes.POINTER(linint1D_data)),
            ('t0', ctypes.POINTER(linint1D_data)),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Neutral1D",
            struct=Neutral1D.Struct(),
            )
        self._species: Tuple[str]
        self._n0: unyt.unyt_array
        self._t0: unyt.unyt_array
        self._rhogrid: unyt.unyt_array

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.n0[0].x_min,
                self._struct_.n0[0].x_max,
                self._struct_.n0[0].n_x
                ) * unyt.dimensionless
        if self._format == Format.HDF5:
            nrho, rho0, rho1 = self._read_hdf5("nrho", "rhomin", "rhomax")
            return np.linspace(rho0, rho1, nrho)
        return self._rhogrid.copy()

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
            return data * unyt.m**(-3)
        if self._format == Format.HDF5:
            return self._read_hdf5("density")
        return self._density.copy()

    @property
    def species(self):
        """Names of the neutral species."""
        if self._staged:
            nspecies = self._from_struct_("n_species", shape=())
            anum = self._from_struct_("anum", shape=(nspecies,))
            znum = self._from_struct_("znum", shape=(nspecies,))
            return [physlib.properties2species(anum[i], znum[i])
                    for i in range(nspecies)]
        if self._format == Format.HDF5:
            anum, znum = self._read_hdf5("anum"), self._read_hdf5("znum")
            return [physlib.properties2species(anum[i], znum[i])
                    for i in range(anum.size)]
        return self._species.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        data["anum"], data["znum"] = [], []
        for species in data["species"]:
            s = physlib.species2properties(species)
            data["anum"].append(s.anum)
            data["znum"].append(s.znum)
        del data["species"]
        for grid in ["rhogrid"]:
            name = grid.replace("grid", "")
            data["n" + name] = data[grid].size
            data[name + "min"] = data[grid][0]
            data[name + "max"] = data[grid][-1]
            del data[grid]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "rhogrid":self.rhogrid,
            "species":self.species,
            "density":self.density,
            "temperature":self.temperature,
        }
        return data

    def stage(self):
        init = LIBASCOT.N0_1D_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ndpointer(ctypes.c_int32),
            ndpointer(ctypes.c_int32),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            ns = len(self.species)
            anum, znum = np.zeros((ns,), dtype="int32"), np.zeros((ns,), dtype="int32")
            for i, s in enumerate(self.species):
                species = physlib.species2properties(s)
                anum[i], znum[i] = species.anum, species.znum
            if init(
                ctypes.byref(self._struct_),
                self.rhogrid.size,
                self.rhogrid[0].v,
                self.rhogrid[-1].v,
                ns,
                anum,
                znum,
                self.density.v,
                self.temperature.to("J").v,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._density
                del self._temperature
            self._staged = True

    def unstage(self):
        free = LIBASCOT.N0_1D_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._density = self.density
                self._temperature = self.temperature
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateNeutral1DMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create Neutral1D input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_neutral1d(
            self,
            rhogrid: utils.ArrayLike | None = None,
            species: List[str] | Tuple[str] | None = None,
            density: utils.ArrayLike | None = None,
            temperature: utils.ArrayLike | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> Neutral1D:
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
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        inputdata : ~a5py.data.neutral.Neutral1D
            Freshly minted input data object.
        """
        parameters = _variants.parse_parameters(
            rhogrid, species, density, temperature,
        )
        if parameters["species"] is not None:
            ns = parameters["species"].size
            for s in parameters["species"]:
                try:
                    physlib.species2properties(s)
                except KeyError as e:
                    raise e from None
        else:
            ns = 2

        default_rhogrid = np.linspace(0., 1., 3)
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rhogrid", "density", "temperature", "species"],
            units=["1", "m**(-3)", "eV", ""],
            shape=[(nrho,), (nrho, ns), (nrho, ns), (ns,)],
            dtype=["f8", "f8", "f8", "s"],
            default=[default_rhogrid, np.ones((nrho,ns)), np.ones((nrho,ns)),
                     np.array(["H1", "H2"])],
        )
        for abscissa in ["rhogrid"]:
            utils.check_abscissa(parameters[abscissa], abscissa)
        meta = _variants.new_metadata("Neutral1D", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
