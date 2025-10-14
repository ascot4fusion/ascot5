"""Defines Neutral3D arbitrary neutral density input class and the corresponding
factory method.
"""
import ctypes
from typing import Tuple, List, Optional

import unyt
import numpy as np

from a5py import utils
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


@Leaf.register
class Neutral3D(InputVariant):
    """Arbitrary neutral density in 3D."""

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Neutral3D",
            struct=Neutral3D.Struct(),
            )
        self._species: Tuple[str]
        self._rgrid: unyt.unyt_array
        self._phigrid: unyt.unyt_array
        self._zgrid: unyt.unyt_array
        self._density: unyt.unyt_array
        self._temperature: unyt.unyt_array

    @property
    def rgrid(self):
        """The uniform grid in R in which data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.n0[0].x_min,
                self._struct_.n0[0].x_max,
                self._struct_.n0[0].n_x
                ) * unyt.m
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nr", "rmin", "rmax")
            return np.linspace(r0, r1, nr)
        return self._rgrid.copy()

    @property
    def phigrid(self):
        """The uniform grid in phi in which data is tabulated."""
        if self._staged:
            return (np.linspace(
                self._struct_.n0[0].y_min,
                self._struct_.n0[0].y_max,
                self._struct_.n0[0].n_y
                ) * unyt.rad).to("deg")
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nphi", "phimin", "phimax")
            return np.linspace(r0, r1, nr)
        return self._phigrid.copy()

    @property
    def zgrid(self):
        """The uniform grid in z in which data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.n0[0].z_min,
                self._struct_.n0[0].z_max,
                self._struct_.n0[0].n_z
                ) * unyt.m
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        return self._zgrid.copy()

    @property
    def temperature(self):
        """Species-wise temperature."""
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
        """Species-wise density."""
        if self._staged:
            nspecies = self._from_struct_("n_species", shape=())
            data = self._from_struct_("n0", idx=0)
            for i in range(1, nspecies):
                data = np.stack((data, self._from_struct_("n0", idx=i)), axis=3)
            if nspecies == 1:
                data = np.expand_dims(data, axis=3)
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
        for grid in ["rgrid", "phigrid", "zgrid"]:
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
            "rgrid":self.rgrid,
            "zgrid":self.zgrid,
            "phigrid":self.phigrid,
            "species":self.species,
            "density":self.density,
            "temperature":self.temperature,
        }
        return data

    def stage(self):
        init = LIBASCOT.N0_3D_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
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
                self.rgrid.size,
                self.rgrid[0].v,
                self.rgrid[-1].v,
                self.phigrid.size,
                self.phigrid[0].to("rad").v,
                self.phigrid[-1].to("rad").v,
                self.zgrid.size,
                self.zgrid[0].v,
                self.zgrid[-1].v,
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
        free = LIBASCOT.N0_3D_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._density = self.density
                self._temperature = self.temperature
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create Neutral3D input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_neutral3d(
            self,
            rgrid: utils.ArrayLike | None = None,
            phigrid: utils.ArrayLike | None = None,
            zgrid: utils.ArrayLike | None = None,
            species: List[str] | Tuple[str] | None = None,
            density: utils.ArrayLike | None = None,
            temperature: utils.ArrayLike | None = None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> Neutral3D:
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
        parameters = _variants.parse_parameters(
            rgrid, phigrid, zgrid, species, density, temperature,
        )
        if parameters["species"] is not None:
            ns = parameters["species"].size
            for s in parameters["species"]:
                try:
                    physlib.species2properties(s)
                except KeyError as e:
                    raise e from None
        else:
            ns = 1
        default_rgrid, default_phigrid, default_zgrid = (
            np.linspace(1, 2, 45),
            np.linspace(0, 360, 31)[:-1],
            np.linspace(-1, 1, 90),
        )
        nr = (default_rgrid.size if parameters["rgrid"] is None
              else parameters["rgrid"].size)
        nphi = (default_phigrid.size if parameters["phigrid"] is None
              else parameters["phigrid"].size)
        nz = (default_zgrid.size if parameters["zgrid"] is None
              else parameters["zgrid"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rgrid", "phigrid", "zgrid", "density", "temperature",
                   "species"],
            units=["m", "deg", "m", "m**(-3)", "eV", ""],
            shape=[(nr,), (nphi,), (nz,), (nr,nphi,nz,ns), (nr,nphi,nz,ns),
                   (ns,)],
            dtype=["f8", "f8", "f8", "f8", "f8", "s"],
            default=[default_rgrid, default_phigrid, default_zgrid,
                     np.ones((nr,nphi,nz,ns)), np.ones((nr,nphi,nz,ns)),
                     np.array(["H1"])],
        )
        meta = _variants.new_metadata("Neutral3D", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
