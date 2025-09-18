"""Defines particle (gyro-orbit) input class :class:`ParticleMarker` and the
corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from ..cstructs import free
from ..access import _variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils, physlib
from ...exceptions import AscotIOException

from .state import MarkerState

class ParticleMarker(InputVariant):
    """Marker input in particle (6D) phase-space."""


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="ParticleMarker",
            struct=None,
            )

    @property
    def species(self) -> str:
        """Marker species."""
        return self._species

    @property
    def ids(self) -> np.ndarray:
        """Unique identifier for each marker."""
        if self._format == Format.HDF5:
            return self._read_hdf5("ids")

        out = np.zeros((len(self._struct_),))
        for i in range(out.size):
            out[i] = self._struct_[i].id
        return out

    @property
    def r(self) -> unyt.unyt_array:
        r"""Field-line :math:`R` coordinate."""
        if self._format == Format.HDF5:
            return self._read_hdf5("r")

        out = np.zeros((len(self._struct_),)) * unyt.m
        for i in range(out.size):
            out[i] = self._struct_[i].r
        return out

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Field-line :math:`\phi` coordinate."""
        if self._format == Format.HDF5:
            return self._read_hdf5("phi")

        out = np.zeros((len(self._struct_),)) * unyt.rad
        for i in range(out.size):
            out[i] = self._struct_[i].phi
        return out.to("deg")

    @property
    def z(self) -> unyt.unyt_array:
        r"""Field-line :math:`z` coordinate."""
        if self._format == Format.HDF5:
            return self._read_hdf5("z")

        out = np.zeros((len(self._struct_),)) * unyt.m
        for i in range(out.size):
            out[i] = self._struct_[i].z
        return out

    @property
    def direction(self) -> unyt.unyt_array:
        r"""Field-line direction (negative means opposite to the magnetic field
        vector).
        """
        if self._format == Format.HDF5:
            return self._read_hdf5("direction")

        out = np.zeros((len(self._struct_),)) * unyt.dimensionless
        for i in range(out.size):
            out[i] = self._struct_[i].ppar
        return out

    @property
    def time(self) -> unyt.unyt_array:
        """Field-line time."""
        if self._format == Format.HDF5:
            return self._read_hdf5("time")

        out = np.zeros((len(self._struct_),)) * unyt.s
        for i in range(out.size):
            out[i] = self._struct_[i].time
        return out

    @property
    def n(self) -> int:
        """Number of markers."""
        if self._format == Format.HDF5:
            return self._read_hdf5("ids").size
        return len(self._struct_)

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "species":self.species,
            "r":self.r,
            "phi":self.phi,
            "z":self.z,
            "vr":self.vr,
            "vphi":self.vphi,
            "vz":self.vz,
            "weight":self.weight,
            "time":self.time,
        }
        return data

    def stage(self):
        pass

    def unstage(self):
        pass


# pylint: disable=too-few-public-methods
class CreateParticleMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create :class:`ParticleMarker` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_particlemarker(
            self,
            species: str | None = None,
            ids: unyt.unyt_array | None = None,
            charge: unyt.unyt_array | None = None,
            r: unyt.unyt_array | None = None,
            phi: unyt.unyt_array | None = None,
            z: unyt.unyt_array | None = None,
            vr: unyt.unyt_array | None = None,
            vphi: unyt.unyt_array | None = None,
            vz: unyt.unyt_array | None = None,
            weight: unyt.unyt_array | None = None,
            time: Optional[unyt.unyt_array] | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> ParticleMarker:
        r"""Create marker input in particle (6D) phase-space.

        This is the recommended marker input, especially for simulations that
        trace the whole particle orbit to avoid information loss due to
        guiding-center to particle transformations.

        This input can be used in all simulation modes.

        Parameters
        ----------
        species : str
            Marker species.
        ids : array_like (n,)
            Unique identifier for each marker (must be a positive integer).
        charge : array_like (n,)
            Charge state.
        r : array_like (n,)
            Particle R coordinate.
        phi : array_like (n,)
            Particle phi coordinate.
        z : array_like (n,)
            Particle z coordinate.
        vr : array_like (n,)
            Velocity R-component.
        vphi : array_like (n,)
            Velocity phi-component.
        vz : array_like (n,)
            Velocity z-component.
        weight : array_like (n,)
            How many physical particles are represented by this marker.
        time : array_like (n,), optional
            Time instant (e.g. with respect to the plasma pulse) when the marker
            is created.

            Relevant mainly for time-dependent simulations.
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
        inputdata : ParticleMarker
            Freshly minted input data object.
        """
        parameters = _variants.parse_parameters(
            species, ids, charge, r, phi, z, vr, vphi, vz, weight, time,
        )
        n = 1 if parameters["ids"] is None else parameters["ids"].size
        _variants.validate_required_parameters(
            parameters,
            names=["ids", "charge", "r", "phi", "z", "vr", "vphi", "vz",
                   "weight", "time", "species"],
            units=["1", "e", "m", "deg", "m", "m/s", "m/s", "m/s",
                   "markers/s", "s", ""],
            shape=(n,),
            dtype=["i8", "i4", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8",
                   "s"],
            default=[
                1, 1, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, "H1",
            ],
        )
        try:
            physlib.species2properties(parameters["species"])
        except KeyError as e:
            raise e from None
        meta = _variants.new_metadata("ParticleMarker", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            if parameter != "species":
                getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj

        parameters = _variants.parse_parameters(
            ids, r, phi, z, direction, time,
        )
        n = 1 if parameters["ids"] is None else parameters["ids"].size
        _variants.validate_required_parameters(
            parameters,
            names=["ids", "r", "phi", "z",],
            units=["1", "m", "deg", "m",],
            shape=(n,),
            dtype=["i8", "f8", "f8", "f8",],
            default=[1, 1.0, 0.0, 0.0,],
        )
        _variants.validate_optional_parameters(
            parameters,
            names=["direction", "time"],
            units=["1", "s"],
            shape=(n,),
            dtype=["f8", "f8",],
            default=[np.ones((n,)), np.ones((n,)),],
        )
        meta = _variants.new_metadata("FieldlineMarker", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        obj._struct_ = (MarkerState.Structure * n)()
        parameters.update({
            "id":parameters["ids"], "ppar":parameters["direction"]
            })
        del parameters["ids"]
        del parameters["direction"]
        for key in parameters.keys():
            for i in range(n):
                setattr(obj._struct_[i], key, parameters[key][i])

        if store_hdf5:
            obj._export_hdf5()
        return obj