"""Defines particle (gyro-orbit) input class :class:`ParticleMarker` and the
corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from ..cstructs import free
from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils, physlib
from ...exceptions import AscotIOException

from .cstructs import (INPUT_PARTICLE_TYPE, input_particle,
                       allocate_input_particles)


class ParticleMarker(InputVariant):
    """Marker input in particle (6D) phase-space."""

    Struct = input_particle


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="ParticleMarker",
            struct=ctypes.POINTER(input_particle)(),
            )
        self._species: str
        self._ids: unyt.unyt_array
        self._charge: unyt.unyt_array
        self._r: unyt.unyt_array
        self._phi: unyt.unyt_array
        self._z: unyt.unyt_array
        self._vr: unyt.unyt_array
        self._vphi: unyt.unyt_array
        self._vz: unyt.unyt_array
        self._weight: unyt.unyt_array
        self._time: unyt.unyt_array

    @property
    def species(self) -> str:
        """Marker species."""
        return self._species

    @property
    def ids(self) -> unyt.unyt_array:
        """Unique identifier for each marker."""
        return self._ids.copy()

    @property
    def charge(self) -> unyt.unyt_array:
        """Charge state."""
        return self._charge.copy()

    @property
    def r(self) -> unyt.unyt_array:
        r"""Particle :math:`R` coordinate."""
        if self._format == Format.HDF5:
            nrho, rho0, rho1 = self._read_hdf5("nrho", "rhomin", "rhomax")
            return np.linspace(rho0, rho1, nrho)
        return self._r.copy()

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Particle :math:`\phi` coordinate."""
        return self._phi.copy()

    @property
    def z(self) -> unyt.unyt_array:
        r"""Particle :math:`z` coordinate."""
        return self._z.copy()

    @property
    def vr(self) -> unyt.unyt_array:
        r"""Particle velocity :math:`R` component."""
        return self._vr.copy()

    @property
    def vphi(self) -> unyt.unyt_array:
        r"""Particle velocity :math:`\phi` component."""
        return self._vphi.copy()

    @property
    def vz(self) -> unyt.unyt_array:
        r"""Particle velocity :math:`z` component."""
        return self._vz.copy()

    @property
    def weight(self) -> unyt.unyt_array:
        r"""Particle weight."""
        return self._weight.copy()

    @property
    def time(self) -> unyt.unyt_array:
        r"""Particle time."""
        return self._time.copy()

    @property
    def n(self) -> int:
        """Number of markers."""
        return self.ids.size

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
        n = self.n
        self._struct_ = allocate_input_particles(n)
        particle_type = INPUT_PARTICLE_TYPE["p"]
        properties = physlib.species2properties(self.species)
        anum, znum, mass = properties.anum, properties.znum, properties.mass
        for i in range(n):
            self._struct_[i].type = particle_type
            self._struct_[i].p.anum = anum
            self._struct_[i].p.znum = znum
            self._struct_[i].p.mass = mass
        for fieldname in ["ids", "r", "phi", "z", "vr", "vphi", "vz", "weight",
                          "time"]:
            structname = fieldname if fieldname != "ids" else "id"
            field = getattr(self, fieldname)
            for i in range(n):
                setattr(self._struct_[i].p, structname, field[i])

            if self._format == Format.MEMORY:
                delattr(self, "_"+fieldname)
            self._staged = True

    def unstage(self):
        if self._staged:
            if self._format is Format.MEMORY:
                self._density = self.density
                self._temperature = self.temperature
            free(ctypes.byref(self._struct_))
            self._staged = False


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
        parameters = variants.parse_parameters(
            species, ids, charge, r, phi, z, vr, vphi, vz, weight, time,
        )
        n = 1 if parameters["ids"] is None else parameters["ids"].size
        variants.validate_required_parameters(
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
        meta = variants.new_metadata("ParticleMarker", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            if parameter is not "species":
                getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
