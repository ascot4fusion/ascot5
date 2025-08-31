"""Defines magnetic field line marker input class :class:`FieldlineMarker` and
the corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils
from ...exceptions import AscotIOException

from .cstructs import particle_state

class FieldlineMarker(InputVariant):
    """Marker input in magnetic field line (3D) phase-space."""

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="FieldlineMarker",
            struct=None,
            )

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
            "r":self.r,
            "phi":self.phi,
            "z":self.z,
            "direction":self.direction,
            "time":self.time,
        }
        return data

    def stage(self):
        pass

    def unstage(self):
        pass


# pylint: disable=too-few-public-methods
class CreateFieldlineMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create `FieldlineMarker` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_fieldlinemarker(
            self,
            ids: unyt.unyt_array | None = None,
            r: unyt.unyt_array | None = None,
            phi: unyt.unyt_array | None = None,
            z: unyt.unyt_array | None = None,
            direction: unyt.unyt_array | None = None,
            time: Optional[unyt.unyt_array] | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> FieldlineMarker:
        r"""Create markers representing magnetic field lines.

        This input can only be used in magnetic field-line tracing.

        Parameters
        ----------
        ids : array_like (n,)
            Unique identifier for each marker (must be a positive integer).
        r : array_like (n,)
            Field line radial coordinate.
        phi : array_like (n,)
            Field line toroidal coordinate.
        z : array_like (n,)
            Field line axial coordinate.
        direction : array_like (n,)
            The direction of the field line: +1 for co-parallel, -1 for
            counter-parallel.
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
        inputdata : FieldlineMarker
            Freshly minted input data object.
        """
        parameters = variants.parse_parameters(
            ids, r, phi, z, direction, time,
        )
        n = 1 if parameters["ids"] is None else parameters["ids"].size
        variants.validate_required_parameters(
            parameters,
            names=["ids", "r", "phi", "z",],
            units=["1", "m", "deg", "m",],
            shape=(n,),
            dtype=["i8", "f8", "f8", "f8",],
            default=[1, 1.0, 0.0, 0.0,],
        )
        variants.validate_optional_parameters(
            parameters,
            names=["direction", "time"],
            units=["1", "s"],
            shape=(n,),
            dtype=["f8", "f8",],
            default=[np.ones((n,)), np.ones((n,)),],
        )
        meta = variants.new_metadata("FieldlineMarker", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        obj._struct_ = (particle_state * n)()
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
