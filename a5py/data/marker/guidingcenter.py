"""Defines guiding center marker input class :class:`GuidingCenterMarker` and
the corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils
from ...exceptions import AscotIOException


class GuidingCenterMarker(InputVariant):
    """Marker input in guiding-center (5D) phase-space."""

    


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="GuidingCenterMarker",
            struct=GuidingCenterMarker.Struct(),
            )
        self._rhogrid: unyt.unyt_array
        self._dvdrho: unyt.unyt_array

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._format == Format.HDF5:
            nrho, rho0, rho1 = self._read_hdf5("nrho", "rhomin", "rhomax")
            return np.linspace(rho0, rho1, nrho)
        if self._format == Format.CSTRUCT:
            return self._from_struct_("dV", shape=(1,), units="m")
        return self._rhogrid.copy()

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
        """Return a dictionary with sufficient data to duplicate this instance.

        Returns
        -------
        data : dict[str, np.ndarray or unyt.unyt_array]
            Data that can be passed to create_etc to duplicate this instance.
        """
        data = {
            "dvdrho":self.dvdrho,
        }
        return data


# pylint: disable=too-few-public-methods
class CreateGuidingcenterMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create `GuidingCenterMarker` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_gcmarker(
            self,
            species: str | None = None,
            ids: unyt.unyt_array | None = None,
            charge: unyt.unyt_array | None = None,
            r: unyt.unyt_array | None = None,
            phi: unyt.unyt_array | None = None,
            z: unyt.unyt_array | None = None,
            ekin: unyt.unyt_array | None = None,
            pitch: unyt.unyt_array | None = None,
            gyroangle: unyt.unyt_array | None = None,
            weight: unyt.unyt_array | None = None,
            time: Optional[unyt.unyt_array] | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> GuidingCenterMarker:
        r"""Create marker input in guiding-center (5D) phase-space.

        This input can be used in all simulation modes as long as the marker
        is not neutral, but for gyro-orbit simulations the
        :class:`ParticleMarker` input is recommended.

        Parameters
        ----------
        species : str
            Marker species.
        ids : array_like (n,)
            Unique identifier for each marker (must be a positive integer).
        charge : array_like (n,)
            Charge state.
        r : array_like (n,)
            Guiding center R coordinate.
        phi : array_like (n,)
            Guiding center phi coordinate.
        z : array_like (n,)
            Guiding center z coordinate.
        ekin : array_like (n,)
            Kinetic energy.
        pitch : array_like (n,)
            Pitch.
        gyroangle : array_like (n,)
            Gyro angle.

            Usually this doesn't matter as the gyro angle is not relevant for
            guiding center dynamics. This is only used if the simulation
            includes a guiding-center to particle transformation.
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
        inputdata : GuidingCenterMarker
            Freshly minted input data object.
        """
        parameters = variants.parse_parameters(
            ids, charge, r, phi, z, ekin, pitch, gyroangle, weight,
        )
        default_rhogrid = np.linspace(0., 1., 3)
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        variants.validate_required_parameters(
            parameters,
            names=["ids", "charge", "r", "phi", "z", "ekin", "pitch",
                   "gyroangle", "weight",],
            units=["1", "1", "m", "deg", "m", "eV", "1", "rad",],
            shape=(n,),
            dtype="f8",
            default=[np.array([0., 1.]), np.zeros((2,)), 1.],
        )
        meta = variants.new_metadata("GuidingCenterMarker", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
