"""Moment diagnostics.
Collects 1D radial profiles of distribution moments
"""
import ctypes
from typing import Optional

import unyt
import numpy as np

from a5py.libascot import DataStruct
from a5py.data.hist import HistAxis
from a5py.data.options import MomentParams, MOMENTS


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the ``DiagMoment`` struct in diag_data.h."""

    _fields_ = [
        ("axis", HistAxis),
        ("values", ctypes.POINTER(ctypes.c_double)),
    ]


class Moment():
    """Interface to moment diagnostics data.

    This class provides access to the moment data, whether it is stored in
    a file or in memory.
    """

    def __init__(self):
        self._file = None
        self._cdata = None
        self._params = None

    @property
    def rho(self) -> unyt.unyt_array:
        """Rho bin edges."""
        return unyt.unyt_array(
            np.linspace(0.0, 1.0, self._params.rho_bins + 1), "dimensionless"
            )

    @property
    def values(self) -> np.ndarray:
        """Collected moments, shape (rho_bins, number of moments).

        The columns are in the same order as in ``MOMENTS``.
        """
        if self._file is None:
            return self._cdata.values_ref
        return self._file.read("values")

    def save(self, file):
        """Write the data to file and release the C data."""
        file.write("values", self.values)
        del self._cdata
        self._file = file

    def extract(
            self,
            *moments: str,
            volume: Optional[np.ndarray | unyt.unyt_array]=None,
            ) -> unyt.unyt_array | tuple[unyt.unyt_array, ...]:
        """Return radial profiles of the collected moments.

        Parameters
        ----------
        *moments : str
            Names of the moments or "rho" for the rho bin
            centers.
        volume : array_like (rho_bins,), optional
            Volume of each rho shell [m^3].

            If given, the collected sums are divided by it to give the physical
            moments (e.g. particle density). Otherwise the collected sums are
            returned as is.

        Returns
        -------
        *val : unyt_array (rho_bins,) or tuple(unyt_array)
            The moment profiles in the order they were queried.
        """
        names = list(MOMENTS.keys())
        invalid = [m for m in moments if m not in names and m != "rho"]
        if invalid:
            raise ValueError(
                f"Unknown moment(s): {invalid}. Possible values are: "
                f"{['rho'] + names}."
                )

        if volume is not None:
            volume = np.asanyarray(volume)
            if volume.shape != (self._params.rho_bins,):
                raise ValueError(
                    f"volume must have shape ({self._params.rho_bins},), "
                    f"got {volume.shape}."
                    )
            if not isinstance(volume, unyt.unyt_array):
                volume = volume * unyt.m**3

        values = self.values
        out = []
        for moment in moments:
            if moment == "rho":
                edges = self.rho
                out.append((edges[1:] + edges[:-1]) / 2)
                continue
            val = unyt.unyt_array(values[:, names.index(moment)], MOMENTS[moment])
            if volume is not None:
                val = val / volume
            out.append(val)
        return out[0] if len(out) == 1 else tuple(out)

    @classmethod
    def from_params(cls, params: MomentParams) -> Optional["Moment"]:
        """Initialize moment diagnostics from simulation parameters.

        Parameters
        ----------
        params : MomentParams
            Moment diagnostics parameters.

        Returns
        -------
        moment : Moment or None
            The initialized diagnostics or None if the diagnostics is disabled.
        """
        if params.rho_bins == 0:
            return None

        cdata = Struct()
        cdata.axis.min = 0.0
        cdata.axis.max = 1.0
        cdata.axis.n = params.rho_bins

        cdata.values_ref = np.zeros((params.rho_bins, len(MOMENTS)))
        cdata.values = cdata.values_ref.ctypes.data_as(
            ctypes.POINTER(ctypes.c_double)
            )
        cdata.values_ref.flags.writeable = False

        moment = cls()
        moment._cdata = cdata
        moment._params = params
        return moment
