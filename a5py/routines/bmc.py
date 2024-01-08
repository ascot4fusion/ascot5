"""Class for BMC post-processing.
"""
import numpy as np
import unyt
import a5py.routines.plotting as a5plt
from a5py.routines.distmixin import DistMixin
from a5py.ascot5io import Marker, State

class BMCMixin(DistMixin):
    """Mixin class with post-processing results related to BMC.
    """

    def _require(self, *args):
        """Check if required data is present and raise exception if not.

        This is a helper function to quickly check that the data is available.

        Parameters
        ----------
        *args : `str`
            Name(s) of the required data.

        Raises
        ------
        AscotNoDataException
            Raised if the required data is not present.
        """
        for arg in args:
            if not hasattr(self, arg):
                raise AscotNoDataException(
                    "Data for \"" +  arg + "\" is required but not present.")

    def getdist(self, dist, exi=False, ekin_edges=None, pitch_edges=None,
                plotexi=False):
        """Return distribution function.

        Parameters
        ----------
        dist : {"5d", "6d", "rho5d", "rho6d", "com"}
            Name of the distribution.

            If ``dist`` is ``None``, returns a list of available distributions.
        exi : bool, optional
            Convert the momentum space to energy-pitch.

            The distribution is normalized to conserve the particle number.

            Not applicable for ``dist`` = {"6d", "rho6d", "com"}.
        ekin_edges : int or array_like, optional
            Number of bins or bin edges in the energy abscissa if ``exi=True``.
        pitch_edges : int or array_like, optional
            Number of bins or bin edges in the pitch abscissa if ``exi=True``.
        plotexi : bool, optional
            Visualize the transformation from ppar-perp to energy-pitch if
            if ``exi=True``.

            Use this option to adjust energy and pitch abscissae to your liking.

        Returns
        -------
        data : :class:`DistData` or list [str]
            The distribution data object or a list of available distributions if
            ``dist`` is ``None``.
        """
        if dist is None:
            dists = ["5d", "6d", "rho5d", "rho6d", "com"]
            for d in dists:
                try:
                    self._require("_dist" + d)
                except AscotNoDataException:
                    dists.remove(d)
            return dists

        if dist == "5d":
            self._require("_dist5d")
            distout = self._dist5d.get()
        elif dist == "6d":
            self._require("_dist6d")
            distout = self._dist6d.get()
        elif dist == "rho5d":
            self._require("_distrho5d")
            distout = self._distrho5d.get()
        elif dist == "rho6d":
            self._require("_distrho6d")
            distout = self._distrho6d.get()
        elif dist == "com":
            self._require("_distcom")
            distout = self._distcom.get()
        else:
            raise ValueError("Unknown distribution")

        mass = self.options.read()["BMC_MASS"] * unyt.amu
        return self._getdist(distout, mass, exi=exi, ekin_edges=ekin_edges,
                             pitch_edges=pitch_edges, plotexi=plotexi)