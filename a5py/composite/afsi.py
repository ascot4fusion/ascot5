class AfsiRun():

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
        """Return 5D distribution function of one of the fusion products.

        Parameters
        ----------
        dist : {"prod1", "prod2"}
            Which product to return.
        exi : bool, optional
            Convert the momentum space to energy-pitch.

            The distribution is normalized to conserve the particle number.
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
        data : :class:`DistData`
            The distribution data object.
        """
        if dist == "prod1":
            self._require("_prod1dist5d", "_reaction")
            distout = self._prod1dist5d.get()
            mass = self._reaction.get()[2]
        elif dist == "prod2":
            self._require("_prod2dist5d", "_reaction")
            distout = self._prod2dist5d.get()
            mass = self._reaction.get()[3]
        else:
            raise ValueError("dist must be either 'prod1' or 'prod2'")

        return self._getdist(distout, mass, exi=exi, ekin_edges=ekin_edges,
                             pitch_edges=pitch_edges, plotexi=plotexi)
