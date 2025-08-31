import numpy as np

from a5py.exceptions import AscotNoDataException

import a5py.routines.plotting as a5plt
#from a5py.data import Dist
#from a5py.data.dist import DistMoment
import a5py.physlib as physlib

class DistMixin():

    def _getdist(self, distout, mass, exi=False, ekin_edges=None,
                 pitch_edges=None, plotexi=False):
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
        if exi:
            if "ppar" in distout.abscissae and "pper":
                if plotexi: dist0 = distout._copy()
                if ekin_edges is None:
                    ekin_edges  = int(distout.abscissa("ppar").size / 2)
                if pitch_edges is None:
                    pitch_edges = distout.abscissa("pperp").size
                distout = Dist.ppappe2ekinpitch(
                    distout, mass, ekin_edges=ekin_edges,
                    pitch_edges=pitch_edges)
            else:
                raise ValueError(
                    "Energy-pitch transformation only valid for 5d")

        if exi and plotexi:
            integrate = {}
            for k in dist0.abscissae:
                if k not in ["ppar", "pperp"]:
                    integrate[k] = np.s_[:]
            dist0.integrate(**integrate)
            integrate = {}
            for k in distout.abscissae:
                if k not in ["ekin", "pitch"]:
                    integrate[k] = np.s_[:]
            dist1 = distout.integrate(copy=True, **integrate)

            g = physlib.gamma_energy(mass, dist1.abscissa_edges("ekin"))
            pnorm_edges = physlib.pnorm_gamma(mass, g).to("kg*m/s")
            pitch_edges = dist1.abscissa_edges("pitch")

            fig = a5plt.plt.figure()
            ax1 = fig.add_subplot(3,1,1)
            ax2 = fig.add_subplot(3,1,2, projection='polar')
            ax3 = fig.add_subplot(3,1,3)

            dist0.plot(axes=ax1)
            a5plt.momentumpolargrid(pnorm_edges, pitch_edges, axes=ax1)
            a5plt.momentumpolarplot(pnorm_edges, pitch_edges,
                                    dist1.distribution(), axes=ax2)
            dist1.plot(axes=ax3)
            a5plt.plt.show()

        return distout
