import numpy as np

from a5py.exceptions import AscotNoDataException

import a5py.plotting.plotting as a5plt
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

    mass = np.mean(self.getstate("mass"))
    return self._getdist(distout, mass, exi=exi, ekin_edges=ekin_edges,
                        pitch_edges=pitch_edges, plotexi=plotexi)

def getdist_moments(self, dist, *moments, volmethod="prism"):
    """Calculate moments of distribution.

    Parameters
    ----------
    dist : str or :class:`DistData`
        Distribution from which moments are calculated.
    *moments : str
        Moments to be calculated.
    volmethod : {"mc", "prism"}, optional
        Method used to calculate the volume.

        It is a good idea to verify that same results are obtained with both
        methods.

    Returns
    -------
    out : :class:`DistMoment`
        Distribution object containing moments as ordinates.
    """
    # Initialize the moment object and evaluate the volume.
    if all([x in dist.abscissae for x in ["rho", "theta", "phi"]]):
        rhodist = True
        nrho   = dist.abscissa_edges("rho").size
        ntheta = dist.abscissa_edges("theta").size
        nphi   = dist.abscissa_edges("phi").size
        volume, area, r, phi, z = self._root._ascot.input_rhovolume(
            method=volmethod, tol=1e-1, nrho=nrho, ntheta=ntheta, nphi=nphi,
            return_area=True, return_coords=True)
        volume[volume == 0] = 1e-8 # To avoid division by zero
        out = DistMoment(
            dist.abscissa_edges("rho"), dist.abscissa_edges("theta"),
            dist.abscissa_edges("phi"), r, phi, z, area, volume, rhodist)
    elif all([x in dist.abscissae for x in ["r", "phi", "z"]]):
        rhodist = False
        r = dist.abscissa_edges("r")
        p = dist.abscissa_edges("phi").to("rad").v
        z = dist.abscissa_edges("z")
        v1, v2 = np.meshgrid( r[1:] - r[:-1], z[1:] - z[:-1])
        area = v1 * v2
        v1, v2, v3 = np.meshgrid(
            p[1:] - p[:-1], r[1:]**2 - r[:-1]**2, z[1:] - z[:-1])
        volume = 0.5 * v1 * v2 * v3
        r, phi, z = np.meshgrid(
            dist.abscissa("r"), dist.abscissa("phi"), dist.abscissa("z"))
        out = DistMoment(
            dist.abscissa_edges("r"), dist.abscissa_edges("phi"),
            dist.abscissa_edges("z"), r, phi, z, area, volume, rhodist)
    else:
        raise ValueError(
            "Distribution has neither (rho, theta, phi) nor (R, phi, z)")

    ordinates = {}
    mass = np.mean(self.getstate("mass"))
    if "density" in moments:
        Dist.density(dist, out)
    if "chargedensity" in moments:
        Dist.chargedensity(dist, out)
    if "energydensity" in moments:
        Dist.energydensity(mass, dist, out)
    if "toroidalcurrent" in moments:
        Dist.toroidalcurrent(self._root._ascot, mass, dist, out)
    if "parallelcurrent" in moments:
        Dist.parallelcurrent(self._root._ascot, mass, dist, out)
    if "currentdrive" in moments:
        Dist.parallelcurrent(self._root._ascot, mass, dist, out, drive=True)
    if "pressure" in moments:
        Dist.pressure(mass, dist, out)
    if "powerdep" in moments:
        Dist.powerdep(self._root._ascot, mass, dist, out)
    if "electronpowerdep" in moments:
        Dist.electronpowerdep(self._root._ascot, mass, dist, out)
    if "ionpowerdep" in moments:
        Dist.ionpowerdep(self._root._ascot, mass, dist, out)
    #if "jxbtorque" in moments:
    #    Dist.jxbtorque(self._root._ascot, mass, dist, out)
    #if "colltorque" in moments:
    #    Dist.collTorque(self._root._ascot, mass, dist, out)
    #if "canmomtorque" in moments:
    #    Dist.canMomentTorque(dist, out)
    return out

def getdist_list(self, show=True):
    """List all available distributions and moments.

    Parameters
    ----------
    show : bool, optional
        Print the output.

    Returns
    -------
    dists : list [str]
        List of available distributions.
    moms : list [(str, str)]
        List of distribution moments that can be evaluated and their
        meaning.
    """
    dists = []
    if hasattr(self, "_dist5d"):    dists.append("5d")
    if hasattr(self, "_distrho5d"): dists.append("rho5d")
    if hasattr(self, "_dist6d"):    dists.append("6d")
    if hasattr(self, "_distrho6d"): dists.append("rho6d")
    if hasattr(self, "_distcom"):   dists.append("com")

    moms = []
    if "5d" not in dists and "rho5d" not in dists:
        if show:
            print("Available distributions:")
            for d in dists:
                print(d)
            print("\nMoments available only for 5d/rho5d dists.")
        return dists, moms

    moms = [
        ("density", "Particle number density"),
        ("chargedensity", "Charge density"),
        ("energydensity", "Energy density"),
        ("toroidalcurrent", "Toroidal current"),
        ("parallelcurrent", "Parallel current"),
        ("pressure", "Pressure"),
        ("powerdep", "Total deposited power"),
        ("ionpowerdep", "Power deposited to ions"),
        ("electronpowerdep", "Power deposited to electrons"),
        #("jxbtorque", "j_rad x B_pol torque"),
        #("colltorque", "Torque from collisions"),
        #("canmomtorque", "Torque from change in can. tor. ang. momentum"),
    ]
    if show:
        print("Available distributions:")
        for d in dists:
            print(d)
        print("Available Moments:")
        for d in moms:
            print(d[0] + " : " + d[1])
    return dists, moms
