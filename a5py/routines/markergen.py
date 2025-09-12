"""Generate marker populations.
"""
import numpy as np
import unyt

#from a5py.data.marker import Marker
from a5py.data.dist import DistData
from a5py import physlib

class MarkerGenerator():
    """Tool for processing distributions into markers.

    Attributes
    ----------
    _ascot : :class:`Ascot`
        The Ascot object this instance belongs to.
    """

    def __init__(self, ascot):
        self._ascot = ascot

    def generate(self, nmrk, mass, charge, anum, znum, particledist,
                 markerdist=None, mode='gc', minweight=0, return_dists=False):
        """Generate weighted markers from marker and particle distributions.

        This function takes two 5D distributions that must have identical grids.
        The marker distribution is a (normalized) probability distribution from
        which markers are sampled. The particle distribution presents physical
        particles, e.g. it could be an alpha particle birth rate. Once markers
        are sampled from the marker distribution, then each marker in a same
        cell is given same weight, so that those markers represent all
        the physical particles in that particular cell.

        As markers are sampled from the marker distribution, their initial
        coordinates are randomly chosen from the phase space region contained
        within the cell.

        Parameters
        ----------
        nmrk : int
            Number of markers to be generated.
        mass : float
            Particle species` mass.
        charge : float
            Particle species` charge.
        anum : int
            Particle species` atomic mass number.
        znum : int
            Particle species` charge number.
        particledist : :class:`Dist`
            5D distribution from which marker weight is sampled.
        markerdist : :class:`Dist`, optional
            5D distribution from which marker initial coordinates are sampled.

            If not given, the coordinates are sampled from ``particledist``.
        mode : {'prt', 'gc'}, optional
            Decides whether the distributions (and returned markeres) are in
            particle or guiding center phase-space.
        minweight : float, optional
            Minimum weight a marker must have or else it is rejected.
        return_dists : bool, optional
            Return ``markerdist`` and ``particledist`` calculated from
            the generated markers.

            These distributions hould converge to inputs as the marker number is
            increased. The exception is ``particledist`` which will always get
            zero values in regions where ``markerdist`` is zero.

        Returns
        -------
        mrk : dict
            Marker input data on a dictionary.
        mrkdist : :class:`Dist`
            Distribution of initialized markers which ideally should be close to
            the one given as an input unless markers were rejected
            by ``minweight``.
        prtdist : :class:`Dist`
            Weighted distribution of initialized markers which ideally should be
            close to the one given as an input, unless some phase-space regions
            in ``markerdist`` were left intentionally empty even though they
            would contain physical particles.
        """
        if markerdist is None:
            markerdist = particledist
        # Number of markers successfully generated
        ngen      = 0
        # Cell indices of generated markers
        icell     = np.zeros((nmrk,), dtype="i8")

        # Generate a number random for each marker, and when that marker is put
        # in the first cell where rand > threshold.
        threshold = np.append(0, np.cumsum(markerdist.histogram().ravel()))
        threshold /= threshold[-1]
        while ngen < nmrk:
            if ngen == 0: rejected = np.s_[:]
            icell[rejected] = \
                np.digitize( np.random.rand(nmrk-ngen,), bins=threshold ) - 1

            # Each marker is given a weight that corresponds to number of
            # physical particles in that cell, divided by the number of markers
            # in that cell
            _, idx, counts = \
                np.unique(icell, return_inverse=True, return_counts=True)
            weight = particledist.histogram().ravel()[icell] / counts[idx]

            # Reject based on the minweight
            rejected = weight < minweight
            ngen = np.sum(~rejected)

        # Shuffle markers just in case the order they were created is biased
        idx = np.arange(nmrk)
        np.random.shuffle(idx)
        icell  = icell[idx]
        weight = weight[idx].ravel()

        # Init marker species
        mrk = Marker.generate(mode, n=nmrk)
        mrk["anum"][:]   = anum
        mrk["znum"][:]   = znum
        mrk["mass"][:]   = mass
        mrk["charge"][:] = charge
        mrk["weight"][:] = weight

        # Randomize initial coordinates
        ic1, ic2, ic3, ip1, ip2 = \
            np.unravel_index(icell, markerdist.distribution().shape)
        def randomize(edges, idx):
            """Picks a random value between [edges[idx+1], edges[idx]]
            """
            return edges[idx] \
                + (edges[idx+1] - edges[idx]) * np.random.rand(idx.size,)

        if set(['r', 'phi', 'z']).issubset(markerdist.abscissae):
            mrk["r"]   = randomize(markerdist.abscissa_edges("r"),   ic1)
            mrk["phi"] = randomize(markerdist.abscissa_edges("phi"), ic2)
            mrk["z"]   = randomize(markerdist.abscissa_edges("z"),   ic3)
        elif set(['rho', 'theta', 'phi']).issubset(markerdist.abscissae):
            rhos       = randomize(markerdist.abscissa_edges("rho"),   ic1)
            thetas     = randomize(markerdist.abscissa_edges("theta"), ic2)
            mrk["phi"] = randomize(markerdist.abscissa_edges("phi"),   ic3)

            mrk["r"], mrk["z"] = self._ascot.input_rhotheta2rz(rhos, thetas, \
                                        mrk["phi"], 0*unyt.s)
        else:
            raise ValueError("Spatial abscissae basis not found form the "\
                        "distribution.")

        # Both (ppar, pperp) and (ekin, pitch) are accepted
        if "ppar" in markerdist.abscissae:
            ppa = randomize(markerdist.abscissa_edges("ppar"),  ip1)
            ppe = randomize(markerdist.abscissa_edges("pperp"), ip2)

            if mode == 'gc':
                pnorm = np.sqrt(ppa**2 + ppe**2)
                mrk["energy"] = physlib.energy_momentum(mass, pnorm).to("eV")
                mrk["pitch"]  = ppa / pnorm
                mrk['zeta']   = 2 * np.pi * np.random.rand(mrk['n'])
            elif mode == 'prt':
                br, bphi, bz = self._ascot.input_eval(
                    mrk['r'], mrk['phi'], mrk['z'], mrk['time'],
                    'br', 'bphi', 'bz')
                bhat = np.array([br, bphi, bz]) \
                    / np.sqrt(br**2 + bphi**2 + bz**2).v
                e1 = np.zeros(bhat.shape)
                e1[2,:] = 1
                e2 = np.cross(bhat.T, e1.T).T
                e1 = e2 / np.sqrt(np.sum(e2**2, axis=0))
                e2 = np.cross(bhat.T, e1.T).T

                zeta = 2 * np.pi * np.random.rand(mrk['n'])
                perphat = -np.sin(zeta)*e1-np.cos(zeta)*e2
                pvec = bhat * ppa + perphat * ppe
                mrk['vr']   = pvec[0,:] / mrk['mass']
                mrk['vphi'] = pvec[1,:] / mrk['mass']
                mrk['vz']   = pvec[2,:] / mrk['mass']

        else:
            if mode == 'gc':
                mrk["energy"] = \
                    randomize(markerdist.abscissa_edges("ekin"), ip1)
                mrk["pitch"]  = \
                    randomize(markerdist.abscissa_edges("pitch"), ip2)
                mrk['zeta']   = 2 * np.pi * np.random.rand(mrk['n'])
            elif mode == 'prt':
                energy = randomize(markerdist.abscissa_edges("ekin"), ip1)
                pitch  = randomize(markerdist.abscissa_edges("pitch"), ip2)
                gamma = physlib.gamma_energy(mass, energy)
                pnorm = physlib.pnorm_gamma(mass, gamma)
                ppa = pitch * pnorm
                ppe = np.sqrt(1.0 - pitch**2) * pnorm
                br, bphi, bz = self._ascot.input_eval(
                    mrk['r'], mrk['phi'], mrk['z'], mrk['time'],
                    'br', 'bphi', 'bz')
                bhat = np.array([br, bphi, bz]) \
                    / np.sqrt(br**2 + bphi**2 + bz**2).v
                e1 = np.zeros(bhat.shape)
                e1[2,:] = 1
                e2 = np.cross(bhat.T, e1.T).T
                e1 = e2 / np.sqrt(np.sum(e2**2, axis=0))
                e2 = np.cross(bhat.T, e1.T).T

                zeta = 2 * np.pi * np.random.rand(mrk['n'])
                perphat = -np.sin(zeta)*e1-np.cos(zeta)*e2
                pvec = bhat * ppa + perphat * ppe
                mrk['vr']   = pvec[0,:] / mrk['mass']
                mrk['vphi'] = pvec[1,:] / mrk['mass']
                mrk['vz']   = pvec[2,:] / mrk['mass']

        if not return_dists:
            return mrk

        # Create output distributions
        vol = markerdist.phasespacevolume()
        mrkdist = markerdist._copy()
        d = mrkdist._distribution.ravel().v * 0
        np.add.at(d, icell, 1)
        mrkdist._distribution = d.reshape(vol.shape) / vol.units

        prtdist = particledist._copy()
        d = prtdist._distribution.ravel().v * 0
        np.add.at(d, icell, mrk["weight"].v)
        prtdist._distribution = d.reshape(vol.shape) / vol

        return mrk, mrkdist, prtdist

    def rhoto5d(self, rho, prob, r_edges, phi_edges, z_edges, mom1_edges,
                mom2_edges):
        """Maps a 1D profile to a normalized 5D dist.

        The resulting probability density has only spatial variance.

        Parameters
        ----------
        rho : array_like, (nrho,)
            The radial rho grid where the profile is given.
        prob : array_like, (nrho,)
            Distribution value at the grid centers.
        r_edges : array_like
            R abscissa edges for the output distribution.
        phi_edges : array_like
            Phi abscissa edges for the output distribution.
        z_edges : array_like
            Z abscissa edges for the output distribution.
        mom1_edges : array_like
            Either ppar or ekin abscissa edges for the output distribution.
        mom2_edges : array_like
            Either pperp or pitch abscissa edges for the output distribution.

        Returns
        -------
        markerdist : :class:`Dist`
            Normalized 5D distribution from which markers can be sampled.
        """
        d = np.zeros((r_edges.size-1, phi_edges.size-1, z_edges.size-1,
                      mom1_edges.size-1, mom2_edges.size-1))
        if mom1_edges.units == mom2_edges.units:
            dist = DistData(d, r=r_edges, phi=phi_edges, z=z_edges,
                            ppar=mom1_edges, pperp=mom2_edges)
        else:
            dist = DistData(d, r=r_edges, phi=phi_edges, z=z_edges,
                            ekin=mom1_edges, pitch=mom2_edges)
        rhorpz = self._ascot.input_eval(
            dist.abscissa("r"), dist.abscissa("phi"),
            dist.abscissa("z"), 0*unyt.s, "rho", grid=True)[:,:,:,0]

        # To avoid private plasma we have to set everythin outside separatrix
        # to zero
        r, z = self._ascot.input_rhotheta2rz(
            np.ones((360,)), np.linspace(0,2*np.pi,360)*unyt.rad, 0*unyt.rad,
            0*unyt.s)
        R,P,Z = np.meshgrid(dist.abscissa("r"), dist.abscissa("phi"),
                            dist.abscissa("z"), indexing="ij")

        rhorpz = np.interp(rhorpz, rho, prob, left=0.0, right=0.0)
        rhorpz[Z < np.amin(z)] = 0.0
        rhorpz[Z > np.amax(z)] = 0.0
        rhorpz /= np.sum(rhorpz)
        rhorpz = np.tile(rhorpz.T,(mom2_edges.size-1,mom1_edges.size-1,1,1,1)).T
        dist._distribution = rhorpz / dist.phasespacevolume().units

        return dist

def gen_markerdist(a5, mass, charge, energy, rgrid, zgrid, kgrid,
                   rhoksidist=None, plot=False):
    """
    Generate marker (r,zksi) distribution from a given distribution.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        rgrid : array_like (n,1) <br>
            R coordinate grid.
        zgrid : array_like (n,1) <br>
            z coordinate grid.
        kgrid : array_like (n,1) <br>
            ksi coordinate grid.
        rhoksidist : tuple, optional <br>
            Tuple (rhogrid, ksigrid, weights) if markers are to be generated
            based weighted (rho_omp, ksi_omp) distribution. See phasespace.py.

    Returns:
        Tuple (rgrid,zgrid,ksigrid,markerdist) where first three are the grid
        axes (same as was given) and fourth element is the marker distribution
        normalized to one.
    """
    free = False
    if isinstance(a5, str):
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    # Marker distribution is generated with a Monte Carlo method by generating
    # chuncks of markers until the change in each cell is less than the
    # tolerance.
    Ndraw = 100000
    tol = 1e-1
    err = tol+1
    newdist = np.zeros((rgrid.size-1, zgrid.size-1, kgrid.size-1))
    probdensity = np.append( [0], np.cumsum(rhoksidist[2].ravel()) )
    while(err > tol):
        olddist = np.copy(newdist)

        if rhoksidist is not None:
            rhogrid = rhoksidist[0]
            ksigrid = rhoksidist[1]

            rhovals, ksivals = np.meshgrid(rhogrid[:-1], ksigrid[:-1],
                                           indexing='ij')

            prob = np.random.rand(Ndraw)
            idx  = np.searchsorted(probdensity, prob, side='left')-1

            weights = np.ones(idx.shape)
            rhovals = rhovals.ravel()[idx] + np.random.rand(Ndraw) \
                      * (rhogrid[1]-rhogrid[0])
            ksivals = ksivals.ravel()[idx] + np.random.rand(Ndraw) \
                      * (ksigrid[1]-ksigrid[0])

            newdist += phasespace.maprhoksi2rzk(a5, mass, charge, energy,
                                                rhovals, ksivals,
                                                rgrid, zgrid, kgrid,
                                                weights=weights)

        if np.sum(olddist.ravel()) > 0:
            ratio = np.absolute(
                ( newdist.ravel() / np.sum(newdist.ravel()) )
                / ( olddist.ravel() / np.sum(olddist.ravel())) )
            err = np.mean(np.absolute( ratio[np.isfinite(ratio)] - 1))
            print("Iteration error: " + str(err) + " Tolerance: " + str(tol))

    if free:
        a5.free(bfield=True)


    # Normalize
    newdist /= np.sum(newdist.ravel())

    if plot:
        plotdist((rgrid, zgrid, kgrid, newdist))

    return (rgrid, zgrid, kgrid, newdist)


def gen_DTalphadist(a5, rgrid, zgrid, ksigrid, plot=False):
    """
    Generate DT fusion alpha particle distribution based on NRL formula.

    This function assumes the HDF5 file has a plasma where first two ion species
    are Deuterium and Tritium.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        rgrid : array_like (n,1) <br>
            R coordinate grid.
        zgrid : array_like (n,1) <br>
            z coordinate grid.
        ksigrid : array_like (n,1) <br>
            ksi coordinate grid.

    Returns:
        Tuple (rgrid,zgrid,ksigrid,P_DT) where first three are the grid axes
        (same as was given) and fourth element is the number of alpha particles
        (E = 3.5 MeV) born in second in each cell.
    """
    free = False
    if isinstance(a5, str):
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True, plasma=True)

    # Prepare grid
    dr = (rgrid[1] - rgrid[0])
    dz = (zgrid[1] - zgrid[0])
    rvals = rgrid[:-1] + dr / 2
    zvals = zgrid[:-1] + dz / 2
    ksivals = ksigrid[:-1] + (ksigrid[1] - ksigrid[0]) / 2
    R, Z, Ksi = np.meshgrid(rvals, zvals, ksivals, indexing='ij')

    # Evaluate ion temperature and D and T density
    ti   = a5.evaluate(R=R, phi=0, z=Z, t=0, quantity="ti1").reshape(
        rgrid.size-1, zgrid.size-1, ksigrid.size-1) / (1e3*const.e)
    n1   = a5.evaluate(R=R, phi=0, z=Z, t=0, quantity="ni1").reshape(
        rgrid.size-1, zgrid.size-1, ksigrid.size-1)
    n2   = a5.evaluate(R=R, phi=0, z=Z, t=0, quantity="ni2").reshape(
        rgrid.size-1, zgrid.size-1, ksigrid.size-1)

    # Fusion alpha power density (NRL)
    sigmaDT = 3.68e-12 * np.power(ti, -2.0/3) \
              * np.exp(-19.94*np.power(ti, -1.0/3) )
    PDT     = 5.6e-13*n1*n2*sigmaDT / 1e6 # W / m^3

    # Convert power density to total alpha power in each cell
    PDT = PDT * 2*np.pi*R * dr * dz / (ksigrid.size-1) # W/bin

    # Convert power to particle birth rate and replace NaNs with zeros.
    dist = PDT / (3.5e6*const.e)
    dist[np.isnan(dist)] = 0

    if free:
        a5.free(bfield=True, plasma=True)

    if plot:
        plotdist((rgrid, zgrid, ksigrid, dist))

    return (rgrid, zgrid, ksigrid, dist)
