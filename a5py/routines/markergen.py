"""Generate marker populations.
"""
import numpy as np
import unyt

from a5py.ascot5io.marker import Marker
from a5py.ascot5io.dist import DistData
from a5py import physlib
from a5py.physlib.units import parseunits
from findiff import Diff
import warnings
from scipy.interpolate import CloughTocher2DInterpolator
import matplotlib.pyplot as plt
from scipy.integrate import simpson

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
            rejected = weight <= minweight
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

    @parseunits(density='m**-3', temperature='keV', mass='amu',
                charge='e', anum='1', znum='1', rhopol='1',
                strip=False)
    def from_maxwellian(self, nmrk, mass, charge, anum, znum,
                        rhopol, density, temperature,
                        mode, vth_cutoff: float=2.0):
        r"""
        Generate markers from a Maxwellian distribution function.

        Provided the density and temperatures, this will generate the
        markers that sample said distribution function that is to be
        described as:
        .. math::

            f(\vec{r}, \vec{v}) = n(\rho_{pol})\left(\frac{m}{2\pi kT}\right)^{3/2}
            \exp\left(-\frac{mv^2}{2kT(\rho_{pol})}\right)
        
        where :math:`\vec{r}` is the position vector, :math:`\vec{v}` is the
        velocity vector, :math:`n(\rho_{pol})` is the density at the given
        :math:`\rho_{pol}`, :math:`m` is the mass of the particle, :math:`k`
        is the Boltzmann constant, and :math:`T` is the temperature of the
        distribution. 
        
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
        rhopol : float
            The normalized poloidal radius where the density and temperature
            are given.
        density : float
            The density of the distribution at the given :math:`\rho_{pol}`.
        temperature : float
            The temperature of the distribution at the given
            :math:`\rho_{pol}`.
        mode : {'prt', 'gc'}
            Decides whether the distributions (and returned markers) are in
            particle or guiding center phase-space.
        vth_cutoff : float, optional
            The cutoff value for the thermal velocity. This is used to
            determine the maximum velocity of the particles. The default
            value is 2.0, which means that the maximum velocity is 2 times
            the thermal velocity.
        vth_sampling : float, optional
            Serves as the cut-off to separate between the importance sampling
            and the uniform sampling. The default value is np.inf, which means
            that the importance sampling is used for all the particles. If
            the value is set to 0.0, then the uniform sampling is used for all
            the particles. The value is the multiplier of the local thermal
            velocity.
        
        Returns
        -------
        mrk : dict
            Marker input data on a dictionary.       
        """
        # Checking the input parameters
        if mode not in ['prt', 'gc']:
            raise ValueError("mode must be either 'prt' or 'gc'.")
        if nmrk <= 0:
            raise ValueError("nmrk must be a positive integer.")
        if anum.value <= 0:
            raise ValueError("anum must be a positive integer.")
        if znum.value <= 0:
            raise ValueError("znum must be a positive integer.")
        if mass.value <= 0:
            raise ValueError("mass must be a positive number.")
        if vth_cutoff <= 0:
            raise ValueError("vth_cutoff must be a positive number.")
        
        # Checking that the density and temperature are positive
        density = np.atleast_1d(density)
        temperature = np.atleast_1d(temperature)
        rhopol = np.atleast_1d(rhopol)
    
        if np.any(density.value <= 0):
            raise ValueError("Density must be a positive number.")
        if np.any(temperature.value <= 0):
            raise ValueError("Temperature must be a positive number.")
        if np.any(rhopol.value < 0):
            raise ValueError("rhopol must be a positive number.")
        
        if rhopol.size != density.size or rhopol.size != temperature.size:
            raise ValueError("rhopol, density and temperature must have the same size.")
        
        # Getting the separatrix from the input
        r, z = self._ascot.input_rhotheta2rz(
            np.ones((360,)), np.linspace(0,2*np.pi, 360)*unyt.rad, 0*unyt.rad,
            0*unyt.s)
        
        rmax = r.max()
        zmax = z.max()
        rmin = r.min()
        zmin = z.min()

        # Initialize arrays for the coordinates and flags
        rvals = np.zeros(nmrk) * unyt.m
        phivals = np.zeros(nmrk) * unyt.rad
        zvals = np.zeros(nmrk) * unyt.m
        flags = np.ones(nmrk, dtype=bool)

        while np.any(flags):
            # Generate random values for the (R, phi, z) coordinates
            # rvals[flags] = np.random.uniform(rmin, rmax, np.sum(flags)) * unyt.m
            phivals[flags] = np.random.uniform(0, 2 * np.pi, np.sum(flags)) * unyt.rad
            zvals[flags] = np.random.uniform(zmin, zmax, np.sum(flags)) * unyt.m

            # To sample the radial coordinates, we need to include the 
            # Jacobian, which is the radial coordinate itself.
            uniform = np.random.uniform(0, 1, np.sum(flags))
            rvals[flags] = np.sqrt(rmin**2 + uniform * (rmax**2 - rmin**2))

            # Obtain the rhopol at which the particles are located
            rho = self._ascot.input_eval(rvals, phivals, zvals, 0 * unyt.s,
                         "rho", grid=False).value

            # Get the density and temperature for the given rho
            idens = np.interp(rho, rhopol, density.value) * density.units
            itemp = np.interp(rho, rhopol, temperature.value) * temperature.units

            # Update flags for particles with zero density or temperature
            flags = (idens == 0) | (itemp.to('eV').value == 0) | (rho > 1.05)

        # Computing the thermal velocity
        vth = np.sqrt(itemp / mass).to('m/s')

        mrk = dict(r=rvals, phi=phivals, z=zvals, 
                   vr=np.zeros(nmrk) * unyt.m / unyt.s,
                   vphi=np.zeros(nmrk) * unyt.m / unyt.s,
                   vz=np.zeros(nmrk) * unyt.m / unyt.s,
                   anum=anum, znum=znum, mass=mass,
                   charge=charge, time=np.zeros(nmrk) * unyt.s,
                   weight=np.ones(nmrk))
        
        # Generating the velocities. Uniform distribution to conver
        # properly tails of the distribution.
        mrk['vr'] = np.random.normal(0, 1, nmrk) * vth
        mrk['vphi'] = np.random.normal(0, 1, nmrk) * vth
        mrk['vz'] = np.random.normal(0, 1, nmrk) * vth
        vabs2 = mrk['vr']**2 + mrk['vphi']**2 + mrk['vz']**2

        # Generating the weights
        weights = idens #* (mass / (2 * np.pi * itemp))**(3/2) \
                        #* np.exp(-vabs2 / vth**2)
        # We need to normalize the weights so that they sum to 1.
        weights /= np.sum(weights) # Normalize the weights to sum to 1

        # We normalize the weights so they represent the number of 
        # particles in the distribution.
        rgrid = np.linspace(rmin, rmax, 512)
        zgrid = np.linspace(zmin, zmax, 513)
        rhorz = self._ascot.input_eval(
            rgrid, 0 * unyt.rad, zgrid, 0 * unyt.s,
            "rho", grid=True)
        dnsrz = np.interp(rhorz, rhopol, density.value).squeeze() * unyt.m**-3
        
        # We now integrate the dnsrz over the grid to get the
        # number of particles in the distribution.
        dr = rgrid[1] - rgrid[0]
        dz = zgrid[1] - zgrid[0]
        Nphys = 2 * np.pi * simpson(simpson(dnsrz * rgrid[:, None], zgrid), rgrid)
        weights *= Nphys
        mrk['weight'] = weights #* unyt.particles / (unyt.m**5 * unyt.rad / unyt.s**3)
        
        # Getting the guiding center coordinates.
        mrkgc = self._ascot.prt2gc(rprt=mrk['r'], phiprt=mrk['phi'], zprt=mrk['z'],
                                   vr=mrk['vr'], vphi=mrk['vphi'], vz=mrk['vz'],
                                   mass=mass, charge=charge
                                  )
        
        # Getting the magnetic field at the guiding center position.
        bnorm = self._ascot.input_eval(
            mrkgc['rgc'], mrkgc['phigc'], mrkgc['zgc'], 0 * unyt.s,
            'bnorm', grid=False)
        
        mrkgc['pperpgc'] = np.sqrt(2 * mrkgc['mugc'] * bnorm * mass).to('kg*m/s')

        if mode.lower() == 'gc':
            mrkout = dict(r=mrkgc['rgc'], phi=mrkgc['phigc'], z=mrkgc['zgc'],
                          ppara=mrkgc['pparagc'], pperp=mrkgc['pperpgc'],
                          mass=mass, charge=charge, anum=anum, znum=znum,
                          weight=mrk['weight'], time=mrk['time'])
        elif mode.lower() == 'prt':
            mrkout = mrk.copy()

        # We have to transform anum, znum, mass and charge to arrays.
        mrkout['anum'] = np.full(nmrk, anum.value) * anum.units
        mrkout['znum'] = np.full(nmrk, znum.value) * znum.units
        mrkout['mass'] = np.full(nmrk, mass.value) * mass.units
        mrkout['charge'] = np.full(nmrk, charge.value) * charge.units

        
        # We generate the distribution.
        pperpmax = np.nanmax(mrkgc['pperpgc']) * 1.1
        pperpmin = np.nanmin(mrkgc['pperpgc']) * 0.9
        pparmax = np.nanmax(mrkgc['pparagc']) * 1.1
        pparmin = np.nanmin(mrkgc['pparagc']) * 1.1
        npara = 64
        nperp = 65
        nr = 51
        nz = 52
        nphi = 1

        r_edges = np.linspace(rmin*0.95, rmax*1.05, nr+1)
        z_edges = np.linspace(zmin*1.05, zmax*1.05, nz+1)
        phi_edges = np.linspace(0, 2 * np.pi, nphi+1) * unyt.rad
        ppar_edges = np.linspace(pparmin, pparmax, npara+1)
        pperp_edges = np.linspace(pperpmin, pperpmax, nperp+1)

        # Getting the distribution
        d = np.zeros((nr, nphi, nz, npara, nperp))
        idx_r = np.digitize(mrkgc['rgc'].to('m').value, r_edges.to('m').value) - 1
        idx_z = np.digitize(mrkgc['zgc'].to('m').value, z_edges.to('m').value) - 1
        idx_phi = np.zeros_like(idx_r)
        idx_ppar = np.digitize(mrkgc['pparagc'].to('kg*m/s').value, ppar_edges.to('kg*m/s').value) - 1
        idx_pperp = np.digitize(mrkgc['pperpgc'].to('kg*m/s').value, pperp_edges.to('kg*m/s').value) - 1

        idx_flags = (idx_r >= 0) & (idx_r < nr) & \
                     (idx_z >= 0) & (idx_z < nz) & \
                     (idx_phi >= 0) & (idx_phi < nphi) & \
                     (idx_ppar >= 0) & (idx_ppar < npara) & \
                     (idx_pperp >= 0) & (idx_pperp < nperp)
        idx_r = idx_r[idx_flags]
        idx_z = idx_z[idx_flags]
        idx_phi = idx_phi[idx_flags]
        idx_ppar = idx_ppar[idx_flags]
        idx_pperp = idx_pperp[idx_flags]
        w0 = mrk['weight'][idx_flags]

        d[idx_r, idx_phi, idx_z, idx_ppar, idx_pperp] += w0

        coords = {'r': r_edges, 'phi': phi_edges, 'z': z_edges,
                 'ppar': ppar_edges, 'pperp': pperp_edges}
        dist = DistData(d, **coords)

        return mrkout, dist

    @parseunits(density='m**-3', temperature='keV', mass='amu',
                charge='e', anum='1', znum='1', rhopol='1',
                strip=False)
    def from_maxwellian_source(self, nmrk, mass, charge, anum, znum,
                               rhopol, density, temperature,
                               mode, vth_cutoff: float=2.0):
        r"""
        Generate markers from a Maxwellian distribution function source.

        Provided the density and temperatures, this will generate the
        markers that sample said distribution function that is to be
        described as:
        .. math::

            f(\vec{r}, \vec{v}) = n(\rho_{pol})\left(\frac{m}{2\pi kT}\right)^{3/2}
            \exp\left(-\frac{mv^2}{2kT(\rho_{pol})}\right)
        
        where :math:`\vec{r}` is the position vector, :math:`\vec{v}` is the
        velocity vector, :math:`n(\rho_{pol})` is the density at the given
        :math:`\rho_{pol}`, :math:`m` is the mass of the particle, :math:`k`
        is the Boltzmann constant, and :math:`T` is the temperature of the
        distribution.

        This is different from the `from_maxwellian` method in that in the sense
        that this generate markers that represent a source of particles, i.e., fluxes
        in the phase-space, enabling slowing down calculations with the time-independent
        structure of ASCOT5.
        
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
        rhopol : float
            The normalized poloidal radius where the density and temperature
            are given.
        density : float
            The density of the distribution at the given :math:`\rho_{pol}`.
        temperature : float
            The temperature of the distribution at the given
            :math:`\rho_{pol}`.
        mode : {'prt', 'gc'}
            Decides whether the distributions (and returned markers) are in
            particle or guiding center phase-space.
        vth_cutoff : float, optional
            The cutoff value for the thermal velocity. This is used to
            determine the maximum velocity of the particles. The default
            value is 2.0, which means that the maximum velocity is 2 times
            the thermal velocity.
        vth_sampling : float, optional
            Serves as the cut-off to separate between the importance sampling
            and the uniform sampling. The default value is np.inf, which means
            that the importance sampling is used for all the particles. If
            the value is set to 0.0, then the uniform sampling is used for all
            the particles. The value is the multiplier of the local thermal
            velocity.
        
        Returns
        -------
        mrk : dict
            Marker input data on a dictionary.       
        """
        # Checking the input parameters
        if mode not in ['prt', 'gc']:
            raise ValueError("mode must be either 'prt' or 'gc'.")
        if nmrk <= 0:
            raise ValueError("nmrk must be a positive integer.")
        if anum.value <= 0:
            raise ValueError("anum must be a positive integer.")
        if znum.value <= 0:
            raise ValueError("znum must be a positive integer.")
        if mass.value <= 0:
            raise ValueError("mass must be a positive number.")
        if vth_cutoff <= 0:
            raise ValueError("vth_cutoff must be a positive number.")
        
        # Checking that the density and temperature are positive
        density = np.atleast_1d(density)
        temperature = np.atleast_1d(temperature)
        rhopol = np.atleast_1d(rhopol)

        if np.any(density.value <= 0):
            raise ValueError("Density must be a positive number.")
        if np.any(temperature.value <= 0):
            raise ValueError("Temperature must be a positive number.")
        if np.any(rhopol.value < 0):
            raise ValueError("rhopol must be a positive number.")
        
        if rhopol.size != density.size or rhopol.size != temperature.size:
            raise ValueError("rhopol, density and temperature must have the same size.")
        
        # Getting the separatrix from the input
        r, z = self._ascot.input_rhotheta2rz(
            np.ones((360,)), np.linspace(0,2*np.pi, 360)*unyt.rad, 0*unyt.rad,
            0*unyt.s)
        
        
        bdata = self._ascot.data.bfield.active.read()
        
        rmax = bdata['rmax'][0]
        zmax = bdata['zmax'][0]
        rmin = bdata['rmin'][0]
        zmin = bdata['zmin'][0]

        # Initialize arrays for the coordinates and flags
        rvals = np.zeros(nmrk) * unyt.m
        phivals = np.zeros(nmrk) * unyt.rad
        zvals = np.zeros(nmrk) * unyt.m
        flags = np.ones(nmrk, dtype=bool)

        while np.any(flags):
            # Generate random values for the (R, phi, z) coordinates
            rvals[flags] = np.random.uniform(rmin, rmax, np.sum(flags)) * unyt.m
            phivals[flags] = np.random.uniform(0, 2 * np.pi, np.sum(flags)) * unyt.rad
            zvals[flags] = np.random.uniform(zmin, zmax, np.sum(flags)) * unyt.m

            # To sample the radial coordinates, we need to include the 
            # Jacobian, which is the radial coordinate itself.
            # uniform = np.random.uniform(0, 1, np.sum(flags))
            # rvals[flags] = np.sqrt(rmin**2 + uniform * (rmax**2 - rmin**2))

            # Obtain the rhopol at which the particles are located
            rho = self._ascot.input_eval(rvals, phivals, zvals, 0 * unyt.s,
                         "rho", grid=False).value

            # Get the density and temperature for the given rho
            idens = np.interp(rho, rhopol, density.value) * density.units
            itemp = np.interp(rho, rhopol, temperature.value) * temperature.units

            # Update flags for particles with zero density or temperature
            flags = (idens == 0) | (itemp == 0) # | (rho > 1.0)

        # Computing the thermal velocity
        vth = np.sqrt(itemp / mass).to('m/s')

        mrk = dict(r=rvals, phi=phivals, z=zvals, 
                   vr=np.zeros(nmrk) * unyt.m / unyt.s,
                   vphi=np.zeros(nmrk) * unyt.m / unyt.s,
                   vz=np.zeros(nmrk) * unyt.m / unyt.s,
                   anum=anum, znum=znum, mass=mass,
                   charge=charge, time=np.zeros(nmrk) * unyt.s,
                   weight=np.ones(nmrk))
        
        # Generating the velocities. Uniform distribution to conver
        # properly tails of the distribution.
        mrk['vr'] = np.random.normal(0, 1, nmrk) * vth
        mrk['vphi'] = np.random.normal(0, 1, nmrk) * vth
        mrk['vz'] = np.random.normal(0, 1, nmrk) * vth
        vabs2 = mrk['vr']**2 + mrk['vphi']**2 + mrk['vz']**2
        vabs = np.sqrt(vabs2)
        ienergy = 0.5 * mass * vabs2

        # ---- Generating the weights
        f0 = idens # The velocity dependence is included in the way the particles are sampled

        # Let's divide by the average volume in the grid.
        total_vol = 2 * np.pi * (rmax**2 - rmin**2) * (zmax - zmin) * unyt.m**3
        f0 *= total_vol / nmrk

        # We need to get the derivatives of the density and temperature
        # with respect to Psi.
        # The d/drho = d/dPsi * 1/(2* rhopol * PsiN)
        # where PsiN is the difference between the normalized poloidal
        # flux at the axis and at the separatrix.
        ddrho = Diff(0, acc=4)
        ddrho.set_grid({0: rhopol[1] - rhopol[0]})
        dT = ddrho(np.log(temperature)) / (2 * rhopol)
        dn = ddrho(np.log(temperature)) / (2 * rhopol)

        # We check whether the first rhopol is zero, in which case, 
        # we simply set the first value to the second value.
        if rhopol[0] == 0:
            dT[0] = dT[1]
            dn[0] = dn[1]

        # We now need to compute the change between the density and temperature
        psi0 = bdata['psi0']
        psisep = bdata['psi1']
        psiN = psisep - psi0
        dT /= psiN
        dn /= psiN

        # We need now evaluate the derivatives at the positions of the markers.
        tmp = self._ascot.input_eval(mrk['r'], mrk['phi'], mrk['z'], 0 * unyt.s,
                         "rho", "br", "bz", grid=False)
        irho = tmp[0]
        ibr = tmp[1]
        ibz = tmp[2]

        iln_n = np.interp(irho.value, rhopol, dn.value) / unyt.Wb
        iln_T = np.interp(irho.value, rhopol, dT.value) / unyt.Wb

        # The first contribution to the weight are the density gradients:
        w1 = iln_n + (ienergy/itemp - 1.5) * iln_T
        w1 *= rvals * (- ibz * mrk['vr'] + ibr * mrk['vz'])

        # The rest of the contributions come from the slowing down. To achieve
        # that we will evaluate the diffusion coefficients, Dpara, Dperp and K
        # in a grid (Psi, v) and compute the respective derivatives.
        # We get the values of R along the Theta=0, aka, the outer midplane.
        z = bdata['axisz'] * unyt.m
        R = np.linspace(bdata['axisr'], rmax, 500) * unyt.m
        v = np.linspace(vabs.min()*0.9, vabs.max()*1.1, 150)

        rho_on_R = self._ascot.input_eval(
            R, 0 * unyt.rad, z, 0 * unyt.s, 'rho', grid=False).value

        tmp = self._ascot.input_eval_collcoefs(mass, charge,
                                               R, 0*unyt.rad, z, 0*unyt.s,
                                               v, 'dpara', 'dperp', 'f')
        dpara = np.sum(tmp[0], axis=0)
        dperp = np.sum(tmp[1], axis=0)
        k = np.sum(tmp[2], axis=0) / v[None, :]

        # Getting derivatives along the velocity.
        dv = Diff(1, acc=4)
        dv.set_grid({1: v[1] - v[0]})
        d_dpara = dv(dpara)
        d_dperp = dv(dperp)
        d2_dpara = dv(d_dpara)
        dk = dv(k)

        # Evaluating the diffusion and drag coefficients at the
        # positions of the markers.
        grr, gvv = np.meshgrid(rho_on_R, v, indexing='ij')
        points = np.array([grr.ravel(), gvv.ravel()]).T
        dpara_intrp = CloughTocher2DInterpolator(points, dpara.ravel(), 
                                                 fill_value=0.0)
        dperp_intrp = CloughTocher2DInterpolator(points, dperp.ravel(), 
                                                 fill_value=0.0)
        k_intrp = CloughTocher2DInterpolator(points, k.ravel(), 
                                                 fill_value=0.0)
        d_dpara_intrp = CloughTocher2DInterpolator(points, d_dpara.ravel(), 
                                                 fill_value=0.0)
        d_dperp_intrp = CloughTocher2DInterpolator(points, d_dperp.ravel(), 
                                                 fill_value=0.0)
        d2_dpara_intrp = CloughTocher2DInterpolator(points, d2_dpara.ravel(), 
                                                 fill_value=0.0)
        dk_intrp = CloughTocher2DInterpolator(points, dk.ravel(), 
                                                 fill_value=0.0)
        
        dpara = dpara_intrp((irho.value, vabs.value)) * dpara.units
        dperp = dperp_intrp((irho.value, vabs.value)) * dperp.units
        k = k_intrp((irho.value, vabs.value)) * k.units
        dk = dk_intrp((irho.value, vabs.value)) * dk.units
        d_dpara = d_dpara_intrp((irho.value, vabs.value)) * d_dpara.units
        d_dperp = d_dperp_intrp((irho.value, vabs.value)) * d_dperp.units
        d2_dpara = d2_dpara_intrp((irho.value, vabs.value)) * d2_dpara.units


        # Drag contribution.
        # wdrag = v * dK/dv + (3 - E/T) * K
        wdrag = vabs * dk + (3 - ienergy / itemp) * k

        # Diffusion contributions.
        # 1. -1/v * (dDpara/dv + (Dperp - Dpara)/ (2*v))
        wdiff1 = -1/vabs * (d_dpara + (dperp - dpara) / (2 * vabs))

        # 2. -1/2 * d^2 Dpara/dv^2
        wdiff2 = -0.5 * d2_dpara

        # 3. -1/(2v) * d(Dperp - Dpara)/dv
        wdiff3 = -0.5 / vabs * (d_dperp - d_dpara)

        # 4. E/T * 1/v * (dDpara/dv + (Dperp - Dpara)/v)
        wdiff4 = (ienergy / itemp) / vabs * (d_dpara + (dperp - dpara) / vabs)

        # 5. - m/(2T) * (2E/T D_para - D_para + 2D_perp)
        wdiff5 = - mass / (2 * itemp) * (2 * ienergy / itemp * dpara - dpara + 2 * dperp)

        # We now compute the total weight.
        w2 = (wdrag + wdiff1 + wdiff2 + wdiff3 + wdiff4 + wdiff5)

        weights = f0 * (w1 + w2) # The weights are the product of the density and the contributions

        # Now, the tricky part: the weight normalization. We assume that
        # the amount of particles are sufficiently large that its integral 
        # represents the total flux in the phase-space.
        wN = np.sum(weights) # Should be close to zero
        print(f"Total weight: {wN:.3e} particles/s")
        # weights /= wN # Normalizing the weights to the total flux
        mrk['weight'] = weights

        # ---- Getting the guiding center coordinates.
        mrkgc = self._ascot.prt2gc(rprt=mrk['r'], phiprt=mrk['phi'], zprt=mrk['z'],
                                   vr=mrk['vr'], vphi=mrk['vphi'], vz=mrk['vz'],
                                   mass=mass, charge=charge
                                  )
        
        # Getting the magnetic field at the guiding center position.
        bnorm = self._ascot.input_eval(
            mrkgc['rgc'], mrkgc['phigc'], mrkgc['zgc'], 0 * unyt.s,
            'bnorm', grid=False)
        
        mrkgc['pperpgc'] = np.sqrt(2 * mrkgc['mugc'] * bnorm * mass).to('kg*m/s')

        if mode.lower() == 'gc':
            mrkout = dict(r=mrkgc['rgc'], phi=mrkgc['phigc'], z=mrkgc['zgc'],
                          ppar=mrkgc['pparagc'], pperp=mrkgc['pperpgc'],
                          mass=mass*np.ones(nmrk), 
                          charge=charge*np.ones(nmrk), 
                          anum=anum*np.ones(nmrk), 
                          znum=znum*np.ones(nmrk),
                          weight=weights, 
                          time=mrk['time'])
        elif mode.lower() == 'prt':
            mrkout = mrk.copy()
        
        # We generate the distribution.
        pperpmax = np.nanmax(mrkgc['pperpgc']) * 1.1
        pperpmin = np.nanmin(mrkgc['pperpgc']) * 0.9
        pparmax = np.nanmax(mrkgc['pparagc']) * 1.1
        pparmin = np.nanmin(mrkgc['pparagc']) * 1.1
        npara = 64
        nperp = 65
        nr = 31
        nz = 32
        nphi = 1

        r_edges = np.linspace(rmin*0.95, rmax*1.05, nr+1) * unyt.m
        z_edges = np.linspace(zmin*1.05, zmax*1.05, nz+1) * unyt.m
        phi_edges = np.linspace(0, 2 * np.pi, nphi+1) * unyt.rad
        ppar_edges = np.linspace(pparmin, pparmax, npara+1)
        pperp_edges = np.linspace(pperpmin, pperpmax, nperp+1)

        # Getting the distribution
        d = np.zeros((nr, nphi, nz, npara, nperp))
        idx_r = np.digitize(mrkgc['rgc'].to('m').value, r_edges.to('m').value) - 1
        idx_z = np.digitize(mrkgc['zgc'].to('m').value, z_edges.to('m').value) - 1
        idx_phi = np.zeros_like(idx_r)
        idx_ppar = np.digitize(mrkgc['pparagc'].to('kg*m/s').value, ppar_edges.to('kg*m/s').value) - 1
        idx_pperp = np.digitize(mrkgc['pperpgc'].to('kg*m/s').value, pperp_edges.to('kg*m/s').value) - 1

        idx_flags = (idx_r >= 0) & (idx_r < nr) & \
                    (idx_z >= 0) & (idx_z < nz) & \
                    (idx_ppar >= 0) & (idx_ppar < npara) & \
                    (idx_pperp >= 0) & (idx_pperp < nperp)
        
        idx_r = idx_r[idx_flags]
        idx_z = idx_z[idx_flags]
        idx_phi = idx_phi[idx_flags]
        idx_ppar = idx_ppar[idx_flags]
        idx_pperp = idx_pperp[idx_flags]
        # w0 = np.maximum(0.0, mrk['weight'][idx_flags])

        d[idx_r, idx_phi, idx_z, idx_ppar, idx_pperp] += mrk['weight'][idx_flags]

        coords = {'r': r_edges, 'phi': phi_edges, 'z': z_edges,
                 'ppar': ppar_edges, 'pperp': pperp_edges}
        dist = DistData(d, **coords)

        return mrkout, dist
        
