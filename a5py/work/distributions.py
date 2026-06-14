"""Tools for working with distributions."""

import numpy as np

class DistMoment:
    """Class that stores moments calculated from a distribution.
    """

    def __init__(self, x1, x2, x3, r, phi, z, area, volume, rhodist):
        """Initialize moment storage.

        The real space abscissa edges define where the moments are defined.
        The coordinates of the bin edges are required to calculate moments.
        Area and volume are needed for normalization and to transform from
        phase-space to real space.

        Parameters
        ----------
        x1 : array_like
            First real-space coordinate abscissa edges (R or rho).
        x2 : array_like
            First real-space coordinate abscissa edges (R or rho).
        x3 : array_like
            First real-space coordinate abscissa edges (R or rho).
        r : array_like
            R coordinates corresponding to bin centers.
        phi : array_like
            phi coordinates corresponding to bin centers.
        z : array_like
            z coordinates corresponding to bin centers.
        area : array_like
            Poloidal plane area of each bin.
        volume : array_like
            Volume of each bin.
        rhodist : bool
            Flag indicating whether the moment is in (rho,theta,phi) or
            (R,phi,z) coordinates.
        """
        if rhodist:
            self.rho   = x1
            self.theta = x2
            self.phi   = x3
        else:
            self.r     = x1
            self.phi   = x2
            self.z     = x3

        self.rc        = r
        self.phic      = phi
        self.zc        = z
        self.area      = area
        self.volume    = volume
        self.rhodist   = rhodist

    def ordinate(self, ordinate, toravg=False, polavg=False):
        """Return stored moment.

        Parameters
        ----------
        ordinate : str
            Name of the moment.
        toravg : bool, optional
            Return toroidal average of the ordinate.
        polavg : bool, optional
            Return poloidal average of the ordinate.

            Only valid if ``rhodist=True``. Both ``toravg`` and ``polavg`` can
            be set simultaneously, in which a radial profile is returned.

        Returns
        -------
        data : array_like
            Ordinate data.
        """
        name = "_ordinate_" + ordinate
        if not hasattr(self, name):
            raise ValueError("Unknown ordinate: " + ordinate)
        ordinate = getattr(self, name)

        if toravg:
            if self.rhodist:
                ordinate = np.sum(ordinate * self.volume, axis=2) \
                    / np.sum(self.volume, axis=2)
            else:
                ordinate = np.sum(ordinate * self.volume, axis=1) \
                    / np.sum(self.volume, axis=1)
        if polavg:
            if self.rhodist:
                volume = np.sum(self.volume, axis=2) if toravg else self.volume
                ordinate = np.sum(ordinate * volume, axis=1) \
                    / np.sum(volume, axis=1)
            else:
                raise ValueError("Cannot take poloidal average of non-rhodist")
        return ordinate

    def add_ordinates(self, **ordinates):
        """Add moments.

        Parameters
        ----------
        **ordinates : array_like
            Names and data for each ordinate to be added.
        """
        for ordinate, val in ordinates.items():
            if hasattr(self, "_ordinate_" + ordinate):
                raise ValueError("Ordinate %s is already present" % ordinate)
            setattr(self, "_ordinate_" + ordinate, val)

    def list_ordinates(self):
        """List all moments
        """
        ordinates = []
        for k in self.__dict__.keys():
            if "_ordinate_" in k:
                ordinates.append(k[10:])
        return ordinates

    def plot(self, ordinate, axes=None, cax=None, logscale=False):
        """Plot radial or (R,z) profile of a distribution moment.

        The plotted profile is the average of (theta, phi) or phi depending
        on whether the input is calculated from a rho distribution or not.

        Parameters
        ----------
        ordinate : str
            Name of the moment to be plotted.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        logscale: bool, optional
            Whether the plot is in logarithmic scale.
        """
        if self.rhodist:
            ylabel = ordinate
            ordinate = self.ordinate(ordinate, toravg=True, polavg=True)
            ylabel += " [" + str(ordinate.units) + "]"
            a5plt.mesh1d(self.rho, ordinate,
                         xlabel="Normalized poloidal flux",
                         ylabel=ylabel, axes=axes, logscale=logscale)
        else:
            clabel = ordinate
            ordinate = self.ordinate(ordinate, toravg=True)
            clabel += " [" + str(ordinate.units) + "]"
            a5plt.mesh2d(self.r, self.z, ordinate, axesequal=True,
                         xlabel="R [m]", ylabel="z [m]", clabel=clabel,
                         axes=axes, cax=cax, logscale=logscale)

class Dist():

    def get(self):
        """Return the distribution data.

        Returns
        -------
        dist : class:`DistData`
            Distribution data.
        """
        with self as f:
            histogram = np.sum(f["ordinate"][:], axis=0) * unyt.particles
            abscissa_edges = {}
            for i in range(int(f["abscissa_ndim"][:])):
                abscissa = f["abscissa_vec_0"+str(i+1)]
                name     = abscissa.attrs["name_0"+str(i)].decode("utf-8")
                unit     = abscissa.attrs["unit_0"+str(i)].decode("utf-8")
                try:
                    abscissa_edges[name] = abscissa[:] * unyt.Unit(unit)
                except:
                    unit = unit.replace(" ", "*")
                    abscissa_edges[name] = abscissa[:] * unyt.Unit(unit)

        return DistData(histogram, **abscissa_edges)

    @staticmethod
    def density(dist, moment):
        """Calculate number density.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        integrate = {}
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
        dist = dist.integrate(copy=True, **integrate)
        moment.add_ordinates(density=dist.histogram() / moment.volume)

    @staticmethod
    def chargedensity(dist, moment):
        """Calculate charge density.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        integrate = {}
        dist = dist.integrate(copy=True, charge=dist.abscissa("charge"))
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
        dist.integrate(**integrate)
        moment.add_ordinates(chargedensity=dist.histogram() / moment.volume)

    @staticmethod
    def energydensity(mass, dist, moment):
        """Calculate energy density.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        integrate = {}
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi", "pperp", "ppar"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "z", "phi", "pperp", "ppar"]:
                    integrate[k] = np.s_[:]
        dist = dist.integrate(copy = True, **integrate)
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"), dist.abscissa("pperp"))
        pnorm = np.sqrt(ppa.ravel()**2 + ppe.ravel()**2)
        ekin  = (physlib.gamma_momentum(mass, pnorm) - 1) * mass * unyt.c**2
        dist._multiply(ekin.reshape(ppa.shape).T, "ppar", "pperp")
        dist.integrate(ppar=np.s_[:], pperp=np.s_[:])
        moment.add_ordinates(
            energydensity=dist.histogram().to("J") / moment.volume)

    @staticmethod
    def pressure(mass, dist, moment):
        """Calculate pressure.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        integrate = {}
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi", "pperp", "ppar"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "z", "phi", "pperp", "ppar"]:
                    integrate[k] = np.s_[:]
        dist = dist.integrate(copy = True, **integrate)
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"),
                               dist.abscissa("pperp"))
        pnorm = np.sqrt(ppa.ravel()**2 + ppe.ravel()**2)
        vnorm = physlib.velocity_momentum(mass, pnorm).reshape(ppa.shape)
        vpa = ppa * vnorm / pnorm.reshape(ppa.shape)
        vpe = ppe * vnorm / pnorm.reshape(ppa.shape)
        d = dist._copy()
        d._multiply(vpa.T, "ppar", "pperp")
        d.integrate(ppar=np.s_[:], pperp=np.s_[:])
        upa = d.histogram()
        d = dist._copy()
        d.integrate(ppar=np.s_[:], pperp=np.s_[:])
        n = d.histogram()
        upa /= n
        upa[n==0] = 0
        d = dist._copy()
        d._multiply(vpa.T**2, "ppar", "pperp")
        d.integrate(ppar=np.s_[:], pperp=np.s_[:])
        vpa2 = d.histogram()
        d = dist._copy()
        d._multiply(vpa.T, "ppar", "pperp")
        d.integrate(ppar=np.s_[:], pperp=np.s_[:])
        vpa = d.histogram()
        Ppa = mass*(vpa2 - 2*vpa*upa + upa**2) / 3
        d = dist._copy()
        d._multiply(vpe.T**2, "ppar", "pperp")
        d.integrate(ppar=np.s_[:], pperp=np.s_[:])
        vpe2 = d.histogram()
        d = dist._copy()
        d._multiply(vpe.T, "ppar", "pperp")
        d.integrate(ppar=np.s_[:], pperp=np.s_[:])
        vpe = d.histogram()
        Ppe = mass*vpe2 / 3
        moment.add_ordinates(
            pressure=np.sqrt(Ppa**2+2*Ppe**2).to("J") / moment.volume
        )

    @staticmethod
    def toroidalcurrent(ascot, mass, dist, moment):
        """Calculate toroidal current density.

        Parameters
        ----------
        ascot : :class:`Ascot`
            Ascot object for interpolating input data.
        mass : float
            Test particle mass.
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        dist = dist.integrate(copy=True, charge=dist.abscissa("charge"),
                              ppar=dist.abscissa("ppar"))
        bphi, bnorm = ascot.input_eval(
            moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm")
        bphi  = bphi.reshape(moment.volume.shape)
        bnorm = bnorm.reshape(moment.volume.shape)

        integrate = {}
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
        dist.integrate(**integrate)
        dist._distribution *= bphi / (bnorm * mass)
        moment.add_ordinates(
            toroidalcurrent=(dist.histogram() / moment.volume).to("A/m**2"))

    @staticmethod
    def parallelcurrent(ascot, mass, dist, moment, drive=False):
        """Calculate parallel current density.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        drive : bool, optional
            Calculate current drive instead.
        """
        Z = dist.abscissa("charge")[0] / unyt.e
        dist = dist.integrate(copy=True, charge=dist.abscissa("charge"),
                              ppar=dist.abscissa("ppar"))
        integrate = {}
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
        dist.integrate(**integrate)
        dist._distribution *= 1/mass

        if(drive):
            rho, zeff, ne, te, rminor = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0,
                "rho", "zeff", "ne", "te", "rminor"
                )
            clogee = 31.3 - np.log(np.sqrt(ne/unyt.m**3)/(te/unyt.eV))
            q, _, _, ftrap = ascot.input_eval_safetyfactor(
                rho, return_ftrap=True)
            aspectratio = np.sqrt(moment.rc/rminor)
            F = Z/zeff
            nu_eff = (  1.779e-55 * unyt.J**2 * unyt.m**2 * q * moment.rc * ne
                      * clogee / (te**2 * aspectratio**(-3/2)))
            X = ftrap / (1 + ( 1 - 0.1*ftrap ) * np.sqrt(nu_eff)
                           + 0.5 * ( 1 -  ftrap ) * nu_eff / zeff )
            zeffplusone = zeff+1
            G = (  ( 1 + 1.4 / zeffplusone ) * X
                 - (1.9 / zeffplusone) * X**2
                 + (0.3 / zeffplusone) * X**3
                 + (0.2 / zeffplusone) * X**4
                 )
            jpar = (dist.histogram() / moment.volume).to("A/m**2")
            moment.add_ordinates(currentdrive=(1 - F *(1 - G))*jpar)
        else:
            moment.add_ordinates(
                parallelcurrent=(dist.histogram() / moment.volume).to("A/m**2"))

    @staticmethod
    def powerdep(ascot, mass, dist, moment):
        """Calculate collisional power deposition to plasma.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        dist = dist._copy()
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"),
                               dist.abscissa("pperp"))
        pnorm = np.sqrt(ppa.ravel()**2 + ppe.ravel()**2)
        vnorm = physlib.velocity_momentum(mass, pnorm)
        for i, qa in enumerate(dist.abscissa("charge")):
            k = ascot.input_eval_collcoefs(
                mass, qa, moment.rc.ravel(), moment.phic.ravel(),
                moment.zc.ravel(), np.zeros(moment.rc.ravel().shape)*unyt.s,
                vnorm, "k", grid=True)
            k = -np.sum(k, axis=0) # Minus because k is from plasma to particle
            k = k.ravel().reshape(dist._distribution[:,:,:,:,:,i,0].shape)
            dist._distribution[:,:,:,:,:,i,0] *= k.v

        dist._distribution *= k.units * mass
        dist.integrate(charge=np.s_[:], time=np.s_[:])
        dist._multiply(vnorm.reshape(ppa.shape).T, "ppar", "pperp")
        dist.integrate(ppar=np.s_[:], pperp=np.s_[:])
        moment.add_ordinates(
            powerdep=(dist.histogram() / moment.volume ).to("W/m**3"))

    @staticmethod
    def electronpowerdep(ascot, mass, dist, moment):
        """Calculate power deposition to electrons.

        Parameters
        ----------
        ascot : :class:`Ascot`
            Ascot object for interpolating input data.
        mass : float
            Test particle mass.
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        dist = dist._copy()
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"),
                               dist.abscissa("pperp"))
        pnorm = np.sqrt(ppa.ravel()**2 + ppe.ravel()**2)
        vnorm = physlib.velocity_momentum(mass, pnorm)
        for i, qa in enumerate(dist.abscissa("charge")):
            k = ascot.input_eval_collcoefs(
                mass, qa, moment.rc.ravel(), moment.phic.ravel(),
                moment.zc.ravel(), np.zeros(moment.rc.ravel().shape)*unyt.s,
                vnorm, "k", grid=True)
            k = -k[0,:,:] # Minus because k is from plasma to particle
            k = k.ravel().reshape(dist._distribution[:,:,:,:,:,i,0].shape)
            dist._distribution[:,:,:,:,:,i,0] *= k.v

        dist._distribution *= k.units * mass
        dist.integrate(charge=np.s_[:], time=np.s_[:])
        dist._multiply(vnorm.reshape(ppa.shape).T, "ppar", "pperp")
        dist.integrate(ppar=np.s_[:], pperp=np.s_[:])
        moment.add_ordinates(
            electronpowerdep=(dist.histogram() / moment.volume ).to("W/m**3"))

    @staticmethod
    def ionpowerdep(ascot, mass, dist, moment):
        """Calculate power deposition to ions.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        dist = dist._copy()
        ppa, ppe = np.meshgrid(dist.abscissa("ppar"),
                               dist.abscissa("pperp"))
        pnorm = np.sqrt(ppa.ravel()**2 + ppe.ravel()**2)
        vnorm = physlib.velocity_momentum(mass, pnorm)
        for i, qa in enumerate(dist.abscissa("charge")):
            k = ascot.input_eval_collcoefs(
                mass, qa, moment.rc.ravel(), moment.phic.ravel(),
                moment.zc.ravel(), np.zeros(moment.rc.ravel().shape)*unyt.s,
                vnorm, "k", grid=True)
            k = -np.sum(k[1:], axis=0) # Minus because k is from plasma to prt
            k = k.ravel().reshape(dist._distribution[:,:,:,:,:,i,0].shape)
            dist._distribution[:,:,:,:,:,i,0] *= k.v

        dist._distribution *= k.units * mass
        dist.integrate(charge=np.s_[:], time=np.s_[:])
        dist._multiply(vnorm.reshape(ppa.shape).T, "ppar", "pperp")
        dist.integrate(ppar=np.s_[:], pperp=np.s_[:])
        moment.add_ordinates(
            ionpowerdep=(dist.histogram() / moment.volume ).to("W/m**3"))

    @staticmethod
    def jxbtorque(ascot, mass, dist, moment):
        """Calculate j_rad x B_pol toroidal torque.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        dist = dist._copy()
        gradbr, gradbphi, gradbz, curlbr, curlbphi, curlbz, br, bphi, bz = \
            ascot.input_eval(moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(), 0*unyt.s,
                             "gradbr", "gradbphi", "gradbz", "curlbr",
                             "curlbphi", "curlbz", "br", "bphi", "bz")
        bvec = unyt.unyt_array([br, bphi, bz]).T
        bnorm = np.sqrt(np.sum(bvec**2))
        curlb = unyt.unyt_array([gradbr, gradbphi, gradbz]).T
        gradb = unyt.unyt_array([curlbr, curlbphi, curlbz]).T
        curlbhat = ( np.cross(gradb, bvec/bnorm)*unyt.T/unyt.m - curlb ) / bnorm
        for iq, q in enumerate(dist.abscissa("charge")):
            for ippa, ppa in enumerate(dist.abscissa("ppar")):
                Bstar = bvec + (ppa / q) * curlbhat
                for ippe, ppe in enumerate(dist.abscissa("pperp")):
                    gamma = physlib.gamma_momentum(mass, np.sqrt(ppa**2+ppe**2))
                    dr = (ppa / (gamma*mass)) * Bstar[:,0] \
                        / (np.sum(Bstar * bvec / bnorm, axis=1))
                    dz = (ppa / (gamma*mass)) * Bstar[:,2] \
                        / (np.sum(Bstar * bvec / bnorm, axis=1))
                    deltapsi = bz * moment.rc.ravel() * dr - br * moment.rc.ravel() * dz
                    deltapsi = deltapsi.reshape(moment.volume.shape)
                    dist._distribution[:,:,:,ippa, ippe, iq, 0] *= (deltapsi * q).v
        integrate = {}
        dist._distribution *= (deltapsi * q).units
        if moment.rhodist:
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
        else:
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
        dist.integrate(**integrate)
        moment.add_ordinates(
            jxbtorque=(dist.histogram() / moment.volume).to("N*m/m**3"))

    @staticmethod
    def collTorque(ascot, mass, dist, moment):
        """Calculate power deposition to ions.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        if not "ekin" in dist.abscissae or not "pitch" in dist.abscissae:
            raise ValueError("Distribution must be in energy-pitch basis.")
        if moment.rhodist:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            pitch = dist.abscissa("pitch")
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = np.sum(coefs["K"],axis=0).reshape((*moment.volume.shape,va.size))
                nu = np.sum(coefs["nu"],axis=0).reshape((*moment.volume.shape,va.size))

                dpitch = -nu
                dPpara = mass*K + mass*va*dpitch
                units = dPpara.units
                dist._distribution[:,:,:,:,0,0] *= dPpara.v

            bphi, bnorm = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm")
            bphi   = bphi.reshape(moment.volume.shape)
            bnorm  = bnorm.reshape(moment.volume.shape)

            #dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist.integrate(ekin=np.s_[:], charge=np.s_[:], time=np.s_[:])
            dist._distribution *= (bphi/bnorm) *moment.rc *units
        else:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            pitch = dist.abscissa("pitch")
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = np.sum(coefs["K"],axis=0).reshape((*moment.volume.shape,va.size))
                nu = np.sum(coefs["nu"],axis=0).reshape((*moment.volume.shape,va.size))

                dpitch = -nu
                dPpara = mass*K + mass*va*dpitch
                units = dPpara.units
                dist._distribution[:,:,:,:,0,0] *= dPpara.v

            bphi, bnorm = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm")
            bphi   = bphi.reshape(moment.volume.shape)
            bnorm  = bnorm.reshape(moment.volume.shape)
            dist.integrate(ekin=np.s_[:], charge=np.s_[:], time=np.s_[:])
            #dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])

            r = np.transpose(moment.rc, (1,0,2))
            dist._distribution *= -(bphi/bnorm) *r *units

        moment.add_ordinates(collTorque=dist.histogram().to("J") / moment.volume)

    @staticmethod
    def canMomentTorque(dist, moment):
        """Calculate power deposition to ions.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        if not "ekin" in dist.abscissae or not "pitch" in dist.abscissae:
            raise ValueError("Distribution must be in energy-pitch basis.")
        if moment.rhodist:
            # Placeholder
            dist = dist.integrate(copy=True, charge=dist.abscissa("charge"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
        else:
            # Placeholder
            dist = dist.integrate(copy=True, charge=dist.abscissa("charge"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
        moment.add_ordinates(canMomentTorque=-dist.histogram() / moment.volume)

    @staticmethod
    def ppappe2ekinpitch(dist, mass, ekin_edges=10, pitch_edges=10):
        """Convert ppa and ppe distribution abscissae to energy and pitch.

        The conversion is done by first converting the input distribution to
        "markers" - one for each ppa-ppe cell and they have momentum equal
        to the midpoint of the cell. Then these "markers" are rebinned in
        the energy-pitch distribution. This transformation preserves the
        particle number unlike the alternative which would be to interpolate
        the particle density, but which otherwise could be more accurate.

        Energy is in electronvolts and pitch is ppa/(ppa^2 + ppe^2)^0.5. The
        transformation is relativistic.

        Parameters
        ----------
        dist : :class:`DistData`
            A ppar-pperp distribution.
        masskg : float
            Mass of the species (required for energy conversion).

            Note that distribution is assumed to consist of markers with equal
            mass for this conversion to work.
        ekin_edges : array_like or int, optional
            Energy grid edges or the number of bins in the new distribution.

            If edges are not given explicitly, the maximum energy is calculated
            from ppar or pperp abscissa (whichever gives the highest) and the
            minimum is set to zero.
        xi_edges : array_like or int, optional
            Pitch grid edges or the number of bins in the new distribution.

            If edges are not given explicitly, the pitch is binned in the
            [-1, 1] interval.

        Returns
        -------
        exidist : :class:`DistData`
            Energy-pitch distribution which is otherwise equivalent to the
            one given as input.
        """
        if isinstance(ekin_edges, int):
            p2max = np.maximum(dist.abscissa_edges("ppar")[-1]**2,
                               dist.abscissa_edges("ppar")[0]**2,
                               dist.abscissa_edges("pperp")[-1]**2)
            ekinmax = physlib.energy_momentum(mass, np.sqrt(p2max)).to("eV")
            ekin_edges = np.linspace(0*unyt.eV, ekinmax, ekin_edges)
        if isinstance(pitch_edges, int):
            pitch_edges = np.linspace(-1, 1, pitch_edges)

        try:
            ekin_edges.units
        except AttributeError:
            ekin_edges *= unyt.eV
        try:
            pitch_edges.units
        except AttributeError:
            pitch_edges *= unyt.dimensionless

        # Create a new empty distribution where ppar and pperp are replaced
        # by ekin and pitch
        dim = []
        abscissa_edges = {}
        for k in dist.abscissae:
            if k == "ppar":
                dim.append(ekin_edges.size-1)
                abscissa_edges["ekin"] = ekin_edges
            elif k == "pperp":
                dim.append(pitch_edges.size-1)
                abscissa_edges["pitch"] = pitch_edges
            else:
                dim.append(dist.abscissa_edges(k).size-1)
                abscissa_edges[k] = dist.abscissa_edges(k)
        exdist = DistData(np.zeros(dim)*unyt.particles, **abscissa_edges)

        # Transform E-xi grid to points in (ppa,ppa) space that are used in
        # interpolation.
        xi, ekin = np.meshgrid(
            exdist.abscissa("pitch"), exdist.abscissa("ekin") )
        pnorm = physlib.pnorm_gamma(
            mass, physlib.gamma_energy(mass, ekin.ravel()) )
        ppa = ( xi.ravel() * pnorm ).to("amu*m/s").v
        ppe = ( np.sqrt(1 - xi.ravel()**2 ) * pnorm).to("amu*m/s").v

        ## NOTE: Jacobian and some other quantities would be needed only for
        ## the interpolation method, which was replaced by this histogram
        ## approach. They are kept here in case there is need to use the
        ## interpolation method in some cases.

        # Coordinate transform Jacobian: dppa dppe = |jac| dE dxi
        # Jacobian for transform (ppa, ppe) -> (p, xi) is p / sqrt(1-xi^2)
        # because jac = dppa / dp  = xi, dppe / dp  = sqrt(1-xi^2)
        #               dppa / dxi = p,  dppe / dxi = -xi p / sqrt(1-xi^2),
        # and the Jacobian for (p, xi) -> (E, xi) is m gamma / p.
        #
        # Therefore the combined Jacobian is:
        # ( m gamma / p ) / sqrt(1-xi*xi).
        jac = (mass + ekin / unyt.c**2) / np.sqrt(1 - xi**2)

        # Quantities needed in iteration on each loop
        ippa     = dist.abscissae.index("ppar")
        units    = dist.distribution().units
        exishape = (exdist.abscissa("ekin").size, exdist.abscissa("pitch").size)
        ppar     = dist.abscissa("ppar").to("amu*m/s").v
        pperp    = dist.abscissa("pperp").to("amu*m/s").v

        ppa0,ppe0 = np.meshgrid(dist.abscissa("ppar"), dist.abscissa("pperp"), indexing="ij")
        ppa0 = ppa0.ravel()
        ppe0 = ppe0.ravel()
        pnorm0 = np.sqrt(ppa0**2 + ppe0**2)
        pitch0 = ppa0 / pnorm0
        ekin0 = ((physlib.gamma_momentum(mass, pnorm0) - 1) \
                 * mass * unyt.c**2).to("eV")

        # Use itertools to conveniently make N "for" loops into a single loop
        ranges = []
        for a in dist.abscissae:
            if a != "ppar" and a != "pperp":
                ranges.append(range(dist.abscissa(a).size))

        #
        ie = np.digitize(ekin0,  ekin_edges)-1
        ip = np.digitize(pitch0, pitch_edges)-1
        mask = np.logical_or.reduce([ie < 0, ip < 0, ie >= ekin_edges.size-1, ip >= pitch_edges.size-1])
        ie = ie[~mask]
        ip = ip[~mask]
        vol = exdist.phasespacevolume()
        hist = dist.histogram()

        # We loop over all dimensions except ppar and pperp
        for itr in itertools.product(*ranges):

            # Consctruct indices to slice (ppa, ppa) at this iteration
            idx = [slice(None)] * (len(itr) + 2)
            idx[:ippa]   = itr[:ippa] # Coordinates before ppa
            idx[ippa+2:] = itr[ippa:] # Coordinates after ppa and ppe
            idx = tuple(idx)

            # Interpolate distribution at (ekin, xi) grid points
            #f = RectBivariateSpline(
            #    ppar, pperp, np.squeeze(dist.distribution()[idx]), kx=1, ky=1)
            #d = np.reshape(f.ev(ppa, ppe), exishape) * units

            # Apply jacobian and store the values to exidist
            #exdist._distribution[idx] = d * jac

            #d = np.histogram2d(
            #    ekin0, pitch0, bins=(ekin_edges, pitch_edges),
            #    weights=np.squeeze(dist.histogram()[idx]).T.ravel())[0]
            #exdist._distribution[idx] = d / vol.units
            a = np.zeros(exdist._distribution[idx].shape)
            np.add.at(a, (ie,ip), hist[idx].v.ravel()[~mask])
            exdist._distribution[idx] = a / vol.units

        exdist._distribution /= vol.v

        return exdist

def plot(self, axes=None, cax=None, logscale=False):
    """Plot distribution in 1D or 2D.

    This method assumes that the input distribution has been integrated,
    sliced, and interpolated so that only one or two dimensions have
    a size above one.

    Parameters
    ----------
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    logscale: bool, optional
        Whether the plot is in logarithmic scale.
    """
    x, y = None, None
    for key in self.abscissae:
        val = self.abscissa_edges(key)
        if val.size > 2:
            if x is None:
                x = val
                xlabel = key + " [" + str(x.units) + "]"
            elif y is None:
                y = val
                ylabel = key + " [" + str(y.units) + "]"
            else:
                raise ValueError(
                    "The distribution has more than two dimensions with "
                    + "size greater than one")
    if x is None: raise ValueError("The distribution is zero dimensional")

    ordinate = np.squeeze(self.distribution())
    if y is None:
        ylabel = "f" + " [" + str(ordinate.units) + "]"
        a5plt.mesh1d(x, ordinate, xlabel=xlabel, ylabel=ylabel, axes=axes,
                        logscale=logscale)
    else:
        axesequal = x.units == y.units
        clabel = "f" + " [" + str(ordinate.units) + "]"
        a5plt.mesh2d(x, y, ordinate, axesequal=axesequal, xlabel=xlabel,
                        ylabel=ylabel, clabel=clabel, axes=axes, cax=cax,
                        logscale=logscale)
