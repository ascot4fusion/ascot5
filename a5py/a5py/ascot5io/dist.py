"""Distribution output.
"""
import numpy as np
import unyt

import itertools
from scipy.interpolate import griddata, RectBivariateSpline
from a5py import physlib

from .coreio.treedata import DataContainer

class DistData():
    """Distribution data object.

    Attributes
    ----------
    abscissae : list [str]
        Name of each abscissa in same order as they appear in the distribution.
    _distribution : array_like
        N-dimensional array, where N is the number of abscissae, containing
        the distribution data.
    **abscissa_edges : array_like
        Edges of each abscissa stored as attributes with name
        "_<name of the abscissa>" e.g. "_r".
    """

    def __init__(self, histogram, **abscissa_edges):
        """Initialize distribution data from given abscissae and histogram.

        Parameters
        ----------
        histogram : array_like
            N-dimensional array, where N is the number of abscissae, containing
            the distribution data as a histogram.
        **abscissa_edges : array_like
            Name and edges of each abscissa.
        """
        self.abscissae = list(abscissa_edges.keys())
        for k in abscissa_edges:
            setattr(self, "_" + k, abscissa_edges[k])

        self._distribution = histogram / self.phasespacevolume()

    def _copy(self):
        """Make a (deep) copy of this object.

        Returns
        -------
        copy : :class:`DistData`
            Copy of this object.
        """
        abscissa_edges = {}
        for name in self.abscissae:
            abscissa_edges[name] = getattr(self, "_" + name).copy()
        histogram = self.histogram().copy()

        return DistData(histogram, **abscissa_edges)

    def abscissa(self, name):
        """Return abscissa values.

        Parameters
        ----------
        name : str
            Name of the abscissa.

        Returns
        -------
        abscissa : array_like
            Abscissa values, i.e., the values at the center of each bin.

        Raises
        ------
        ValueError
            If ``name`` does not correspond to any abscissa.
        """
        abscissa = "_" + name
        if not hasattr(self, abscissa):
            raise ValueError("Unknown abscissa: " + name)
        edges = getattr(self, abscissa)
        return ( edges[1:] + edges[:-1] ) / 2

    def abscissa_edges(self, name):
        """Return abscissa edges.

        Parameters
        ----------
        name : str
            Name of the abscissa.

        Returns
        -------
        abscissa : array_like
            Abscissa bin edges.

        Raises
        ------
        ValueError
            If ``name`` does not correspond to any abscissa.
        """
        abscissa = "_" + name
        if not hasattr(self, abscissa):
            raise ValueError("Unknown abscissa: " + name)
        return getattr(self, abscissa).copy()

    def distribution(self):
        """Return the distribution function.

        Returns
        -------
        dist : array_like
            N-dimensional array, where N is the number of abscissae, containing
            the distribution data.
        """
        return self._distribution

    def histogram(self):
        """Return the distribution as a histogram (particles per bin).

        Returns
        -------
        dist : array_like
            N-dimensional array, where N is the number of abscissae, containing
            the distribution data as a histogram.
        """
        return self._distribution * self.phasespacevolume()

    def slice(self, copy=False, **abscissae):
        """Take slices of distribution.

        Parameters
        ----------
        copy : bool, optional
            Retain original distribution and return a copy which is sliced.
        **abscissae : slice
            Name of the coordinate and corresponding slice to be taken.

        Returns
        -------
        dist : :class:`DistData`
            The sliced distribution if ``copy=True``.
        """
        dist = self if not copy else self._copy()

        for k, s in abscissae.items():
            dim = dist.abscissae.index(k)
            idx = np.histogram(dist.abscissa(k)[s], dist.abscissa_edges(k))[0]
            idx = np.nonzero(idx)[0]

            endpoints = dist.abscissa_edges(k)[[0,-1]]
            midpoints = dist.abscissa(k)[s]
            edges = np.zeros((midpoints.size+1,))
            edges[0]  = endpoints[0]
            edges[-1] = endpoints[-1]
            if edges.size > 2:
                edges[1:-1] = ( midpoints[1:] + midpoints[:-1] ) / 2
            setattr(dist, "_" + k, edges)

            dist._distribution = dist._distribution.take(indices=idx, axis=dim)

        if copy: return dist

    def integrate(self, copy=False, **abscissae):
        """Integrate distribution along the given dimension.

        Parameters
        ----------
        copy : bool, optional
            Retain original distribution and return a copy which is integrated.
        **abscissae : slice
            Name of the coordinate and corresponding slice which is integrated.

        Returns
        -------
        dist : :class:`DistData`
            The integrated distribution if ``copy=True``.
        """
        for k in abscissae.keys():
            if not hasattr(self, "_" + k):
                raise ValueError("Unknown abscissa: " + k)

        dist = self if not copy else self._copy()
        dim = 0
        for k in tuple(dist.abscissae):
            if k not in abscissae.keys():
                dim +=1
                continue

            idx = abscissae[k]
            if isinstance(idx, slice):
                ds   = np.diff(dist.abscissa_edges(k))
                mask = np.zeros(ds.shape)
                mask[idx] = 1
                ds  *= mask
            else:
                ds = np.diff(dist.abscissa_edges(k)) * idx

            dist._distribution, s = np.average(
                dist.distribution(), axis=dim, weights=ds, returned=True)

            # np.average seems to strip units from ds most of the time but not
            # always...
            if hasattr(s, "units"):
                dist._distribution *= s
            else:
                dist._distribution *= s * ds.units

            dist.abscissae.remove(k)
            delattr(dist, "_" + k)

        if copy: return dist

    def interpolate(self, **coordinates):
        """Perform (N-dim) linear interpolation on the distribution.

        Parameters
        ----------
        **coordinates : slice
            Name of the coordinate and value at which the distribution is
            interpolated.

        Returns
        -------
        dist : array_like
            The distribution interpolated at given point.
        """
        for k in coordinates.keys():
            if not hasattr(self, "_" + k):
                raise ValueError("Unknown abscissa: " + k)
            val = coordinates[k]
            edges = self.abscissa_edges(k)
            dim = edges.units.dimensions
            try:
                # Try to get argument units
                valdim = val.units.dimensions
            except AttributeError:
                # Argument doesn't have units, assign and add warning
                val = val*unit
                valdim = dim
                warnings.warn(AscotUnitWarning)
            if valdim != dim:
                raise ValueError(
                    "\"%s\" has incorrect dimensions: expected %s but got %s" %
                    (name, dim, valdim))

            if val < edges[0] or val > edges[-1]:
                raise ValueError(
                    "Coordinate out of bounds [%.2f, %.2f]: %s = %.2f",
                    (edges[0], edges[-1], k, val) )

        dim = 0
        dist = self.distribution().copy()
        for k in tuple(self.abscissae):
            if k not in coordinates.keys():
                dim += 1
                continue

            val = coordinates[k]
            edges = self.abscissa_edges(k)
            i = np.argmax(val < edges) - 1
            t = ( val - edges[i] ) / ( edges[i+1] - edges[i] )

            dist0 = dist.take(indices=i,   axis=dim)
            dist1 = dist.take(indices=i+1, axis=dim)

            dist = t * dist0 + (1-t) * dist1
        return dist

    def phasespacevolume(self):
        """Calculate phase-space volume of each bin.

        Returns
        -------
        vol : array_like
            Phase space volumes with units and same shape as the distribution.
        """
        vol = 1
        for name in self.abscissae:
            edges = self.abscissa_edges(name)
            dV    = (edges[1:] - edges[:-1])
            vol   = np.multiply.outer(vol, dV)
        return vol

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

class Dist(DataContainer):

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
        if moment.rhodist:
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
            dist = dist.integrate(copy=True, **integrate)
        else:
            integrate = {}
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
        if moment.rhodist:
            dist = dist.integrate(copy=True, charge=dist.abscissa("charge"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
        else:
            dist = dist.integrate(copy=True, charge=dist.abscissa("charge"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
        moment.add_ordinates(chargedensity=dist.histogram() / moment.volume)

    @staticmethod
    def energydensity(dist, moment):
        """Calculate energy density.

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
            dist = dist.integrate(copy=True, ekin=dist.abscissa("ekin"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]  
            dist.integrate(**integrate)
        else:
            dist = dist.integrate(copy=True, ekin=dist.abscissa("ekin"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
        moment.add_ordinates(energydensity=dist.histogram().to("J") / moment.volume)

    @staticmethod
    def pressure(ascot, mass, dist, moment):
        """Calculate pressure.

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
            vnorm2 = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) ) **2
            dist = dist.integrate(copy=True, pitch=np.s_[:], ekin=vnorm2, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= mass / 3.0
            
        else:
            vnorm2 = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) ) **2
            dist = dist.integrate(copy=True, pitch=np.s_[:], ekin=vnorm2, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= mass / 3.0
            
        moment.add_ordinates(pressure=dist.histogram().to("J") / moment.volume)
        
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
        if not "ppar" in dist.abscissae or not "pperp" in dist.abscissae:
            raise ValueError("Distribution must be in ppar-pperp basis.")
        if moment.rhodist:
            dist = dist.integrate(copy=True,
                                  charge=dist.abscissa("charge"),
                                  ppar=dist.abscissa("ppar"))
            bphi, bnorm = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm")
            bphi  = bphi.reshape(moment.volume.shape)
            bnorm = bnorm.reshape(moment.volume.shape)
            
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
            dist._distribution *= bphi / (bnorm * mass)
        else:
            dist = dist.integrate(copy=True, 
                                  charge=dist.abscissa("charge") ,
                                  ppar=dist.abscissa("ppar"))
            bphi, bnorm = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm")
            bphi  = bphi.reshape(moment.volume.shape)
            bnorm = bnorm.reshape(moment.volume.shape)

            integrate = {}
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
            dist._distribution *= bphi / (bnorm * mass)
        moment.add_ordinates(toroidalcurrent=(dist.histogram() / moment.volume).to("A/m**2"))

    @staticmethod
    def parallelcurrent(ascot, mass, dist, moment):
        """Calculate pressure.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        if not "ppar" in dist.abscissae or not "pperp" in dist.abscissae:
            raise ValueError("Distribution must be in ppar-pperp basis.")
        if moment.rhodist:
            dist = dist.integrate(copy=True,
                                  charge=dist.abscissa("charge"),
                                  ppar=dist.abscissa("ppar"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
            dist._distribution *= 1/mass
        else:
            dist = dist.integrate(copy=True, 
                                  charge=dist.abscissa("charge") ,
                                  ppar=dist.abscissa("ppar"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
            dist._distribution *= 1/mass
            
        moment.add_ordinates(parallelcurrent=(dist.histogram() / moment.volume).to("A/m**2"))

    @staticmethod
    def powerdep(ascot, mass, dist, moment):
        """Calculate number density.

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
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = np.sum(coefs["K"],axis=0).reshape((*moment.volume.shape,va.size))
                units = K.units
                dist._distribution[:,:,:,:,0,0] *= K.v
            dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= units * mass
        else:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = np.sum(coefs["K"],axis=0).reshape((*moment.volume.shape,va.size))
                units = K.units
                dist._distribution[:,:,:,:,0,0] *= K.v
            dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= units * mass
        moment.add_ordinates(powerdep=(-dist.histogram() / moment.volume ).to("W/m**3"))

        
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
        if not "ekin" in dist.abscissae or not "pitch" in dist.abscissae:
            raise ValueError("Distribution must be in energy-pitch basis.")
        if moment.rhodist:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = coefs["K"][0].reshape((*moment.volume.shape,va.size))
                units = K.units
                dist._distribution[:,:,:,:,0,0] *= K.v
            dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= units * mass
        else:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = coefs["K"][0].reshape((*moment.volume.shape,va.size))
                units = K.units
                dist._distribution[:,:,:,:,0,0] *= K.v
            dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= units * mass
        moment.add_ordinates(electronpowerdep=(-dist.histogram() / moment.volume).to("W/m**3"))
        
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
        if not "ekin" in dist.abscissae or not "pitch" in dist.abscissae:
            raise ValueError("Distribution must be in energy-pitch basis.")
        if moment.rhodist:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = np.sum(coefs["K"][1:],axis=0).reshape((*moment.volume.shape,va.size))
                units = K.units
                dist._distribution[:,:,:,:,0,0] *= K.v
            dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= units * mass
        else:
            va = physlib.vnorm_gamma(
                physlib.gamma_energy( mass, dist.abscissa("ekin") ) )
            dist = dist.integrate(copy=True, pitch=np.s_[:])
            for qa in dist.abscissa("charge"):
                coefs = ascot.input_eval_collcoefs(
                    mass.to("kg"), qa.to("C"), moment.rc.ravel(), moment.phic.ravel(), moment.zc.ravel(),
                    np.zeros(moment.rc.ravel().shape)*unyt.s, va)
                K = np.sum(coefs["K"][1:],axis=0).reshape((*moment.volume.shape,va.size))
                units = K.units
                dist._distribution[:,:,:,:,0,0] *= K.v
            dist.integrate(ekin=va, charge=np.s_[:], time=np.s_[:])
            dist._distribution *= units * mass
        moment.add_ordinates(ionpowerdep=(-dist.histogram() / moment.volume).to("W/m**3"))


    @staticmethod
    def jxBTorque(ascot, mass, dist, moment):
        """Calculate power deposition to ions.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
        if not "ppar" in dist.abscissae or not "pperp" in dist.abscissae:
            raise ValueError("Distribution must be in ppar-pperp basis.")
        if moment.rhodist:
            bphi, bnorm, br, bz = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm", "br", "bz")
            dbrdr, dbphidr , dbzdr  = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "brdr", "bphidr", "bzdr")
            dbrdphi, dbphidphi , dbzdphi  = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "brdphi", "bphidphi", "bzdphi")
            dbrdz, dbphidz , dbzdz  = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "brdz", "bphidz", "bzdz")
            
            GradBr   = (br*dbrdr + bphi*dbphidr + bz*dbzdr) /bnorm
            GradBphi = (br*dbrdphi + bphi*dbphidphi + bz*dbzdphi) /bnorm
            GradBz   = (br*dbrdz + bphi*dbphidz + bz*dbzdz) /bnorm

            r      = moment.rc
            b = np.array([br/bnorm, bphi/bnorm, bz/bnorm]).T
            GradB = np.array([GradBr, GradBphi, GradBz]).T
            # curlB brings an issue: #
            # "The requested array has an inhomogeneous shape after 1 dimensions."
            curlB = np.array([dbzdphi /r - dbphidz, dbrdz - dbzdr, (bphi-dbrdphi)/r -dbphidr]).T
            GradBcrossB = np.cross(GradB ,b)
            crossb = (curlB - GradBcrossB) / bnorm

            for qa in dist.abscissa("charge"):
                dist._distribution[:,:,:,:,0,0] *= bnorm/qa + dist.abscissa("ppar")*crossb
            
            bpol = np.sqrt(br**2 + bz**2) .reshape(moment.volume.shape)
            r      = moment.rc .reshape(moment.volume.shape)
            q = dist.abscissa("charge")
            
            dist = dist.integrate(copy=True, ppar=dist.abscissa("ppar"))
            integrate = {}
            for k in dist.abscissae:
                if k not in ["rho", "theta", "phi"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
                    
            #dist = dist.integrate(time=np.s_[:])
            dist._distribution *= 1/mass * bpol * r
            
            
        else:
            # Placeholder
            dist = dist.integrate(copy=True, 
                                  charge=dist.abscissa("charge") ,
                                  ppar=dist.abscissa("ppar"))
            bphi, bnorm, br, bz = ascot.input_eval(
                moment.rc, moment.phic, moment.zc, 0*unyt.s, "bphi", "bnorm", "br", "bz")
            btheta = np.sqrt(br**2 + bz**2)
            
            bphi   = bphi.reshape(moment.volume.shape)
            bnorm  = bnorm.reshape(moment.volume.shape)
            btheta = btheta.reshape(moment.volume.shape)
            r      = moment.rc.reshape(moment.volume.shape)
            
            integrate = {}
            for k in dist.abscissae:
                if k not in ["r", "phi", "z"]:
                    integrate[k] = np.s_[:]
            dist.integrate(**integrate)
            dist._distribution *= -(bphi / (bnorm * mass)) * btheta *r
            
            
        moment.add_ordinates(jxBTorque=dist.histogram() / moment.volume)

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

        This function operates by looping through all other coordinates except
        ppa and ppe, and at each loop calculates
        f_Exi(E, xi) = f_ppappe(ppa(E_i, xi_i), ppe(E_i, xi_i)) where E_i and
        xi_i are grid points of the new energy-pitch distribution. Interpolation
        is done bilinearly.

        Energy is in electronvolts and pitch is ppa/(ppa^2 + ppe^2)^0.5. The
        transformation is not relativistic.

        dist : dict_like <br>
            A ppa-ppe distribution. May hold other dimensions as well.
        masskg : float <br>
            Mass of the species (required for energy conversion) in kg. Note
            that distribution is assumed to consist of markers with equal mass.
        E_edges : array_like or int, optional <br>
            Energy grid edges in the new distribution. If not given,
            linspace(0, Emax, 10) will be used where Emax is
            e*0.5*masskg*max(vpa^2, vpe^2). If an integer is given then
            linspace(0, Emax, E_edges) is used.
        xi_edges : array_like or int, optional <br>
            Pitch grid edges in the new distribution. If not given,
            linspace(-1, 1, 10) will be used. If an integer is given then
            linspace(-1, 1, xi_edges) is used.

        Returns:
        Energy-pitch distribution dictionary whose other dimensions are same as
        in input.
        """
        if isinstance(ekin_edges, int):
            p2max = np.maximum(dist.abscissa_edges("ppar")[-1]**2,
                               dist.abscissa_edges("ppar")[0]**2,
                               dist.abscissa_edges("pperp")[-1]**2)
            ekinmax = physlib.energy_momentum(mass, np.sqrt(p2max)).to("eV")
            ekin_edges = np.linspace(0, ekinmax, ekin_edges)
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

        # Create new, empty, distribution where ppar and pperp are replaced
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

        # Transform E-xi grid to points in (vpa,vpa) space that are used in
        # interpolation.
        xi, ekin = np.meshgrid(exdist.abscissa("pitch"), exdist.abscissa("ekin"))
        pnorm = physlib.pnorm_gamma(mass, physlib.gamma_energy(mass, ekin.ravel()))
        ppa = (xi.ravel() * pnorm).to("amu*m/s").v
        ppe = (np.sqrt(1 - xi.ravel()**2) * pnorm).to("amu*m/s").v

        #p2pa, p2pe = np.meshgrid(exdist.abscissa("ppar"), exdist.abscissa("pperp"))
        #pnorm   = np.sqrt(p2pa**2 + p2pe**2)
        #p2pitch = p2pa / pnorm
        #p2ekin  = (np.sqrt(1 + (pnorm/(mass*unyt.c))) - 1)*

        # Coordinate transform Jacobian: dppa dppe = |jac| dE dxi
        # Jacobian for transform (ppa, ppe) -> (p, xi) is p / sqrt(1-xi^2)
        # because jac = dppa / dp  = xi, dppe / dp  = sqrt(1-xi^2)
        #               dppa / dxi = p,  dppe / dxi = -xi p / sqrt(1-xi^2),
        # and the Jacobian for (p, xi) -> (E, xi) is e E0 / c^2 p when
        # E is in electronvolts. Therefore the combined Jacobian is
        # (e E0 / c^2) / sqrt(1-xi*xi).
        #E0  = np.sqrt( (pg*constants.c)**2 + (masskg.v*constants.c**2)**2 )
        #jac = (constants.e * E0 / constants.c**2) / np.sqrt(1 - xig*xig)
        jac = (mass + ekin / unyt.c**2) / np.sqrt(1 - xi**2)

        # Quantities needed in iteration on each loop
        ippa = dist.abscissae.index("ppar")
        units = dist.distribution().units
        exishape = (exdist.abscissa("ekin").size, exdist.abscissa("pitch").size)
        ppar  = dist.abscissa("ppar").to("amu*m/s").v
        pperp = dist.abscissa("pperp").to("amu*m/s").v

        # Use itertools to conveniently make N "for" loops into a single loop
        ranges = []
        for a in dist.abscissae:
            if a != "ppar" and a != "pperp":
                ranges.append(range(dist.abscissa(a).size))

        # We loop over all dimensions except ppar and pperp
        for itr in itertools.product(*ranges):

            # Consctruct indices to slice (ppa, ppa) at this iteration
            idx = [slice(None)] * (len(itr) + 2)
            idx[:ippa]   = itr[:ippa] # Coordinates before ppa
            idx[ippa+2:] = itr[ippa:] # Coordinates after ppa and ppe
            idx = tuple(idx)

            # Interpolate distribution at (ekin, xi) grid points
            f = RectBivariateSpline(
                ppar, pperp, np.squeeze(dist.distribution()[idx]), kx=1, ky=1)
            d = np.reshape(f.ev(ppa, ppe), exishape) * units

            # Apply jacobian and store the values to exidist
            exdist._distribution[idx] = d * jac

        return exdist

class Dist_5D(Dist):

    def read(self):
        """Read 5D distribution from a HDF5 file to a dictionary.
        """
        out = {}
        with self as f:

            # A Short helper function to calculate grid points from grid edges.
            def edges2grid(edges):
                return np.linspace(0.5*(edges[0]+edges[1]),
                                   0.5*(edges[-2]+edges[-1]), num=edges.size-1)

            abscissae = [0] * int(f["abscissa_ndim"][:])
            for i in range(0, len(abscissae)):
                abscissa     = f["abscissa_vec_0"+str(i+1)]
                name         = abscissa.attrs["name_0"+str(i)].decode("utf-8")
                abscissae[i] = name

                unit = abscissa.attrs["unit_0"+str(i)].decode("utf-8")
                out[name + "_edges"] = abscissa[:] * unyt.Unit(unit)
                out[name]            = edges2grid(out[name + "_edges"])
                out["n" + name]      = out[name].size

            out["abscissae"] = abscissae
            out["histogram"] = f["ordinate"][0,:,:,:,:,:,:,:] \
                * unyt.dimensionless

        return out

class Dist_6D(Dist):
    pass

class Dist_rho5D(Dist):
    pass

class Dist_rho6D(Dist):
    pass

class Dist_COM(Dist):

    def read(self):
        """Read 5D distribution from a HDF5 file to a dictionary.
        """
        out = {}
        with self as f:

            # A Short helper function to calculate grid points from grid edges.
            def edges2grid(edges):
                return np.linspace(0.5*(edges[0]+edges[1]),
                                   0.5*(edges[-2]+edges[-1]), num=edges.size-1)

            abscissae = [0] * int(f["abscissa_ndim"][:])
            for i in range(0, len(abscissae)):
                abscissa     = f["abscissa_vec_0"+str(i+1)]
                name         = abscissa.attrs["name_0"+str(i)].decode("utf-8")
                abscissae[i] = name

                unit = abscissa.attrs["unit_0"+str(i)].decode("utf-8")
                out[name + "_edges"] = abscissa[:] #* unyt.Unit(unit)
                out[name]            = edges2grid(out[name + "_edges"])
                out["n" + name]      = out[name].size

            out["abscissae"] = abscissae
            out["histogram"] = f["ordinate"][0,:,:,:] \
                * unyt.dimensionless

        return out

    
    def plot_muptor(self, E, axes=None):
        """
        Plot constant of motion distribution on (mu, Ptor) plane for a given E.
        """
        import matplotlib.pyplot as plt
        data = self.read()
        hist = data["histogram"]
        mu_edges, Ptor_edges = data["mu_edges"], data["ptor_edges"]
        E_vector, E_edges = data["ekin"], data["ekin_edges"]
        if E is None:
            plot = np.sum(hist, axis=1)
        elif type(E) == int:
            i = E
            plot = hist[:,i,:]
        else:
            i = np.argmax(E < E_edges)-1
            plot1 = hist[:,i,:]
            plot2 = hist[:,i+1,:]

            k = E/(E_vector[i+1]-E_vector[i])
            plot = plot1*k + (1-k)*plot2

        h = axes.pcolormesh(Ptor_edges, mu_edges, plot, shading="flat")
        axes.set_ylabel("mu (J/T)")
        axes.set_xlabel("Ptor (Js)")
        plt.colorbar(h, ax=axes)

    def plot_muEkin(self, Ptor, axes=None):
        """
        Plot constant of motion distribution on (mu, Ptor) plane for a given E.
        """
        import matplotlib.pyplot as plt
        data = self.read()
        hist = data["histogram"]
        mu_edges, E_edges = data["mu_edges"], data["ekin_edges"]
        Ptor_vector, Ptor_edges = data["ptor"], data["ptor_edges"]
        if Ptor is None:
            plot = np.sum(hist, axis=2)
        elif type(Ptor) == int:
            i = Ptor
            plot = hist[:,:,i]
        else:
            i = np.argmax(Ptor < Ptor_edges)-1
            plot1 = hist[:,:,i]
            plot2 = hist[:,:,i+1]

            k = Ptor/(Ptor_vector[i+1]-Ptor_vector[i])
            plot = plot1*k + (1-k)*plot2
        
        h = axes.pcolormesh(E_edges, mu_edges, plot, shading="flat")
        axes.set_ylabel("mu (J/T)")
        axes.set_xlabel("E (J)")
        plt.colorbar(h, ax=axes)

        
    def plot_EkinPtor(self, mu, axes=None):
        """
        Plot constant of motion distribution on (mu, Ptor) plane for a given E.
        """
        import matplotlib.pyplot as plt
        data = self.read()
        hist = data["histogram"]
        E_edges, Ptor_edges = data["ekin_edges"], data["ptor_edges"]
        mu_vector, mu_edges = data["mu"], data["mu_edges"]
        if mu is None:
            plot = np.sum(hist, axis=0)
        elif type(mu) == int:
            i = mu
            plot = hist[i,:,:]
        else:
            i = np.argmax(mu < mu_edges)-1
            plot1 = hist[i,:,:]
            plot2 = hist[i+1,:,:]

            k = mu/(mu_vector[i+1]-mu_vector[i])
            plot = plot1*k + (1-k)*plot2
        
        h = axes.pcolormesh(Ptor_edges, E_edges, plot, shading="flat")
        axes.set_ylabel("E (J)")
        axes.set_xlabel("Ptor (Js)")
        plt.colorbar(h, ax=axes)
