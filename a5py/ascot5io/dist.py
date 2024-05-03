"""Distribution output.
"""
import numpy as np
import unyt

import itertools
from scipy.interpolate import griddata, RectBivariateSpline
from a5py import physlib
import a5py.routines.plotting as a5plt

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

    def _multiply(self, fac, *abscissae):
        """Multiply distribution with a value.

        This method is used to calculate moments from the distribution.
        """
        ranges = []; indices = []
        for i, a in enumerate(self.abscissae):
            if a not in abscissae:
                ranges.append(range(self.abscissa(a).size))
                indices.append(i)

        idx = [slice(None)] * len(self.abscissae)
        for itr in itertools.product(*ranges):
            for k, i in enumerate(indices):
                idx[i] = itr[k]
            self._distribution[tuple(idx)] *= fac.v

        self._distribution *= fac.units

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

            if isinstance(s, int):
                if s < 0:
                    N_len = len(dist.abscissa(k))
                    ns = np.s_[N_len+s:N_len+s+2]
                else:
                    ns = np.s_[s:s+2]
            elif isinstance(s, slice):
                N_len = len(dist.abscissa(k))
                if s.stop is None:
                    stop = N_len+1
                elif s.stop < 0:
                    stop = N_len+s.stop+1
                else:
                    stop = s.stop+1
                if s.start is None:
                    start = 0
                elif s.start < 0:
                    start = N_len+s.start
                else:
                    start = s.start
                ns = np.s_[start:stop]

            edges = dist.abscissa_edges(k)[ns]

            setattr(dist, "_" + k, edges)

            dist._distribution = dist._distribution.take(indices=idx, axis=dim)

        if copy: return dist

    def integrate(self, copy=False, **abscissae):
        """Integrate distribution along the given dimension.

        Parameters
        ----------
        copy : bool, optional
            Retain original distribution and return a copy which is integrated.
        **abscissae : slice or array_like
            Name of the coordinate and corresponding slice which is integrated.

            If argument is an array, it must have the same size as the
            corresponding dimension. The integration is then performed as
            int f(x,...) * array(x) * dx.

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
        x = None; y = None;
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
            # Swap pitch and energy
            if "ekin" in xlabel and "pitch" in ylabel:
                temp = x; x = y; y = temp
                temp = xlabel; xlabel = ylabel; ylabel = temp
                ordinate = ordinate.T

            axesequal = x.units == y.units
            clabel = "f" + " [" + str(ordinate.units) + "]"
            a5plt.mesh2d(x, y, ordinate, axesequal=axesequal, xlabel=xlabel,
                         ylabel=ylabel, clabel=clabel, axes=axes, cax=cax,
                         logscale=logscale)

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
        dist._multiply(ekin.reshape(ppa.shape), "ppar", "pperp")
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
        vnorm = physlib.velocity_momentum(mass, pnorm)
        dist._multiply(vnorm.reshape(ppa.shape)**2 * mass / 3, "ppar", "pperp")
        dist.integrate(ppar=np.s_[:], pperp=np.s_[:])
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
    def parallelcurrent(ascot, mass, dist, moment):
        """Calculate parallel current density.

        Parameters
        ----------
        dist : :class:`DistData`
            Distribution from where the moments are calculated.
        moment : class:`DistMoment`
            Moment data used in calculation and where the result is stored.
        """
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
        """Read COM distribution from a HDF5 file to a dictionary.
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
            out["histogram"] = f["ordinate"][0,:,:,:] \
                * unyt.dimensionless

        return out
