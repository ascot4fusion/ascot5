"""Histogram diagnostics.
"""
import ctypes
from typing import Optional, Sequence
from types import MappingProxyType

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

import itertools
from scipy.interpolate import RegularGridInterpolator


from a5py import utils
from a5py.physlib import Species
from a5py.libascot import LIBASCOT, DataStruct
from a5py.exceptions import AscotMeltdownError
from a5py.data.options import HistParams, DIMENSIONS



# pylint: disable=too-few-public-methods
class HistAxis(ctypes.Structure):
    _fields_ = [
        ("min", ctypes.c_double),
        ("max", ctypes.c_double),
        ("n", ctypes.c_size_t),
    ]

# pylint: disable=too-few-public-methods
class Struct(DataStruct):

    _fields_ = [
        ("nbin", ctypes.c_size_t),
        ("strides", ctypes.c_size_t * (len(DIMENSIONS) - 1)),
        ("dimensions", HistAxis * len(DIMENSIONS)),
        ("values", ctypes.POINTER(ctypes.c_double)),
    ]


class Dist():
    """N-dimensional phase-space distribution on an equally spaced grid.

    The distribution is internally stored as a phase-space density,
    i.e. histogram values divided by the corresponding phase-space
    bin volume.

    Parameters
    ----------
    histogram : np.ndarray
        N-dimensional histogram values or particle counts per bin.
    bin_edges : dict[str, np.ndarray]
        Mapping from dimension name to bin edge coordinates.
        Each edge array must have length ``n_bins + 1`` for the
        corresponding histogram axis.
    species : str, optional
        Name of the particle species represented by the distribution.
    charge : int, optional
        Charge state of the particle species. Must be specified if
        ``species`` is given.

    Raises
    ------
    ValueError
        If ``species`` is specified but ``charge`` is not.
    """

    def __init__(
            self, histogram: np.ndarray,
            bin_edges: dict[str, np.ndarray],
            species: Optional[str | Species] = None,
            charge: Optional[int] = None,
            ):
        """Initialize distribution from a histogram data."""
        self._bin_edges = {}
        if len(bin_edges) != histogram.ndim:
            raise ValueError(
                "Bin edges must have the same number of dimensions as the "
                "histogram."
            )
        for (key, edges), exp_size in zip(bin_edges.items(), histogram.shape):

            if edges.ndim != 1:
                raise ValueError(
                    f"Bin edges for dimension '{key}' must be 1-dimensional."
                )

            if edges.size < 2:
                raise ValueError(
                    f"Bin edges for dimension '{key}' must contain at least "
                    "two entries."
                )

            dv = edges[1:] - edges[:-1]
            if not np.allclose(dv, dv[0]):
                raise ValueError(
                    f"Dimension '{key}' is not uniformly spaced."
                )

            if exp_size != edges.size - 1:
                raise ValueError(
                    f"Bin edges for dimension '{key}' must have length "
                    f"{exp_size + 1}."
                )

            self._bin_edges[key] = unyt.unyt_array(edges)

        self._distribution = histogram / self.phasespacevolume
        self._species, self._charge = None, None
        if species is not None:
            self._species = (
                species if isinstance(species, Species)
                else Species.from_string(species)
            )
            self._charge = charge
        if species is not None and charge is None:
            raise ValueError(
                "`charge` must be specified if `species` is not None."
                )

    @property
    def bin_edges(self):
        """Bin edge coordinates for each phase-space dimension.

        Returns
        -------
        dict[str, np.ndarray]
            Read-only mapping from dimension name to bin edge array.
        """
        return MappingProxyType(self._bin_edges)

    @property
    def dimensions(self):
        """Bin center coordinates for each phase-space dimension.

        Returns
        -------
        dict[str, np.ndarray]
            Read-only mapping from dimension name to bin center array.
        """
        dimensions = {}
        for key, val in self.bin_edges.items():
            dimensions[key] = ( val[1:] + val[:-1] ) / 2
        return dimensions

    @property
    def distribution(self):
        """Phase-space density values.

        Returns
        -------
        np.ndarray
            The stored distribution array.
        """
        return self._distribution

    @property
    def histogram(self):
        """Distribution represented as histogram values (i,e, particles) per
        bin.

        Returns
        -------
        np.ndarray
            Histogram values obtained by multiplying the distribution
            by the phase-space bin volume.
        """
        return self._distribution * self.phasespacevolume

    @property
    def phasespacevolume(self):
        """Return the volume of the phase space.

        Returns
        -------
        volume : float
            Volume of the phase space.
        """
        volume = 1
        for dim in self.bin_edges.values():
            dv = dim[1:] - dim[:-1]
            volume = np.multiply.outer(volume, dv)
        return volume

    @property
    def species(self):
        return self._species

    @property
    def charge(self):
        return self._charge

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

    def copy(self):
        """Make a (deep) copy of this object.

        Returns
        -------
        copy : :class:`DistData`
            Copy of this object.
        """
        return Dist(
            histogram=self.histogram,
            bin_edges=self.bin_edges,
            species=self._species,
            charge=self._charge,
        )

    def multiply(
        self,
        multiplier: np.ndarray | float,
        dimensions: Optional[Sequence[str]] = None,
        ):
        """Multiply the distribution in place.

        Parameters
        ----------
        multiplier : array_like or scalar
            The value to multiply the distribution by.

            This can either be a scalar, an array with the same shape as
            the distribution or an array with .
        dimensions : sequence of str, optional
            Names of the dimensions corresponding to the axes of `value`
            when `value.ndim < self.distribution.ndim`.

        Raises
        ------
        ValueError
            If dimensions are inconsistent with the shape of `value`
            or the distribution.
        """
        distribution = self.distribution
        multiplier = np.asarray(multiplier)
        if (
            (multiplier.size == 1 or multiplier.shape == distribution.shape)
             and axes is None
             ):
            self._distribution *= multiplier
            return
        elif axes is None:
            raise ValueError("axes must be specified if b is not a scalar")

        if len(axes) != multiplier.ndim:
            raise ValueError(
                "The number of `dimensions` must match the number of dimensions"
                " of `multiplier`."
                )

        if distribution.ndim < multiplier.ndim:
            raise ValueError("b.ndim cannot exceed a.ndim")

        axes = []
        for i, dimension in enumerate(self.dimensions):
            if dimension in dimensions:
                axes.append(i)

        if len(axes) != len(dimensions):
            raise ValueError(
                "asdadad"
            )

        axes = tuple(axes)
        for bdim, ax in enumerate(axes):
            if multiplier.shape[bdim] != distribution.shape[ax]:
                raise ValueError(
                    f"Shape mismatch: multiplier.shape[{bdim}]={multiplier.shape[bdim]} "
                    f"!= a.shape[{ax}]={distribution.shape[ax]}"
                )

        shape = [1] * distribution.ndim
        for bdim, ax in enumerate(axes):
            shape[ax] = multiplier.shape[bdim]

        self._distribution *= multiplier.reshape(shape)

    def slice(self, slices):
        """Slice distribution.

        Parameters
        ----------
        slices : dict[str, int or slice]
            Dictionary mapping dimension names to slices.
        """
        bin_edges = list(self.bin_edges)

        unknown = set(slices) - set(bin_edges)
        if unknown:
            raise KeyError(
                f"Unknown dimension(s): {unknown}. The dimensions of "
                f"the distribution are: {bin_edges}"
                )

        new_bin_edges = {}
        for dimension, edges in self.bin_edges.items():
            if dimension not in slices:
                new_bin_edges[dimension] = edges
                continue

            sel = slices[dimension]
            if not isinstance(sel, (int, np.integer, slice)):
                raise TypeError(
                    f"Slice for dimension {dimension!r} must be int or "
                    f"slice, got {type(sel).__name__}"
                )

            if isinstance(sel, int):
                sel = slice(sel, sel + 1) if sel >= 0 else slice(sel, None)
                slices[dimension] = sel

            start, stop, step = sel.indices(len(edges) - 1)
            new_bin_edges[dimension] = edges[start : stop + 1 : step]
            if (start < 0 or stop < 0 or start >= len(edges) - 1 or stop >= len(edges) - 1):
                raise ValueError(
                    f"Requested slice is out of bounds for dimension {dimension!r}"
                )

        idx = tuple(slices.get(edges, slice(None)) for edges in bin_edges)
        self._distribution = self.distribution[idx]
        self._bin_edges = new_bin_edges

    def integrate(self, dims: dict[str, np.ndarray | None]) -> None:
        """Integrate distribution over one or more dimensions in-place.

        Parameters
        ----------
        dims : dict[str, ndarray | None]
            Mapping from dimension name to integration weights.

            - ``{"x": None}`` integrates over the full ``x`` dimension.
            - ``{"x": weights}`` performs a weighted integral, where
            ``weights`` has length equal to the number of bins in ``x``.

        Raises
        ------
        KeyError
            If a requested dimension does not exist.
        ValueError
            If a weight array has an invalid shape.
        """
        dist = self._distribution
        bin_edges = dict(self._bin_edges)

        # Determine axes in the current ordering
        dim_names = list(bin_edges)

        # Integrate highest axis first so axis numbers remain valid
        for dim in sorted(dims, key=lambda d: dim_names.index(d), reverse=True):
            if dim not in bin_edges:
                raise KeyError(f"Unknown dimension '{dim}'.")

            axis = dim_names.index(dim)
            edges = bin_edges[dim]

            dx = edges[1] - edges[0]
            nbin = len(edges) - 1

            weights = dims[dim]

            if weights is None:
                factor = dx * np.ones((nbin,))
            else:
                weights = np.asarray(weights)

                if weights.shape != (nbin,):
                    raise ValueError(
                        f"Weights for dimension '{dim}' must have shape "
                        f"({nbin},), got {weights.shape}."
                    )

                factor = weights * dx

            shape = [1] * dist.ndim
            shape[axis] = nbin

            dist = np.sum(
                dist * np.asarray(factor).reshape(shape),
                axis=axis,
            )

            del bin_edges[dim]
            dim_names.pop(axis)

        self._distribution = dist
        self._bin_edges = bin_edges


    def interpolate(self, coords, interpolate_hist=False):
        """Interpolate distribution using linear interpolation.

        The linear interpolation is performed using the values of the
        distribution at the bin centers

        Parameters
        ----------
        coords : dict[str, float]
            Dictionary with dimensions and coordinates at which the distribution
            is to be interpolated.

            It is not necessary to give coordinates for each dimensions. For
            example, if interpolating only in one dimension when the
            distribution is 3D, the result will be 2D array where all values are interpolated in the queried
            dimension only.

        interpolate_hist : bool, optional
            If ``True``, the histogram is interpolated instead.

        Returns
        -------
        out: array_like
            Interpolated value.
        """
        bcoords = np.broadcast_arrays(
            *[np.asarray(coord) for coord in coords.values()]
            )
        query_shape = bcoords[0].shape if bcoords else ()
        broadcasted_coords = {
            name: coord
            for name, coord in zip(coords.keys(), bcoords)
        }

        result = self.distribution if not interpolate_hist else self.histogram
        units = result.units
        result = result.v

        remaining_dimensions = {
            name: axis_index
            for axis_index, name in enumerate(self.bin_edges.keys())
            }

        for name, x in broadcasted_coords.items():
            grid = np.asarray(self.bin_edges[name])
            xmin, xmax = grid[[0, -1]]
            if np.any((x < xmin) | (x > xmax)):
                raise ValueError(
                    f"Coordinate outside range for '{name}': [{xmin}, {xmax}]"
                )

            centers = 0.5 * (grid[:-1] + grid[1:])
            u = (x - centers[0]) / (centers[1] - centers[0])

            i0 = np.floor(u).astype(int)
            i0 = np.clip(i0, 0, len(centers) - 2)

            i1 = i0 + 1

            x0 = centers[i0]
            x1 = centers[i1]

            w = (x - x0) / (x1 - x0)

            axis_index = remaining_dimensions[name]
            v0 = np.take(result, i0, axis=axis_index)
            v1 = np.take(result, i1, axis=axis_index)

            extra_dims = v0.ndim - len(query_shape)
            w_exp = w.reshape(query_shape + (1,) * extra_dims)
            result = (1.0 - w_exp) * v0 + w_exp * v1

            del remaining_dimensions[name]
            for dimension, original_index in remaining_dimensions.items():
                if original_index > axis_index:
                    remaining_dimensions[dimension] -= 1

        return result * units


    def interpolate_points(self, points):
        """Multilinear interpolation on bin-center distribution with
        linear extrapolation to bin edges (inside domain only).
        """

        pts = np.asarray(points, dtype=float)
        single = pts.ndim == 1
        if single:
            pts = pts[None, :]

        ndim = pts.shape[1]
        if ndim != len(self._bin_edges):
            raise ValueError("Dimensionality mismatch.")

        keys = list(self._bin_edges.keys())

        centers = []
        edges = []

        # build centers + edge extents
        for k in keys:
            e = np.asarray(self._bin_edges[k])
            h = e[1] - e[0]  # uniform spacing guaranteed

            c = 0.5 * (e[:-1] + e[1:])
            centers.append(c)
            edges.append((e[0], e[-1], h))

        data = self._distribution

        # --- extend grid by 1 cell on each side per axis ---
        grid = list(centers)
        extended = data

        for axis, (lo, hi, h) in enumerate(edges):

            # linear extrapolation slopes along axis
            sl_left = (np.take(extended, 1, axis=axis) -
                    np.take(extended, 0, axis=axis)) / h
            sl_right = (np.take(extended, -1, axis=axis) -
                        np.take(extended, -2, axis=axis)) / h

            left_slice = np.take(extended, 0, axis=axis) - sl_left * h
            right_slice = np.take(extended, -1, axis=axis) + sl_right * h

            extended = np.concatenate(
                [np.expand_dims(left_slice, axis=axis),
                extended,
                np.expand_dims(right_slice, axis=axis)],
                axis=axis
            )

            # extend coordinate grid accordingly
            c = grid[axis]
            grid[axis] = np.concatenate([[c[0] - h], c, [c[-1] + h]])

        interp = RegularGridInterpolator(
            grid,
            extended,
            method="linear",
            bounds_error=True
        )

        out = interp(pts)
        return out[0] if single else out


class Hist():
    """Interface to histogram diagnostics data.

    This class provides access to the histogram data, whether it is stored in
    a file or in memory.
    """

    def __init__(self):
        self._file = None
        self._cdata = None
        self._params = None

    @property
    def dimensions(self) -> dict[str, unyt.unyt_array]:
        """Get histogram dimensions."""
        out = {}
        for dim in self._params.dimensions:
            out[dim.name] = unyt.unyt_array(np.linspace(
                    dim.min, dim.max, dim.bins + 1,
                    DIMENSIONS[dim.name]
                    ))
        return out

    @property
    def values(self):
        if self._file is None:
            return self._cdata.values_ref
        return self._file.read("values")

    def save(self, file):
        """"""
        file.write("values", self.values)
        del self._cdata
        self._file = file

    def extract(
            self,
            dimensions:
            Sequence[str],
            species: Optional[str]=None,
            charge: Optional[int]=None,
            ) -> Dist:
        """Extract histogram data as a distribution.

        Parameters
        ----------
        dimensions : Sequence[str]
            Dimensions in the resulting distribution.

            If the collected histogram has more dimensions than requested here,
            the histogram is summed over the dimensions that were not requested.
            The dimensions of the output distribution have the same order as the
            query.

        Returns
        -------
        dist : Distribution
            The extracted distribution.
        """
        hist_dimensions = self.dimensions
        hist_names = tuple(hist_dimensions.keys())
        invalid = [
            dim for dim in dimensions
            if dim not in hist_names
        ]
        if invalid:
            raise ValueError(
                f"Unknown dimension(s): {invalid}. Possible values are: "
                f"{list(hist_names)}.",
                )

        keep = [hist_names.index(dim) for dim in dimensions]
        if species is not None:
            keep.append(hist_names.index("charge"))
        reduce_axes = tuple(i for i in range(len(hist_names)) if i not in keep)
        projected = self.values.sum(axis=reduce_axes)
        if species is not None:
            if len(hist_dimensions["charge"].nbin) == 1:
                charge = hist_dimensions["charge"].min
                idx = 0
            projected = projected[..., idx]

        remaining_original = [dim for dim in hist_names if dim in dimensions]
        permutation = [remaining_original.index(dim) for dim in dimensions]
        if permutation != list(range(len(dimensions))):
            projected = np.transpose(projected, permutation)
        else:
            projected = np.transpose(projected)

        projected_axes = {dim: hist_dimensions[dim] for dim in dimensions}
        return Dist(projected, projected_axes, species=species, charge=charge)

    @classmethod
    def from_params(cls, params: HistParams) -> "Hist":
        """Initialize histogram diagnostics from simulation parameters.

        Parameters
        ----------
        params : HistParams
            Histogram diagnostics parameters.
        """
        cdata = Struct()
        idx, dimensions = [], list(DIMENSIONS.keys())
        for i, dimension in enumerate(params.dimensions):
            idx.append(dimensions.index(dimension.name))

        cdata.nbin, dims = 1, []
        max_dim = len(DIMENSIONS)
        cdata.strides[max_dim - 2] = max(cdata.dimensions[max_dim-1].n, 1)
        for i in np.flip(np.arange(max_dim)):
            cdata.dimensions[i].min = 0
            cdata.dimensions[i].max = 1
            cdata.dimensions[i].n = 0
            for j, coord in enumerate(idx):
                if i == coord:
                    nbins = params.dimensions[j].bins
                    cdata.dimensions[i].min = params.dimensions[j].min.v
                    cdata.dimensions[i].max = params.dimensions[j].max.v
                    cdata.dimensions[i].n = nbins
                    cdata.nbin *= nbins
                    dims.append(nbins)
            if i < max_dim - 2:
                cdata.strides[i] = max(cdata.dimensions[i+1].n, 1) * cdata.strides[i + 1]

        cdata.values_ref = np.zeros(np.flip(dims))
        cdata.values = cdata.values_ref.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        cdata.values_ref.flags.writeable = False

        hist = cls()
        hist._cdata = cdata
        hist._params = params
        return hist
