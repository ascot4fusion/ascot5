"""
Basic operations for processing distributions.

File: basic.py
"""

import numpy     as np
import importlib

plt = importlib.util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

from scipy.interpolate import griddata, RectBivariateSpline

def histogram2density(dist):
    """
    Convert histogram to density.

    A histogram tells how many markers there are in a cell whereas density,
    i.e. histogram divided by cell volume, is a physical quantity.

    The density will not necessarily be in SI units, as the volume element is
    defined by the abscissae as vol = (a1[1:] - a1[:-1]) x (a2[1:] - a2[:-1])
    x ... .

    Args:
        dist: dict_like <br>
            Distribution whose histogram ordinate is converted to a density.
    Returns:
        Reference to the dist argument.
    """
    if "density" in dist:
        print("Distribution ordinate is already density.")
        return dist

    # Calculate volumes for each cell
    vol = 1
    for coord in dist["abscissae"]:
        edges = dist[coord + '_edges']
        dv    = (edges[1:] - edges[:-1])
        vol   = np.multiply.outer(vol, dv)

    dist["density"] = dist["histogram"] / vol
    del dist["histogram"]
    return dist

def squeeze(dist, **kwargs):
    """
    Integrate or slice dimensions.

    Integration means a sum sum(dist[a:b,...] * vol) will be calculated and
    the dimensionality of the distribution is reduced by one. Slicing means
    that abscissa is cut as abscissa[a:b] and the ordinate is modified
    accordingly.

    This function operates recursively on each given dimension.

    Args:
        dist : dict_like <br>
            Distribution so be integrated or sliced.
        kwargs <br>
            Name(s) (R, phi, z, vpa, vpe, time, charge) of those dimensions
            along which the distribution is either sliced or integrated over.
            Names are given as keyword value pairs, {str : tuple} or
            {str : array_like}, where a tuple value (a, b) are indices for
            slicing and array [a, b] are indices for integration. A scalar zero
            means whole dimension is integrated.
    Returns:
        Reference to the dist argument.
    """

    # Are kwargs included?
    if len(kwargs.keys()) == 0:
        return dist

    coord  = list(kwargs.keys())[0]
    slices = kwargs[coord]
    del kwargs[coord]

    # Call recursively if additional dimensions still exist in kwargs.
    if len(kwargs.keys()) > 0:
        squeeze(dist, **kwargs)

    dim = dist["abscissae"].index(coord)

    if type(slices) is tuple:
        dist[coord] = dist[coord][slices[0]:slices[1]]
        dist[coord + '_edges'] = dist[coord][slices[0]:slices[1]+1]
        dist['n_' + coord] = dist[coord].size

        # Slice given dimension
        #indx = [Ellipsis] * len(dist["abscissae"])
        #indx[dim] = slice(slices[0], slices[1], 1)
        #dist["density"] = dist["density"][indx]
        dist["density"] = np.take(dist["density"], range(slices[0], slices[1]),
                                  axis=dim)

    if type(slices) is int:
        slices = [0, dist["n_" + coord]]

    if type(slices) is list or not slices:
        edges   = dist[coord + "_edges"]
        weights = (edges[1:] - edges[0:-1])

        if slices is not None:
            mask = np.zeros(weights.shape)
            mask[slices[0]:slices[1]] = 1
            weights = weights*mask

        dist["density"], s = np.average(dist["density"],
                                        axis=dim, weights=weights,
                                        returned=True)
        dist["density"] *= s
        del dist[coord]
        del dist[coord + '_edges']
        del dist['n_' + coord]
        del dist["abscissae"][dim]

    return dist
