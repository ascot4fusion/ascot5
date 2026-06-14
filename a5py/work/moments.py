"""Functions that calculate moments of the distribution.
"""
import numpy as np

def density(dist, volume):
    """Calculate number density.

    Parameters
    ----------
    dist : :class:`DistData`
        Distribution from where the moments are calculated.
    moment : class:`DistMoment`
        Moment data used in calculation and where the result is stored.
    """
    integrate = {}
    if "rho" in dist.dimensions:
        integrate["rho"] = None
        if "theta" in dist.dimensions:
            integrate["theta"] = None
        if "phi" in dist.dimensions:
            integrate["phi"] = None
        if "r" in dist.dimensions or "z" in dist.dimensions:
            raise ValueError()
    else:
        if "r" in dist.dimensions:
            integrate["r"] = None
        if "phi" in dist.dimensions:
            integrate["phi"] = None
        if "z" in dist.dimensions:
            integrate["z"] = None
        if len(integrate) == 0:
            raise ValueError()

    momentum_space = set(dist.dimensions) - set(integrate.keys())
    integrate_these = {}
    for k in momentum_space:
        integrate_these[k] = None
    for i, dim in enumerate(integrate):
        if dist.dimensions[dim].size != volume.shape[i]:
            print(dist.dimensions[dim].size, volume.shape[i])
            raise ValueError()

    dist = dist.copy()
    dist.integrate(integrate_these)
    return dist.histogram / volume
