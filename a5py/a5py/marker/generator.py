"""
Generate marker (and some trivial particle) distributions.

File: marker/generator.py
"""
import numpy as np
import scipy.constants as const

from scipy.interpolate import RegularGridInterpolator

import a5py.marker.phasespace as phasespace
from a5py.ascotpy.ascotpy import Ascotpy

import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

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
    tol = 1e-2
    err = tol+1
    newdist = np.zeros((rgrid.size-1, zgrid.size-1, kgrid.size-1))
    while(err > tol):
        olddist = np.copy(newdist)

        if rhoksidist is not None:
            rhogrid = rhoksidist[0]
            ksigrid = rhoksidist[1]

            rhovals = rhogrid[0] + np.random.rand(Ndraw) \
                      * (rhogrid[-1]-rhogrid[0])
            ksivals = ksigrid[0] + np.random.rand(Ndraw) \
                      * (ksigrid[-1]-ksigrid[0])

            weights = np.zeros(rhovals.shape)
            for i in range(rhogrid.size-1):
                for j in range(ksigrid.size-1):
                    ids = np.logical_and.reduce(
                        [ rhogrid[i] < rhovals, rhovals < rhogrid[i+1],
                          ksigrid[j] < ksivals, ksivals < ksigrid[j+1] ] )
                    weights[ids] = rhoksidist[2][i,j]

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
    PDT     = 5.6e-13*n1*n2*sigmaDT / 1e6

    # Remove weight outside separatrix
    #rho = a5.evaluate(R=R, phi=0, z=Z, t=0, quantity="rho").reshape(
    #    rgrid.size-1, zgrid.size-1, ksigrid.size-1)
    #PDT[rho > 1] = 0

    # Convert power density to total alpha power in each cell
    PDT     = (PDT * 2*np.pi*R * dr * dz / (ksigrid.size-1))

    # Convert power to particle birth rate and replace NaNs with zeros.
    dist = PDT / (3.5e6*const.e)
    dist[np.isnan(dist)] = 0

    if free:
        a5.free(bfield=True, plasma=True)

    if plot:
        plotdist((rgrid, zgrid, ksigrid, dist))

    return (rgrid, zgrid, ksigrid, dist)


def generate(mass, charge, energy, markerdist, particledist, Nmrk,
             minweight=1e-3, plot=False):
    """
    Generate weighted markers from marker and particle distributions.

    Right now only guiding centers are generated but this can be fixed if
    necessary.

    Args:
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        markerdist : tuple <br>
            Tuple (rgrid,zgrid,ksigrid,ordinate) for marker distribution.
        particledist : tuple <br>
            Tuple (rgrid,zgrid,ksigrid,ordinate) for particle distribution.
        Nmrk : int <br>
            Number of markers to be generated.
        minweight : float, optional <br>
            Markers with weight < minweight * average_weight are rejected.

    Returns:
        Guiding center marker data on a dictionary and final marker
        distribution.
    """
    rgrid   = markerdist[0]
    zgrid   = markerdist[1]
    ksigrid = markerdist[2]

    # Interpolate particle distribution to marker distribution grid.
    linint = RegularGridInterpolator(
        (particledist[0][:-1] + (particledist[0][1] - particledist[0][0]) / 2,
         particledist[1][:-1] + (particledist[1][1] - particledist[1][0]) / 2,
         particledist[2][:-1] + (particledist[2][1] - particledist[2][0]) / 2),
        particledist[3], method='linear', fill_value=None)

    R, Z, Ksi  = np.meshgrid(rgrid[:-1], zgrid[:-1], ksigrid[:-1],
                             indexing='ij')
    R   = R.ravel()
    Z   = Z.ravel()
    Ksi = Ksi.ravel()

    particledist = linint( np.array([
        R   + (rgrid[1]   - rgrid[0]  ) / 2,
        Z   + (zgrid[1]   - zgrid[0]  ) / 2,
        Ksi + (ksigrid[1] - ksigrid[0]) / 2]).transpose() )

    # Turn marker distribution into 1D array that holds cumulative
    # probability P(i).
    markerdist = np.append( [0], np.cumsum(markerdist[3].ravel()) )

    # Markers are generated in a sets of Nsample markers. Once the total number
    # of markers is equal or greater than Nmrk, markers are shuffled and the
    # requested number of markers are drawn.
    Nsample = 10000

    indices = np.array([])
    Ntotal = indices.size
    while Ntotal < Nmrk:
        # Draw Nsample uniformly distributed (between 0 and 1) random numbers
        # and choose corresponding cells by finding last index where P is
        # smaller than the drawn number.
        prob = np.random.rand(Nsample)
        idx  = np.searchsorted(markerdist, prob, side='left')

        # Get weights for each marker
        weights = particledist[idx]

        # Remove those markers whose weight would be too small.
        idx = idx[weights > minweight * np.nanmean(weights)]

        indices = np.append(indices, idx)
        Ntotal = indices.size

    # Generate markers within these cells. Each marker represents a phase
    # space element enclosed inside the cell, so actual marker coordinates
    # are drawn from uniform distribution where limits are the cell edges.
    indices = indices.astype(int)
    r   = R[indices]   + np.random.rand(Ntotal) * (rgrid[1]   - rgrid[0]  )
    z   = Z[indices]   + np.random.rand(Ntotal) * (zgrid[1]   - zgrid[0]  )
    ksi = Ksi[indices] + np.random.rand(Ntotal) * (ksigrid[1] - ksigrid[0])

    # Weight is divided equally among markers within a given cell
    weights = particledist[indices]
    _, idx, ncount = np.unique(indices, return_inverse=True, return_counts=True)
    weights = weights/ncount[idx]

    # Shuffle and draw Nmrk markers
    ids = np.arange(Ntotal)
    np.random.shuffle(ids)
    r   = r[ids[:Nmrk]]
    z   = z[ids[:Nmrk]]
    ksi = ksi[ids[:Nmrk]]
    weights = weights[ids[:Nmrk]]

    # Store marker data
    mass   = mass / const.physical_constants["atomic mass constant"][0]
    charge = np.around(charge/const.e)
    anum   = np.around(mass)
    znum   = charge
    energy = energy / const.e
    markers = {
        "n"      : Nmrk,
        "ids"    : np.arange(Nmrk) + 1,
        "mass"   : mass * np.ones((Nmrk,1)),
        "charge" : charge * np.ones((Nmrk,1)),
        "r"      : r,
        "phi"    : 360 * np.random.rand(Nmrk),
        "z"      : z,
        "energy" : energy * np.ones((Nmrk,1)),
        "pitch"  : ksi,
        "zeta"   : 2*np.pi * np.random.rand(Nmrk),
        "anum"   : anum * np.ones((Nmrk,1)),
        "znum"   : znum * np.ones((Nmrk,1)),
        "weight" : weights,
        "time"   : 0 * np.ones((Nmrk,1))
    }

    dist = np.histogramdd((r,z,ksi), bins=(rgrid,zgrid,ksigrid),
                          weights=weights)[0]
    dist = (rgrid,zgrid,ksigrid,dist)
    if plot:
        plotdist(dist)

    return (markers, dist)


def plotdist(dist, axes=None):
    """
    Plot distribution in Rz and in ksi.
    """

    if axes is None:
        fig = plt.figure()

        gs = GridSpec(2,1, height_ratios=[3, 1])
        axesrz  = fig.add_subplot(gs[0,0])
        axesksi = fig.add_subplot(gs[1,0])
    else:
        axesrz  = axes[0]
        axesksi = axes[1]

    if axesrz is not None:
        mesh = axesrz.pcolormesh(dist[0], dist[1],
                                 np.sum(dist[3], axis=2).transpose())
        axesrz.set_aspect('equal', adjustable='box')
        plt.colorbar(mesh, ax=axesrz)

    if axesksi is not None:
        axesksi.bar(dist[2][:-1], np.sum(dist[3], axis=(0,1)),
                    align='edge', width=(dist[2][1] - dist[2][0]))

    if axes is None:
        plt.show(block=False)
