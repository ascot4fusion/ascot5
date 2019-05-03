"""
Populate phase space and map it to real space.

File: transport/phasespace.py
"""
import numpy as np
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy.ascotpy import Ascotpy

def istrapped(a5, mass, charge, energy, pitch, P, mu, rmin):
    """
    Check whether a particle is trapped.
    """

    # Make R grid where to evaluate B and psi
    axis  = a5.evaluate(0, 0, 0, 0, "axis")
    rgrid = np.linspace(rmin, axis["r"], 100)

    bnorm = a5.evaluate(rgrid, 0, axis["z"], 0, "bnorm")
    psi   = a5.evaluate(rgrid, 0, axis["z"], 0, "psi")

    vnorm = np.sqrt(2*energy/mass)

    trapped = np.zeros(mu.shape) == 1
    for i in range(mu.size):
        val1  = vnorm*vnorm - 2*mu*bnorm/mass
        if all(val1 < 0):
            trapped[i] = False
            continue

        val2 = charge * psi + np.sign(pitch) * mass * rgrid * np.sqrt(val1)
        val2[val1 < 0] = P

        if any(val - P > 0) and any(val - P < 0):
            trapped[i] = False
        else:
            trapped[i] = True

    return trapped


def initgrid(a5, mass, charge, energy, rhomin, rhomax, nP, nmu, nrho, nksi,
             padding):
    """
    Create grids.
    """

    rhogrid = np.linspace(rhomin, rhomax, nrho)
    rz   = a5.get_rhotheta_rz( rhogrid, 0, 0, 0 )
    bmin = a5.evaluate(rz[0][-1], 0, rz[1][-1], 0, "bnorm")

    vnorm = np.sqrt(2*energy/mass)
    mumax = mass * vnorm * vnorm / (2*bmin)
    mumin = -mumax

    mugrid = np.linspace(mumin*(1+padding), mumax*(1+padding), nmu)

    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")

    mu    = mugrid
    psi   = psi * np.ones(mu.shape)
    r     = rz[0] * np.ones(mu.shape)
    vpar2 = vnorm * vnorm - 2 * bnorm * np.absolute(mu) / mass
    vpar  = np.sign(mu)*np.sqrt(vpar2)

    P  = mass * r * vpar + charge * psi

    Pmin = np.nanmin(P.ravel())
    Pmax = np.nanmax(P.ravel())

    if Pmin < 0:
        Pmin = Pmin*(1+padding)
    else:
        Pmin = Pmin*(1-padding)

    if Pmax < 0:
        Pmax = Pmax*(1-padding)
    else:
        Pmax = Pmax*(1+padding)

    Pgrid = np.linspace(Pmin, Pmax, nP)

    ksigrid = np.linspace(-1, 1, nksi)

    return (Pgrid, mugrid, rhogrid, ksigrid)


def maprhomu2Pmu(a5, mass, charge, energy, rhovals, muvals,
                 Pgrid, mugrid, weights=None):
    """
    Map rho_omp-mu values to a P-mu grid.
    """
    #rhogrid = np.linspace(0, 1, 100)
    rz    = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    vnorm = np.sqrt(2*energy/mass)

    vpar   = np.sign(muvals) * np.sqrt(vnorm*vnorm
                                       - 2*bnorm*np.absolute(muvals)/mass)
    Pvals  = mass * rz[0] * vpar + charge * psi

    if weights is None:
        weights = np.ones(muvals.shape)

    markers = np.zeros( (Pgrid.size-1, mugrid.size-1) )
    for i in range(Pgrid.size-1):
        for j in range(mugrid.size-1):
            ids = np.logical_and.reduce( [
                Pvals >= Pgrid[i], Pvals < Pgrid[i+1],
                muvals >= mugrid[j], muvals < mugrid[j+1]
            ] )

            markers[i,j] = np.sum(weights[ids])

    return markers


def mapPmu2rzk(a5, mass, charge, energy, Pvals, muvals,
               rgrid, zgrid, ksigrid, weights=None):
    """
    Map mu-rho_omp points
    """
    #Pvals, muvals = np.meshgrid(Pgrid, mugrid)

    vnorm = np.sqrt(2*energy/mass)
    if weights is None:
        weights = np.ones(muvals.shape)

    Pvals   = Pvals.ravel()
    muvals  = muvals.ravel()
    weights = weights.ravel()
    #musign  = np.sign(muvals)

    R, Z, Ksi = np.meshgrid(rgrid, zgrid, ksigrid, indexing='ij')
    psi   = a5.evaluate(R, 0, Z, 0, "psi").reshape(R.shape)
    bnorm = a5.evaluate(R, 0, Z, 0, "bnorm").reshape(R.shape)
    Mu = np.sign(Ksi)*mass * (1-Ksi*Ksi) * vnorm * vnorm / (2*bnorm)
    P  = mass * R * Ksi * vnorm  + charge * psi

    markers = np.zeros( (R.shape[0]-1,R.shape[1]-1,R.shape[2]-1) )
    for i in range(R.shape[0]-1):
        for j in range(R.shape[1]-1):
            for k in range(R.shape[2]-1):
                a1 = (Mu[i,j,k]       - muvals) >= 0
                a2 = (Mu[i+1,j,k]     - muvals) >= 0
                a3 = (Mu[i,j+1,k]     - muvals) >= 0
                a4 = (Mu[i,j,k+1]     - muvals) >= 0
                a5 = (Mu[i+1,j+1,k]   - muvals) >= 0
                a6 = (Mu[i+1,j,k+1]   - muvals) >= 0
                a7 = (Mu[i,j+1,k+1]   - muvals) >= 0
                a8 = (Mu[i+1,j+1,k+1] - muvals) >= 0

                a = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8
                a = np.logical_and.reduce([a > 0, a < 8])

                b1 = (P[i,j,k]       - Pvals) >= 0
                b2 = (P[i+1,j,k]     - Pvals) >= 0
                b3 = (P[i,j+1,k]     - Pvals) >= 0
                b4 = (P[i,j,k+1]     - Pvals) >= 0
                b5 = (P[i+1,j+1,k]   - Pvals) >= 0
                b6 = (P[i+1,j,k+1]   - Pvals) >= 0
                b7 = (P[i,j+1,k+1]   - Pvals) >= 0
                b8 = (P[i+1,j+1,k+1] - Pvals) >= 0

                b = b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8
                b = np.logical_and.reduce([b > 0, b < 8])

                #c1 = np.sign(Ksi[i,j,k])   * musign >= 0
                #c2 = np.sign(Ksi[i,j,k+1]) * musign >= 0
                #c  =  np.logical_or.reduce([c1, c2])

                markers[i,j,k] = np.sum(
                    weights[ np.logical_and.reduce([a, b]) ].ravel() )

    return markers


def maprzk2Pmu(a5, mass, charge, energy, rvals, zvals, ksivals,
               Pgrid, mugrid, weights=None):
    """
    """
    vnorm = np.sqrt(2*energy/mass)
    bnorm = a5.evaluate(rvals, 0, zvals, 0, "bnorm")
    psi   = a5.evaluate(rvals, 0, zvals, 0, "psi")

    muvals = np.sign(ksivals) * (1-ksivals*ksivals) * vnorm * vnorm * mass / (2*bnorm)
    Pvals  = mass * rvals * ksivals * vnorm + charge * psi

    if weights is None:
        weights = np.ones(muvals.shape)

    markers = np.zeros( (Pgrid.size-1, mugrid.size-1) )
    for i in range(Pgrid.size-1):
        for j in range(mugrid.size-1):
            ids = np.logical_and.reduce( [
                Pvals >= Pgrid[i], Pvals < Pgrid[i+1],
                muvals >= mugrid[j], muvals < mugrid[j+1]
            ] )

            markers[i,j] = np.sum(weights[ids])

    return markers


def mapPmu2rhomu():
    """
    """
    pass


def check():
    """
    """
    pass


def populate(fn, mass, charge, energy, rhomin, rhomax, nrho, nmu,
             rgrid, zgrid, ksigrid, weights, Nmrk):
    """
    Fill phase space with markers.
    """
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    nP  = 100
    nmu = 50
    nksi = 20
    nrho = 10

    try:
        Pgrid, mugrid, rhogrid, pitchgrid = initgrid(a5, mass, charge, energy,
                                                     rhomin, rhomax,
                                                     nP, nmu, nrho, nksi, 1e-3)

        rhovals = rhogrid[:-1] + (rhogrid[1] - rhogrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2

        rhovals, muvals = np.meshgrid(rhovals, muvals, indexing='ij')
        markers = np.ones(rhovals.shape)
        markers = maprhomu2Pmu(a5, mass, charge, energy,
                               rhovals.ravel(), muvals.ravel(),
                               Pgrid, mugrid, weights=markers.ravel())

        fig  = plt.figure()
        axes = fig.add_subplot(2,4,1)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = mapPmu2rzk(a5, mass, charge, energy,
                             Pvals.ravel(), muvals.ravel(),
                             rgrid, zgrid, ksigrid, weights=markers.ravel())

        axes = fig.add_subplot(2,4,2)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(markers, axis=2).transpose())
        plt.colorbar(h, ax=axes)

        ksivals = ksigrid[:-1] + (ksigrid[1] - ksigrid[0])/2
        axes = fig.add_subplot(2,4,3)
        h = axes.plot(ksivals, np.sum(np.sum(markers, axis=0), axis=0))

        dist = 1/markers
        dist[weights < 1]    = 0
        dist[np.isnan(dist)] = 0
        dist[np.isinf(dist)] = 0

        n = Nmrk
        while(np.sum(np.floor(n*dist.ravel()))) < Nmrk:
            n = n +1

        dist = np.floor(dist*n)
        dist[dist < 1] = 0

        axes = fig.add_subplot(2,4,4)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(dist, axis=2).transpose())
        plt.colorbar(h, ax=axes)

        n   = int(np.sum(dist.ravel()))
        r   = np.zeros((n,1))
        z   = np.zeros((n,1))
        ksi = np.zeros((n,1))
        w   = np.zeros((n,1))
        n0  = 0
        for i in range(dist.shape[0]):
            for j in range(dist.shape[1]):
                for k in range(dist.shape[2]):
                    n = int(dist[i,j,k])
                    if n == 0:
                        continue

                    r[n0:n0+n]   = rgrid[i]   + np.random.rand(n,1)*(
                        rgrid[1] - rgrid[0])
                    z[n0:n0+n]   = zgrid[j]   + np.random.rand(n,1)*(
                        zgrid[1] - zgrid[0])
                    ksi[n0:n0+n] = ksigrid[k] + np.random.rand(n,1)*(
                        ksigrid[1] - ksigrid[0])
                    w[n0:n0+n]   = weights[i,j,k] / n
                    n = n0+n

        n = r.size
        ids = np.random.shuffle(np.arange(n))
        r   = np.squeeze(r[ids][:Nmrk])
        z   = np.squeeze(z[ids][:Nmrk])
        ksi = np.squeeze(ksi[ids][:Nmrk])
        w   = np.squeeze(w[ids][:Nmrk] * (Nmrk/n))

        markers = maprzk2Pmu(a5, mass, charge, energy, r, z, ksi,
                             Pgrid, mugrid, weights=None)

        print(markers)

        axes = fig.add_subplot(2,4,5)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


def evaluatecoefs():
    """
    Evaluate transport coefficients.
    """
    pass


def run2dmodel():
    """
    """
    pass
