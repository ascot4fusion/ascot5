"""
Populate phase space and map it to real space.

File: transport/phasespace.py
"""
import numpy as np
import scipy.constants as const
from itertools import compress
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy.ascotpy import Ascotpy


def evalPmu(a5, mass, charge, energy, r, z, ksi):
    """
    Evaluate P and mu at the given coordinates.
    """
    psi   = a5.evaluate(r, 0, z, 0, "psi")
    bnorm = a5.evaluate(r, 0, z, 0, "bnorm")
    vnorm = np.sqrt(2*energy/mass)

    mu = (1 - ksi*ksi) * mass * vnorm * vnorm / (2*bnorm)
    P  = mass * r * vnorm * ksi + charge * psi

    return (P, mu)


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


def mapPmu2rhomu(a5, mass, charge, energy, Pvals, muvals, rhogrid, mugrid,
                 weights=None):
    """
    Map P-mu values to a rho_omp-mu grid.
    """
    rz    = a5.get_rhotheta_rz( rhogrid, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    vnorm = np.sqrt(2*energy/mass)

    mu    = mugrid.transpose()
    psi   = psi.reshape(-1,1) * np.ones(mu.shape)
    bnorm = bnorm.reshape(-1,1) * np.ones(mu.shape)
    r     = rz[0].reshape(-1,1) * np.ones(mu.shape)
    vpar2 = vnorm * vnorm - 2 * bnorm * np.absolute(mu) / mass
    vpar  = np.sign(mu)*np.sqrt(vpar2)

    P  = mass * r * vpar + charge * psi

    if weights is None:
        weights = np.ones(muvals.shape)

    markers = np.zeros( (rhogrid.size-1, mugrid.size-1) )
    for i in range(rhogrid.size-1):
        for j in range(mugrid.size-1):
            a = np.array([ P[i,j], P[i+1,j] ])#, P[i,j+1], P[i+1,j+1] ])

            ids = np.logical_and.reduce( [
                Pvals >= np.nanmin(a), Pvals < np.nanmax(a),
                muvals >= mugrid[j], muvals < mugrid[j+1]
            ] )

            markers[i,j] = np.nansum(weights[ids])

    return markers


def checkrhomu2Pmu(fn):
    """
    Check rho_omp-mu & P-mu transformation.
    """
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    mass   = const.physical_constants["alpha particle mass"][0]
    energy = 3.5e6 * const.e
    charge = 2 * const.e

    rhomin = 0.8
    rhomax = 1.0

    nP  = 200
    nmu = 100
    nksi = 20
    nrho = 50
    try:
        fig  = plt.figure()
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

        axes = fig.add_subplot(2,2,1)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = mapPmu2rhomu(a5, mass, charge, energy,
                               Pvals.ravel(), muvals.ravel(),
                               rhogrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,2)
        h = axes.pcolormesh(rhogrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = np.ones(Pvals.shape)
        markers = mapPmu2rhomu(a5, mass, charge, energy,
                               Pvals.ravel(), muvals.ravel(),
                               rhogrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,3)
        h = axes.pcolormesh(rhogrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        rhovals = rhogrid[:-1] + (rhogrid[1] - rhogrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2

        rhovals, muvals = np.meshgrid(rhovals, muvals, indexing='ij')
        markers = maprhomu2Pmu(a5, mass, charge, energy,
                               rhovals.ravel(), muvals.ravel(),
                               Pgrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,4)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


def mapPmu2rzk(a5, mass, charge, energy, Pvals, muvals,
               rgrid, zgrid, ksigrid, weights=None):
    """
    Map P-mu points to real space.
    """

    vnorm = np.sqrt(2*energy/mass)
    if weights is None:
        weights = np.ones(muvals.shape)

    Pvals   = Pvals.ravel()
    muvals  = muvals.ravel()
    weights = weights.ravel()

    R, Z, Ksi = np.meshgrid(rgrid, zgrid, ksigrid, indexing='ij')
    psi   = a5.evaluate(R, 0, Z, 0, "psi").reshape(R.shape)
    bnorm = a5.evaluate(R, 0, Z, 0, "bnorm").reshape(R.shape)
    Mu = np.sign(Ksi)*mass * (1-Ksi*Ksi) * vnorm * vnorm / (2*bnorm)
    P  = mass * R * Ksi * vnorm  + charge * psi

    markers = np.zeros( (R.shape[0]-1,R.shape[1]-1,R.shape[2]-1) )
    newweight = np.zeros(weights.shape)
    for i in range(R.shape[0]-1):
        for j in range(R.shape[1]-1):
            for k in range(R.shape[2]-1):
                vals = [Mu[i,j,k], Mu[i+1,j,k], Mu[i,j+1,k], Mu[i,j,k],
                        Mu[i+1,j+1,k], Mu[i+1,j,k], Mu[i,j+1,k],
                        Mu[i+1,j+1,k]]

                a = np.logical_and.reduce([np.nanmin(vals) <= muvals,
                                           np.nanmax(vals) > muvals])

                vals = [P[i,j,k], P[i+1,j,k], P[i,j+1,k], P[i,j,k],
                        P[i+1,j+1,k], P[i+1,j,k], P[i,j+1,k],
                        P[i+1,j+1,k]]

                b = np.logical_and.reduce([np.nanmin(vals) <= Pvals,
                                           np.nanmax(vals) > Pvals])

                ids = np.logical_and.reduce([a, b])
                if np.sum(ids.ravel()) > 0:
                    newweight[ids] += 1

    weights = weights/newweight
    for i in range(R.shape[0]-1):
        for j in range(R.shape[1]-1):
            for k in range(R.shape[2]-1):
                vals = [Mu[i,j,k], Mu[i+1,j,k], Mu[i,j+1,k], Mu[i,j,k],
                        Mu[i+1,j+1,k], Mu[i+1,j,k], Mu[i,j+1,k],
                        Mu[i+1,j+1,k]]

                a = np.logical_and.reduce([np.nanmin(vals) <= muvals,
                                           np.nanmax(vals) > muvals])

                vals = [P[i,j,k], P[i+1,j,k], P[i,j+1,k], P[i,j,k],
                        P[i+1,j+1,k], P[i+1,j,k], P[i,j+1,k],
                        P[i+1,j+1,k]]

                b = np.logical_and.reduce([np.nanmin(vals) <= Pvals,
                                           np.nanmax(vals) > Pvals])

                ids = np.logical_and.reduce([a, b])
                if np.sum(ids.ravel()) > 0:
                    markers[i,j,k] = np.sum( weights[ids].ravel() )

    return markers


def maprzk2Pmu(a5, mass, charge, energy, rvals, zvals, ksivals,
               Pgrid, mugrid, weights=None):
    """
    Map real space points to P-mu grid.
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


def checkrzk2Pmu(fn):
    """
    Check real-space & P-mu transformation.
    """
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    mass   = const.physical_constants["alpha particle mass"][0]
    energy = 3.5e6 * const.e
    charge = 2 * const.e

    rhomin = 0.8
    rhomax = 1.0

    nP  = 200
    nmu = 100
    nksi = 20
    nrho = 50

    rz   = a5.get_rhotheta_rz( np.array([0.1, 1]), 0, 0, 0 )
    rmax = rz[0][-1]
    rz   = a5.get_rhotheta_rz( np.array([0.1, 1]), 180, 0, 0 )
    rmin = rz[0][-1]
    rz   = a5.get_rhotheta_rz( np.array([0.1, 1]), 90, 0, 0 )
    zmax = rz[1][-1]
    rz   = a5.get_rhotheta_rz( np.array([0.1, 1]), -90, 0, 0 )
    zmin = rz[1][-1]

    nr = 30
    nz = 60
    nk = 50
    rgrid = np.linspace(rmin, rmax, nr)
    zgrid = np.linspace(zmin, zmax, nz)
    kgrid = np.linspace(-1, 1, nk)

    try:
        fig  = plt.figure()
        Pgrid, mugrid, rhogrid, pitchgrid = initgrid(a5, mass, charge, energy,
                                                     rhomin, rhomax,
                                                     nP, nmu, nrho, nksi, 1e-3)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = np.ones(muvals.shape)
        markers = mapPmu2rzk(a5, mass, charge, energy,
                             Pvals.ravel(), muvals.ravel(),
                             rgrid, zgrid, kgrid, markers.ravel())

        axes = fig.add_subplot(2,3,1)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(markers, axis=2).transpose())
        plt.colorbar(h, ax=axes)
        kvals = kgrid[:-1] + (kgrid[1] - kgrid[0])/2
        axes = fig.add_subplot(2,3,2)
        h = axes.plot(kvals, np.sum(np.sum(markers, axis=0), axis=0))

        rvals = rgrid[:-1] + (rgrid[1] - rgrid[0]) / 2
        zvals = zgrid[:-1] + (zgrid[1] - zgrid[0]) / 2
        kvals = kgrid[:-1] + (kgrid[1] - kgrid[0]) / 2
        rvals, zvals, kvals = np.meshgrid(rvals, zvals, kvals, indexing='ij')
        markers = maprzk2Pmu(a5, mass, charge, energy,
                             rvals.ravel(), zvals.ravel(), kvals.ravel(),
                             Pgrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,3,3)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        rvals = rgrid[:-1] + (rgrid[1] - rgrid[0]) / 2
        zvals = zgrid[:-1] + (zgrid[1] - zgrid[0]) / 2
        kvals = kgrid[:-1] + (kgrid[1] - kgrid[0]) / 2
        rvals, zvals, kvals = np.meshgrid(rvals, zvals, kvals, indexing='ij')
        markers = np.ones(rvals.shape)
        markers = maprzk2Pmu(a5, mass, charge, energy,
                             rvals.ravel(), zvals.ravel(), kvals.ravel(),
                             Pgrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,3,4)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = mapPmu2rzk(a5, mass, charge, energy,
                             Pvals.ravel(), muvals.ravel(),
                             rgrid, zgrid, kgrid, markers.ravel())

        axes = fig.add_subplot(2,3,5)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(markers, axis=2).transpose())
        plt.colorbar(h, ax=axes)
        kvals = kgrid[:-1] + (kgrid[1] - kgrid[0])/2
        axes = fig.add_subplot(2,3,6)
        h = axes.plot(kvals, np.sum(np.sum(markers, axis=0), axis=0))

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


def maprhoksi2Pmu(a5, mass, charge, energy, rhovals, ksivals,
                 Pgrid, mugrid, weights=None):
    """
    Map rho_omp-ksi values to a P-mu grid.
    """
    rz    = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    vnorm = np.sqrt(2*energy/mass)

    muvals = np.sign(ksivals)*(1 - ksivals*ksivals) * mass * vnorm * vnorm / (2*bnorm)
    Pvals  = mass * rz[0] * vnorm * ksivals + charge * psi

    mumax = mass * vnorm * vnorm / (2*bnorm)

    if weights is None:
        weights = np.ones(rhovals.shape)

    markers = np.zeros( (Pgrid.size-1, mugrid.size-1) )
    for i in range(Pgrid.size-1):
        for j in range(mugrid.size-1):
            ids = np.logical_and.reduce( [
                Pvals >= Pgrid[i], Pvals < Pgrid[i+1],
                muvals >= mugrid[j], muvals < mugrid[j+1]
            ] )
            markers[i,j] = np.sum(weights[ids])

    return markers


def mapPmu2rhoksi(a5, mass, charge, energy, Pvals, muvals, rhogrid, ksigrid,
                 weights=None):
    """
    Map P-mu values to a rho_omp-ksi grid.
    """
    rz    = a5.get_rhotheta_rz( rhogrid, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    vnorm = np.sqrt(2*energy/mass)

    ksi   = ksigrid.transpose()
    psi   = psi.reshape(-1,1) * np.ones(ksi.shape)
    bnorm = bnorm.reshape(-1,1) * np.ones(ksi.shape)
    r     = rz[0].reshape(-1,1) * np.ones(ksi.shape)
    vpar  = ksi*vnorm

    P  = mass * r * vpar + charge * psi
    mu = np.sign(ksi+1e-8) * (1 - ksi*ksi) * mass * vnorm * vnorm / (2*bnorm)
    mumax = mass * vnorm * vnorm / (2*bnorm)

    if weights is None:
        weights = np.ones(muvals.shape)

    markers = np.zeros( (rhogrid.size-1, ksigrid.size-1) )
    for i in range(rhogrid.size-1):
        for j in range(ksigrid.size-1):
            a = np.array([ P[i,j], P[i+1,j] ])
            b = np.array([ mu[i,j], mu[i,j+1] ])

            ids = np.logical_and.reduce( [
                Pvals >= np.nanmin(a), Pvals < np.nanmax(a)
            ] )

            if ksigrid.size % 2 == 1 and j == (ksigrid.size-1)/2-1 :
                muids = np.logical_and.reduce( [ muvals <=  mu[i,j],
                                                 -muvals < mumax[i,j] ] )
            elif ksigrid.size % 2 == 0 and j == ksigrid.size/2-1:
                ids1 = np.logical_and.reduce( [ muvals <= mu[i,j],
                                                -muvals <= mumax[i,j] ] )
                ids2 = np.logical_and.reduce( [ muvals > mu[i,j+1],
                                                muvals < mumax[i,j] ] )
                muids = np.logical_or.reduce( [ids1, ids2] )
            else:
                muids = np.logical_and.reduce( [ muvals >= np.nanmin(b),
                                                 muvals < np.nanmax(b) ] )

            ids = np.logical_and.reduce( [ids, muids] )

            markers[i,j] = np.nansum(weights[ids])

    return markers


def checkrhoksi2Pmu(fn):
    """
    Check rho_omp-ksi & P-mu transformation.
    """
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    mass   = const.physical_constants["alpha particle mass"][0]
    energy = 3.5e6 * const.e
    charge = 2 * const.e

    rhomin = 0.8
    rhomax = 1.0

    nP  = 200
    nmu = 200
    nksi = 200
    nrho = 200
    try:
        fig  = plt.figure()
        Pgrid, mugrid, rhogrid, ksigrid = initgrid(a5, mass, charge, energy,
                                                   rhomin, rhomax,
                                                   nP, nmu, nrho, nksi, 1e-3)

        rhovals = rhogrid[:-1] + (rhogrid[1] - rhogrid[0]) / 2
        ksivals = ksigrid[:-1] + (ksigrid[1] - ksigrid[0]) / 2

        rhovals, ksivals = np.meshgrid(rhovals, ksivals, indexing='ij')
        markers = np.ones(rhovals.shape)
        markers = maprhoksi2Pmu(a5, mass, charge, energy,
                                rhovals.ravel(), ksivals.ravel(),
                                Pgrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,1)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = mapPmu2rhoksi(a5, mass, charge, energy,
                                Pvals.ravel(), muvals.ravel(),
                                rhogrid, ksigrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,2)
        h = axes.pcolormesh(rhogrid, ksigrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = np.ones(Pvals.shape)
        markers = mapPmu2rhoksi(a5, mass, charge, energy,
                                Pvals.ravel(), muvals.ravel(),
                                rhogrid, ksigrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,3)
        h = axes.pcolormesh(rhogrid, ksigrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        rhovals = rhogrid[:-1] + (rhogrid[1] - rhogrid[0]) / 2
        ksivals = ksigrid[:-1] + (ksigrid[1] - ksigrid[0]) / 2

        rhovals, ksivals = np.meshgrid(rhovals, ksivals, indexing='ij')
        markers = maprhoksi2Pmu(a5, mass, charge, energy,
                                rhovals.ravel(), ksivals.ravel(),
                                Pgrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,2,4)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


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
        markers = mapPmu2rhomu(a5, mass, charge, energy,
                               Pvals.ravel(), muvals.ravel(),
                               rhogrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,4,2)
        h = axes.pcolormesh(rhogrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)
        return

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = mapPmu2rzk(a5, mass, charge, energy,
                             Pvals.ravel(), muvals.ravel(),
                             rgrid, zgrid, ksigrid, markers.ravel())

        axes = fig.add_subplot(2,4,2)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(markers, axis=2).transpose())
        plt.colorbar(h, ax=axes)
        ksivals = ksigrid[:-1] + (ksigrid[1] - ksigrid[0])/2
        axes = fig.add_subplot(2,4,3)
        h = axes.plot(ksivals, np.sum(np.sum(markers, axis=0), axis=0))

        dist = markers
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
        r   = np.zeros((n,))
        z   = np.zeros((n,))
        ksi = np.zeros((n,))
        w   = np.zeros((n,))
        n0  = 0
        for i in range(dist.shape[0]):
            for j in range(dist.shape[1]):
                for k in range(dist.shape[2]):
                    n = int(dist[i,j,k])
                    if n == 0:
                        continue

                    r[n0:n0+n]   = rgrid[i]   + np.random.rand(n)*(
                        rgrid[1] - rgrid[0])
                    z[n0:n0+n]   = zgrid[j]   + np.random.rand(n)*(
                        zgrid[1] - zgrid[0])
                    ksi[n0:n0+n] = ksigrid[k] + np.random.rand(n)*(
                        ksigrid[1] - ksigrid[0])
                    w[n0:n0+n]   = weights[i,j,k] / n
                    n0 = n0+n

        n = r.size

        Nmrk = 20000
        ids = np.arange(n)
        np.random.shuffle(ids)
        r   = np.squeeze(r)[ids[:Nmrk]]
        z   = np.squeeze(z)[ids[:Nmrk]]
        ksi = np.squeeze(ksi)[ids[:Nmrk]]
        w   = np.squeeze(w * (Nmrk/n))[ids[:Nmrk]]

        markers = maprzk2Pmu(a5, mass, charge, energy, r, z, ksi,
                             Pgrid, mugrid, weights=None)

        axes = fig.add_subplot(2,4,5)
        h = axes.pcolormesh(Pgrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        Pvals   = Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2
        muvals  = mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2
        Pvals, muvals = np.meshgrid(Pvals, muvals, indexing='ij')
        markers = mapPmu2rhomu(a5, mass, charge, energy,
                               Pvals.ravel(), muvals.ravel(),
                               rhogrid, mugrid, weights=markers.ravel())

        axes = fig.add_subplot(2,4,6)
        h = axes.pcolormesh(rhogrid, mugrid, markers.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)
