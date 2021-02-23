"""
Transformations in phase-space and to/from real-space.

Phase-space transformations include:
- P,mu -> rho_omp, ksi_omp
- P,mu -> rho_omp, mu
and inverse operations. Here P is canonical toroidal momentum, mu is magnetic
moment, rho_omp is outer mid plane rho and ksi_omp is pitch at outer mid plane.
To clarify, rho_omp and ksi_omp are the values a marker with given P and mu
would have when crossing OMP. Since trapped particles have two crossing points,
the one where ksi_omp has the same sign as the current pitch is chosen. The sign
of pitch is transmitted as a sign of mu, i.e. mu = sign(ksi) * magnetic moment.

The real-space transformation is P,mu -> R,z,ksi and its inverse.

All parameters in this module are in SI units. Ascotpy is required and a file
with a magnetic field. Functions accept either pre-initialized Ascotpy object
or name of the HDF5 file, in which case magnetic field module is initialized
and freed inside the function.

This module is for tokamaks only.

File: transport/ppmappings.py
"""
import numpy as np
import scipy.constants as const
import importlib.util as util

from unyt import c as speed_of_light

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt


def evalpnorm(mass, energy):
    """
    Small helper function to evaluate pnorm.
    """
    gamma = 1 + energy / (mass * speed_of_light.v**2)
    return np.sqrt(gamma**2 - 1) * mass * speed_of_light.v


def evalPmu(a5, mass, charge, energy, r, z, ksi):
    """
    Evaluate P and mu at the given coordinates.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        r : array_like (n,1) <br>
            Major radius coordinates [m].
        z : array_like (n,1) <br>
            z coordinates [m].
        ksi : array_like (n,1) <br>
            Pitch coordinates (v_para/v_tot).

    Returns:
        Tuple (P,mu) where P and mu are (n,1) arrays ([kg*m/s], [J/T]).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    psi   = a5.evaluate(r, 0, z, 0, "psi")
    bphi  = a5.evaluate(r, 0, z, 0, "bphi")
    bnorm = a5.evaluate(r, 0, z, 0, "bnorm")
    pnorm = evalpnorm(mass, energy)

    mu = np.sign(ksi) * ( 1 - ksi**2 ) * pnorm**2 / ( 2 * bnorm * mass )
    P  = r * pnorm * ksi * bphi / bnorm + charge * psi

    if free:
        a5.free(bfield=True)

    return (P, mu)


def istrapped(a5, mass, charge, energy, P, mu, rmin=None):
    """
    Check if a marker with given coordinates is trapped.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        P : array_like (n,1) <br>
            Canonical toroidal momentum coordinates [kg*m/s].
        mu : array_like (n,1) <br>
            Magnetic moment coordinates [J/T].
        rlim : tuple or None <br>
            Major grid minimum value that markers may have. If None, ther value
            is picked at inner midplane separatrix.

    Returns:
        Boolean array with True elements for coordinates that are for trapped
        markers.
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    if rmin is None:
        rz    = a5.get_rhotheta_rz( np.array([0,1]), 180, 0, 0 )
        rmin = rz[0][-1]

    # Make R grid for evaluating bnorm and psi.
    nr = 100
    axis  = a5.evaluate(0, 0, 0, 0, "axis")
    rgrid = np.linspace(rmin, axis["axisr"], nr).ravel()

    bphi  = a5.evaluate(rgrid, 0, axis["axisz"], 0, "bphi").ravel()
    bnorm = a5.evaluate(rgrid, 0, axis["axisz"], 0, "bnorm").ravel()
    psi   = a5.evaluate(rgrid, 0, axis["axisz"], 0, "psi").ravel()
    pnorm = evalpnorm(mass, energy)

    # Find which (P,mu) coordinates correspond to trapped markers. A marker is
    # trapped if a) mu is imaginary for all r < raxis or b) P(r,mu) - P has a
    # root for rmin < r <raxis.
    trapped = np.zeros(mu.shape) == 1
    for i in range(mu.size):
        ppar2  = pnorm*pnorm - 2 * np.absolute(mu[i]) * bnorm * mass
        if all(ppar2 < 0):
            # mu is imaginary for all r < raxis. Therefore marker orbit has no
            # solution there which means the marker is trapped.
            trapped[i] = True
            continue

        # Evaluate P(r,mu)
        val = charge * psi + np.sign(mu[i]) * rgrid * np.sqrt(ppar2) \
              * bphi / bnorm
        val[ppar2 < 0] = np.nan

        if not (any(val - P[i] > 0) and any(val - P[i] < 0)):
            # No root for P(r,mu) - P, marker does not pass IMP and hence is
            # trapped
            trapped[i] = True

    if free:
        a5.free(bfield=True)

    return trapped


def initgridPmu(a5, mass, charge, energy, nP, nmu, padding):
    """
    Initialize (P,mu) grid that covers the whole phase-space.

    This function is useful as both P and mu depend not only on the marker
    parameters, but also on the magnetic background. Here values for P and mu
    are evaluated at OMP, and the extreme values (with some padding) are used
    to define the grid.

    mugrid is given in interval [-mumax, mumax].

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        nP : int <br>
            Number of P grid points.
        nmu : int <br>
            Number of mu grid points.
        padding : array_like (n,1) <br>
            How many percent grids are stretched beyond the extreme values.

    Returns:
        Tuple (Pgrid [kg*m/s], mugrid [J/T]).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    # OMP rz coordinates
    rz = a5.get_rhotheta_rz( np.array([0,1]), 0, 0, 0 )

    # mumax can be found assuming pperp = pnorm and using OMP separatrix value
    # for bmin
    bmin = a5.evaluate(rz[0][-1], 0, rz[1][-1], 0, "bnorm")
    pnorm = evalpnorm(mass, energy)

    mumax = pnorm**2 / ( 2 * bnorm * mass )
    mumin = -mumax

    mugrid = np.linspace(mumin*(1+padding), mumax*(1+padding), nmu)

    # Pmax and min are located at vpar = +-vnorm and psi = psimin or psimax.
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")

    if psi[0] < psi[1]:
        psimin = psi[0]
        rmin   = rz[0][0]
        psimax = psi[1]
        rmax   = rz[0][1]
    else:
        psimin = psi[1]
        rmin   = rz[0][1]
        psimax = psi[0]
        rmax   = rz[0][0]

    bpmin = a5.evaluate(rmin, 0, rz[1][-1], 0, "bphi") / \
            a5.evaluate(rmin, 0, rz[1][-1], 0, "bnorm")
    bpmax = a5.evaluate(rmax, 0, rz[1][-1], 0, "bphi") / \
            a5.evaluate(rmax, 0, rz[1][-1], 0, "bnorm")

    P  = np.array([ rmin * pnorm * bpmin + charge * psimin,
                   -rmin * pnorm * bpmin + charge * psimin,
                    rmax * pnorm * bpmax + charge * psimax,
                   -rmax * pnorm * bpmax + charge * psimax])

    Pmin = np.nanmin(P)
    Pmax = np.nanmax(P)

    if Pmin < 0:
        Pmin = Pmin*(1+padding)
    else:
        Pmin = Pmin*(1-padding)

    if Pmax < 0:
        Pmax = Pmax*(1-padding)
    else:
        Pmax = Pmax*(1+padding)

    Pgrid = np.linspace(Pmin, Pmax, nP)

    if free:
        a5.free(bfield=True)

    return (Pgrid, mugrid)


def maprhomu2Pmu(a5, mass, charge, energy, rhovals, muvals,
                 Pgrid, mugrid, weights=None):
    """
    Map (rho_omp,mu) coordinates to (P,mu) grid.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        rhovals : array_like (n,1) <br>
            Marker rho_omp values.
        muvals : array_like (n,1) <br>
            Marker mu values (including sign).
        Pgrid : array_like (n,1) <br>
            P coordinate grid.
        mugrid : array_like (n,1) <br>
            mu coordinate grid.
        weights : array_like (n,1), optional <br>
            Weights of markers. Uniform weights if None.

    Returns:
        Tuple (ordinate,Pvals,muvals) where ordinate is the resulting histogram
        in (Pgrid,mugrid) grid and Pvals, and muvals are the exact (P,mu) values
        corresponding to (rhovals,muvals).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    rz    = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    bphi  = a5.evaluate(rz[0], 0, rz[1], 0, "bphi")
    pnorm = evalpnorm(mass, energy)

    ppar   = np.sign(muvals) * np.sqrt(pnorm*pnorm
                                       - 2*bnorm*np.absolute(muvals)*mass)
    Pvals  = rz[0] * ppar *bphi / bnorm + charge * psi

    if weights is None:
        weights = np.ones(muvals.shape)

    ordinate = np.histogram2d(Pvals, muvals,
                              bins=[Pgrid.ravel(), mugrid.ravel()],
                              weights=weights)[0]

    if free:
        a5.free(bfield=True)

    return (ordinate, Pvals, muvals)


def mapPmu2rhomu(a5, mass, charge, energy, Pvals, muvals, rhogrid, mugrid,
                 weights=None):
    """
    Map (P,mu) coordinates to (rho_omp,mu) grid.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        Pvals : array_like (n,1) <br>
            Marker P values.
        muvals : array_like (n,1) <br>
            Marker mu values (including sign).
        rhogrid : array_like (n,1) <br>
            rho_omp coordinate grid.
        mugrid : array_like (n,1) <br>
            mu coordinate grid.
        weights : array_like (n,1), optional <br>
            Weights of markers. Uniform weights if None.

    Returns:
        Tuple (ordinate,rhovals,muvals) where ordinate is the resulting
        histogram in (rhogrid,mugrid) grid and rhovals and muvals are values
        each Pvals,muvals pair correspond to (center value of the bin they
        belong to).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    rz    = a5.get_rhotheta_rz( rhogrid, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    bphi  = a5.evaluate(rz[0], 0, rz[1], 0, "bphi")
    pnorm = evalpnorm(mass, energy)

    mu    = mugrid.transpose()
    psi   = psi.reshape(-1,1) * np.ones(mu.shape)
    bnorm = bnorm.reshape(-1,1) * np.ones(mu.shape)
    bphi  = bphi.reshape(-1,1) * np.ones(mu.shape)
    r     = rz[0].reshape(-1,1) * np.ones(mu.shape)
    ppar2 = pnorm * pnorm - 2 * bnorm * np.absolute(mu) * mass
    ppar  = np.sign(mu)*np.sqrt(ppar2)

    P  = r * ppar * bphi / bnorm + charge * psi

    if weights is None:
        weights = np.ones(muvals.shape)

    rhovals   = np.ones(Pvals.shape)
    newmuvals = np.ones(Pvals.shape)
    ordinate  = np.zeros( (rhogrid.size-1, mugrid.size-1) )
    for i in range(rhogrid.size-1):
        for j in range(mugrid.size-1):
            a = np.array([ P[i,j], P[i+1,j] ])

            ids = np.logical_and.reduce( [
                Pvals >= np.nanmin(a), Pvals < np.nanmax(a),
                muvals >= mugrid[j], muvals < mugrid[j+1]
            ] )

            rhovals[ids]   = (rhogrid[i] + rhogrid[i+1]) / 2
            newmuvals[ids] = (mugrid[j] + mugrid[j+1]) / 2
            ordinate[i,j]   = np.nansum(weights[ids])

    if free:
        a5.free(bfield=True)

    return (ordinate, rhovals, newmuvals)


def test_rhomu2Pmu(fn):
    """
    Verify rho_omp-mu & P-mu transformation.

    This function populates (rho,mu) mesh uniformly, transforms it to (P,mu)
    mesh, and transforms the resulting values back to (rho,mu) grid. Therefore
    the initial and final histograms should resemble each other. Process is
    repeated starting from uniform (P,mu) population. Final and intermediate
    histograms are plotted.

    Args:
        fn : str <br>
            Filename for ascot HDF5 file which must have magnetic field present.
    """
    from a5py.ascot5io.ascot5 import Ascot
    from a5py.ascotpy.ascotpy import Ascotpy
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    # Perform check for alpha particles
    mass   = const.physical_constants["alpha particle mass"][0]
    energy = 3.5e6 * const.e
    charge = 2 * const.e

    rhogrid = np.linspace( 0, 1, 100)

    nP  = 100
    nmu = 100
    try:
        fig  = plt.figure()
        Pgrid, mugrid = initgridPmu(a5, mass, charge, energy,
                                    nP, nmu, padding=1e-3)

        # Convert mesh midpoints to value pairs, transform, and plot
        rhovals, muvals = np.meshgrid(
            rhogrid[:-1] + (rhogrid[1] - rhogrid[0]) / 2,
            mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2,
            indexing='ij')

        weights = np.ones(rhovals.shape)
        weights, Pvals, muvals = maprhomu2Pmu(a5, mass, charge, energy,
                                              rhovals.ravel(), muvals.ravel(),
                                              Pgrid, mugrid,
                                              weights=weights.ravel())

        axes = fig.add_subplot(2,2,1)
        h = axes.pcolormesh(Pgrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        # Repeat
        weights = np.ones(Pvals.shape)
        weights = mapPmu2rhomu(a5, mass, charge, energy,
                               Pvals.ravel(), muvals.ravel(),
                               rhogrid, mugrid, weights=weights.ravel())[0]

        axes = fig.add_subplot(2,2,2)
        h = axes.pcolormesh(rhogrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        # Now do the whole stuff otherway around
        Pvals, muvals = np.meshgrid(
            Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2,
            mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2,
            indexing='ij')

        weights = np.ones(Pvals.shape)
        weights, rhovals, muvals = mapPmu2rhomu(a5, mass, charge, energy,
                                                Pvals.ravel(), muvals.ravel(),
                                                rhogrid, mugrid,
                                                weights=weights.ravel())

        axes = fig.add_subplot(2,2,3)
        h = axes.pcolormesh(rhogrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)


        # And repeat
        weights = maprhomu2Pmu(a5, mass, charge, energy,
                               rhovals.ravel(), muvals.ravel(),
                               Pgrid, mugrid, weights=weights.ravel())[0]

        axes = fig.add_subplot(2,2,4)
        h = axes.pcolormesh(Pgrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


def maprhoksi2Pmu(a5, mass, charge, energy, rhovals, ksivals,
                  Pgrid, mugrid, weights=None):
    """
    Map (rho_omp,rho_ksi) coordinates to (P,mu) grid.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        rhovals : array_like (n,1) <br>
            Marker rho_omp values.
        ksivals : array_like (n,1) <br>
            Marker pitch values.
        Pgrid : array_like (n,1) <br>
            P coordinate grid.
        mugrid : array_like (n,1) <br>
            mu coordinate grid.
        weights : array_like (n,1), optional <br>
            Weights of markers. Uniform weights if None.

    Returns:
        Tuple (ordinate,Pvals,muvals) where ordinate is the resulting histogram
        in (Pgrid,mugrid) grid and Pvals, and muvals are the exact (P,mu) values
        corresponding to (rhovals,ksivals).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    ksivals = ksivals
    weights = weights
    rhovals = rhovals

    rz    = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    bphi  = a5.evaluate(rz[0], 0, rz[1], 0, "bphi")
    pnorm = evalpnorm(mass, energy)

    muvals = (1 - 2*(ksivals<0)) * (1 - ksivals*ksivals) \
             * pnorm * pnorm / ( 2 * bnorm * mass )
    Pvals  = rz[0] * pnorm * ksivals * bphi / bnorm + charge * psi

    if weights is None:
        weights = np.ones(rhovals.shape)

    ordinate = np.histogram2d(Pvals,
                              muvals,
                              bins=[Pgrid.ravel(), mugrid.ravel()],
                              weights=weights)[0]

    if free:
        a5.free(bfield=True)

    return (ordinate, Pvals, muvals)


def mapPmu2rhoksi(a5, mass, charge, energy, Pvals, muvals, rhogrid, ksigrid,
                  weights=None):
    """
    Map (P,mu) coordinates to (rho_omp,ksi_omp) grid.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        Pvals : array_like (n,1) <br>
            Marker P values.
        muvals : array_like (n,1) <br>
            Marker mu values (including sign).
        rhogrid : array_like (n,1) <br>
            rho_omp coordinate grid.
        ksigrid : array_like (n,1) <br>
            Pitch coordinate grid.
        weights : array_like (n,1), optional <br>
            Weights of markers. Uniform weights if None.

    Returns:
        Tuple (ordinate,rhovals,ksivals) where ordinate is the resulting
        histogram in (rhogrid,ksigrid) grid and rhovals and ksivals are values
        each Pvals,muvals pair correspond to (center value of the bin they
        belong to).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    rz    = a5.get_rhotheta_rz( rhogrid, 0, 0, 0 )
    psi   = a5.evaluate(rz[0], 0, rz[1], 0, "psi")
    bnorm = a5.evaluate(rz[0], 0, rz[1], 0, "bnorm")
    bphi  = a5.evaluate(rz[0], 0, rz[1], 0, "bphi")
    pnorm = evalpnorm(mass, energy)

    ksi   = ksigrid.transpose()
    psi   = psi.reshape(-1,1) * np.ones(ksi.shape)
    bnorm = bnorm.reshape(-1,1) * np.ones(ksi.shape)
    bphi  = bphi.reshape(-1,1) * np.ones(ksi.shape)
    r     = rz[0].reshape(-1,1) * np.ones(ksi.shape)
    ksi   = np.ones(rz[0].reshape(-1,1).shape).reshape(-1,1)*ksigrid.transpose()
    ppar  = ksi*pnorm

    P  = r * ppar * bphi / bnorm + charge * psi
    mu = (1 - 2*(ksi < 0)) * (1 - ksi*ksi) * pnorm * pnorm / ( 2 * bnorm * mass )
    mumax = pnorm * pnorm / ( 2 * bnorm * mass )

    if weights is None:
        weights = np.ones(muvals.shape)

    rhovals = np.ones(Pvals.shape) + np.nan
    ksivals = np.ones(Pvals.shape) + np.nan
    ordinate = np.zeros( (rhogrid.size-1, ksigrid.size-1) )
    for i in range(rhogrid.size-1):
        for j in range(ksigrid.size-1):
            a = np.array([ P[i,j], P[i+1,j], P[i,j+1], P[i+1,j+1] ])
            b = np.array([ mu[i,j], mu[i,j+1], mu[i,j+1], mu[i+1,j+1]])
            mumaxval = np.nanmax([mumax[i,j], mumax[i+1,j]])

            ids = np.logical_and.reduce( [
                Pvals >= np.nanmin(a), Pvals < np.nanmax(a)
            ] )

            if ksigrid.size % 2 == 0 and j == ksigrid.size/2 - 1:
                a1 = np.logical_and.reduce( [
                    muvals < np.nanmax([mu[i,j], mu[i+1,j]]),
                    muvals >= -mumaxval ] )
                a2 = np.logical_and.reduce( [
                    muvals >= np.nanmin([mu[i,j+1], mu[i+1,j+1]]),
                    muvals <=  mumaxval ] )
                muids = np.logical_or.reduce([a1, a2])
            elif ksigrid.size % 2 == 1 and j == (ksigrid.size-1)/2-1:
                muids = np.logical_and.reduce(
                    [ muvals <  np.nanmax([mu[i,j], mu[i+1,j]]),
                      muvals >= -mumaxval] )
            elif ksigrid.size % 2 == 1 and j == (ksigrid.size-1)/2:
                muids = np.logical_and.reduce(
                    [ muvals >= np.nanmin([mu[i,j+1], mu[i+1,j+1]]),
                      muvals <= mumaxval] )
            else:
                muids = np.logical_and.reduce( [ muvals >= np.nanmin(b),
                                                 muvals < np.nanmax(b) ] )

            ids = np.logical_and.reduce( [ids, muids] )

            rhovals[ids] = (rhogrid[i] + rhogrid[i+1]) / 2
            ksivals[ids] = (ksigrid[j] + ksigrid[j+1]) / 2
            ordinate[i,j] = np.nansum(weights[ids])


    if free:
        a5.free(bfield=True)

    return (ordinate, rhovals, ksivals)


def test_rhoksi2Pmu(fn):
    """
    Verify rho_omp-ksi_omp & P-mu transformation.

    This function populates (rho,ksi) mesh uniformly, transforms it to (P,mu)
    mesh, and transforms the resulting values back to (rho,ksi) grid. Therefore
    the initial and final histograms should resemble each other. Process is
    repeated starting from uniform (P,mu) population. Final and intermediate
    histograms are plotted.

    Args:
        fn : str <br>
            Filename for ascot HDF5 file which must have magnetic field present.
    """
    from a5py.ascot5io.ascot5 import Ascot
    from a5py.ascotpy.ascotpy import Ascotpy
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    mass   = const.physical_constants["alpha particle mass"][0]
    energy = 3.5e6 * const.e
    charge = 2 * const.e

    rhogrid = np.linspace( 0, 1, 200)
    ksigrid = np.linspace(-1, 1, 50)

    nP  = 200
    nmu = 200
    try:
        fig  = plt.figure()
        Pgrid, mugrid = initgridPmu(a5, mass, charge, energy,
                                    nP, nmu, padding=1e-3)

        # Convert mesh midpoints to value pairs, transform, and plot
        rhovals, ksivals = np.meshgrid(
            rhogrid[:-1] + (rhogrid[1] - rhogrid[0]) / 2,
            ksigrid[:-1] + (ksigrid[1] - ksigrid[0]) / 2,
            indexing='ij')
        weights = np.ones(rhovals.shape)
        weights, Pvals, muvals = maprhoksi2Pmu(a5, mass, charge, energy,
                                               rhovals.ravel(), ksivals.ravel(),
                                               Pgrid, mugrid,
                                               weights=weights.ravel())

        axes = fig.add_subplot(2,2,1)
        h = axes.pcolormesh(Pgrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        # Convert back to rho_omp, ksi_omp space and plot
        weights = np.ones(Pvals.shape)
        weights, rhovals, ksivals = mapPmu2rhoksi(a5, mass, charge, energy,
                                Pvals.ravel(), muvals.ravel(),
                                rhogrid, ksigrid, weights=weights.ravel())

        axes = fig.add_subplot(2,2,2)
        h = axes.pcolormesh(rhogrid, ksigrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        # Repeat the process for a uniform (P,mu) population
        Pvals, muvals = np.meshgrid(
            Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2,
            mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2,
            indexing='ij')
        weights = np.ones(Pvals.shape)
        weights, rhovals, ksivals = mapPmu2rhoksi(a5, mass, charge, energy,
                                                  Pvals.ravel(), muvals.ravel(),
                                                  rhogrid, ksigrid,
                                                  weights=weights.ravel())

        axes = fig.add_subplot(2,2,3)
        h = axes.pcolormesh(rhogrid, ksigrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        # And transform back again and plot
        rhovals = rhovals[np.isnan(rhovals) == False]
        ksivals = ksivals[np.isnan(ksivals) == False]
        weights = np.ones(rhovals.shape)
        weights = maprhoksi2Pmu(a5, mass, charge, energy,
                                rhovals.ravel(), ksivals.ravel(),
                                Pgrid, mugrid, weights=weights.ravel())[0]

        axes = fig.add_subplot(2,2,4)
        h = axes.pcolormesh(Pgrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


def maprzk2Pmu(a5, mass, charge, energy, rvals, zvals, ksivals,
               Pgrid, mugrid, weights=None):
    """
    Map (rho_omp,mu) coordinates to (P,mu) grid.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        rvals : array_like (n,1) <br>
            Marker R values.
        zvals : array_like (n,1) <br>
            Marker z values.
        ksivals : array_like (n,1) <br>
            Marker pitch values.
        Pgrid : array_like (n,1) <br>
            P coordinate grid.
        mugrid : array_like (n,1) <br>
            mu coordinate grid.
        weights : array_like (n,1), optional <br>
            Weights of markers. Uniform weights if None.

    Returns:
        Tuple (ordinate,Pvals,muvals) where ordinate is the resulting histogram
        in (Pgrid,mugrid) grid and Pvals, and muvals are the exact (P,mu) values
        corresponding to (rvals,zvals,ksivals).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    bnorm = a5.evaluate(rvals, 0, zvals, 0, "bnorm")
    bphi  = a5.evaluate(rvals, 0, zvals, 0, "bphi")
    psi   = a5.evaluate(rvals, 0, zvals, 0, "psi")
    pnorm = evalpnorm(mass, energy)

    muvals = (1-2*(ksivals<0)) * (1-ksivals*ksivals) \
             * pnorm * pnorm / ( 2 * bnorm * mass )
    Pvals  = rvals * ksivals * pnorm * bphi / bnorm + charge * psi

    if weights is None:
        weights = np.ones(muvals.shape)

    ordinate = np.histogram2d(Pvals, muvals,
                              bins=[Pgrid.ravel(), mugrid.ravel()],
                              weights=weights)[0]

    if free:
        a5.free(bfield=True)

    return (ordinate, Pvals, muvals)


def mapPmu2rzk(a5, mass, charge, energy, Pvals, muvals,
               rgrid, zgrid, ksigrid, weights=None):
    """
    Map (P,mu) coordinates to (R,z,ksi) grid.

    This mapping is not one-to-one as single (P,mu) value corresponds to a whole
    trajectory in (R,z,ksi) real-space. Therefore, this function finds all
    (R,z,ksi) cells a (P,mu) correspond to, and divides the weight evenly among
    them.

    Args:
        a5 : Ascotpy, str <br>
            Ascotpy object or HDF5 filename.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        Pvals : array_like (n,1) <br>
            Marker P values.
        muvals : array_like (n,1) <br>
            Marker mu values (including sign).
        rgrid : array_like (n,1) <br>
            R coordinate grid.
        zgrid : array_like (n,1) <br>
            z coordinate grid.
        ksigrid : array_like (n,1) <br>
            ksi coordinate grid.
        weights : array_like (n,1), optional <br>
            Weights of markers. Uniform weights if None.

    Returns:
        Ordinate. (R,z,ksi) values are not returned as they are not unique.
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    pnorm = evalpnorm(mass, energy)
    if weights is None:
        weights = np.ones(muvals.shape)

    Pvals   = Pvals.ravel()
    muvals  = muvals.ravel()
    weights = weights.ravel()

    # Evaluate what P and mu are at each Rzk-grid node
    R, Z, Ksi = np.meshgrid(rgrid, zgrid, ksigrid, indexing='ij')
    psi   = a5.evaluate(R, 0, Z, 0, "psi").reshape(R.shape)
    bnorm = a5.evaluate(R, 0, Z, 0, "bnorm").reshape(R.shape)
    bphi  = a5.evaluate(R, 0, Z, 0, "bphi").reshape(R.shape)
    Mu = (1-2*(Ksi<0)) * (1-Ksi*Ksi) * pnorm * pnorm / ( 2 * bnorm * mass )
    P  = R * Ksi * pnorm * bphi / bnorm  + charge * psi
    #Pmid = charge * psi

    # Maximum possible mu value at each node
    Mumax = pnorm * pnorm / ( 2 * bnorm * mass )

    # The index search algorithm, which finds 8R,z,ksi) cells (P,mu) correspond
    # to, is run twice: first to find how many cells each (P,mu) correspond to
    # so that their weights can be dividec accordingly, and then to distribute
    # the weight among the cells.
    ncells  = np.zeros(weights.shape)
    for i in range(rgrid.size-1):
        for j in range(zgrid.size-1):
            for k in range(ksigrid.size-1):
                # Check which (Pvals, muvals) values are within this cell by
                # finding those for which Mu_a < mu < Mu_b and P_c < P < P_d
                # where a, b, c, and d are some grid nodes.

                # However, there is a caveat. When ksi passes zero, mu is
                # discontinuous. Therefore we separate those ksi grid points
                # (k and k+1) so we can use Mu_a < mu < Mumax and
                # -Mumax < mu < Mu_b where ksi(k) < 0 and ksi(k+1) > 0.
                vals0 = [Mu[i,j,k], Mu[i+1,j,k], Mu[i,j+1,k], Mu[i+1,j+1,k]]
                vals1 = [Mu[i,j,k+1], Mu[i+1,j,k+1], Mu[i,j+1,k+1],
                         Mu[i+1,j+1,k+1]]

                vals2 = [P[i,j,k], P[i+1,j,k], P[i,j+1,k], P[i+1,j+1,k]]
                vals3 = [P[i,j,k+1], P[i+1,j,k+1], P[i,j+1,k+1], P[i+1,j+1,k+1]]

                Mumaxval = np.nanmax([Mumax[i,j,k], Mumax[i+1,j,k],
                                      Mumax[i,j+1,k], Mumax[i+1,j+1,k]])

                if ksigrid.size % 2 == 0 and k == ksigrid.size/2-1:
                    a1 = np.logical_and.reduce( [ muvals <= np.nanmin(vals0),
                                                  muvals >  -Mumaxval ] )
                    a2 = np.logical_and.reduce( [ muvals <= Mumaxval,
                                                  muvals > np.nanmin(vals1) ] )
                    a = np.logical_or.reduce([a1, a2])

                elif ksigrid.size % 2 == 1 and k == (ksigrid.size-1)/2-1 :
                    a = np.logical_and.reduce( [ muvals < np.nanmax(vals0),
                                                 muvals >= -Mumaxval ] )

                elif ksigrid.size % 2 == 1 and k == (ksigrid.size-1)/2 :
                    a = np.logical_and.reduce( [ muvals >= np.nanmin(vals1),
                                                 muvals <= Mumaxval ] )
                else:
                    a = np.logical_and.reduce(
                        [ muvals >= np.nanmin(np.append(vals0,vals1)),
                          muvals < np.nanmax(np.append(vals0,vals1)) ] )

                b = np.logical_and.reduce(
                    [Pvals >= np.nanmin(np.append(vals2,vals3)),
                     Pvals < np.nanmax(np.append(vals2,vals3)) ])

                # Add this cell to all (P,mu) values that it holds
                ids = np.logical_and.reduce([a, b])
                if np.sum(ids.ravel()) > 0:
                    ncells[ids] += 1

    # Update the weights and run the algorithm again
    weights = weights/ncells
    ordinate = np.zeros( (rgrid.size-1,zgrid.size-1,ksigrid.size-1) )
    for i in range(rgrid.size-1):
        for j in range(zgrid.size-1):
            for k in range(ksigrid.size-1):
                vals0 = [Mu[i,j,k], Mu[i+1,j,k], Mu[i,j+1,k], Mu[i+1,j+1,k]]
                vals1 = [Mu[i,j,k+1], Mu[i+1,j,k+1], Mu[i,j+1,k+1],
                         Mu[i+1,j+1,k+1]]

                vals2 = [P[i,j,k], P[i+1,j,k], P[i,j+1,k], P[i+1,j+1,k]]
                vals3 = [P[i,j,k+1], P[i+1,j,k+1], P[i,j+1,k+1], P[i+1,j+1,k+1]]

                Mumaxval = np.nanmax([Mumax[i,j,k], Mumax[i+1,j,k],
                                      Mumax[i,j+1,k], Mumax[i+1,j+1,k]])

                #if ksigrid.size % 2 == 1 and k == (ksigrid.size-1)/2-1 :
                if ksigrid.size % 2 == 0 and k == ksigrid.size/2-1:
                    a1 = np.logical_and.reduce( [ muvals <= np.nanmin(vals0),
                                                  muvals >=-Mumaxval ] )
                    a2 = np.logical_and.reduce( [ muvals <= Mumaxval,
                                                  muvals >  np.nanmin(vals1) ] )
                    a = np.logical_or.reduce([a1, a2])

                #elif ksigrid.size % 2 == 0 and k == ksigrid.size/2-1:
                elif ksigrid.size % 2 == 1 and k == (ksigrid.size-1)/2-1 :
                    a = np.logical_and.reduce( [ muvals <= np.nanmax(vals0),
                                                 muvals >= -Mumaxval ] )
                #elif ksigrid.size % 2 == 0 and k == ksigrid.size/2:
                elif ksigrid.size % 2 == 1 and k == (ksigrid.size-1)/2 :
                    a = np.logical_and.reduce( [ muvals > np.nanmin(vals1),
                                                 muvals <= Mumaxval ] )
                else:
                    a = np.logical_and.reduce(
                        [ muvals >= np.nanmin(np.append(vals0,vals1)),
                          muvals < np.nanmax(np.append(vals0,vals1)) ] )

                b = np.logical_and.reduce(
                    [ Pvals >= np.nanmin(np.append(vals2,vals3)),
                      Pvals < np.nanmax(np.append(vals2,vals3)) ])

                # Weight of this cell is equal to sum of weights of those (P,mu)
                # values that this cell holds.
                ids = np.logical_and.reduce([a, b])
                if np.sum(ids.ravel()) > 0:
                    ordinate[i,j,k] = np.sum( weights[ids].ravel() )

    if free:
        a5.free(bfield=True)

    return ordinate


def test_rzk2Pmu(fn):
    """
    Verify R-z-ksi & P-mu transformation.

    This function populates (r,z,ksi) mesh uniformly, transforms it to (P,mu)
    mesh, and transforms the resulting values back to (r,z,ksi) grid. Therefore
    the initial and final histograms should resemble each other. Process is
    repeated starting from uniform (P,mu) population. Final and intermediate
    histograms are plotted.

    Args:
        fn : str <br>
            Filename for ascot HDF5 file which must have magnetic field present.
    """
    from a5py.ascot5io.ascot5 import Ascot
    from a5py.ascotpy.ascotpy import Ascotpy
    h5 = Ascot(fn)
    a5 = Ascotpy(fn)
    a5.init(bfield=True)

    mass   = const.physical_constants["alpha particle mass"][0]
    energy = 3.5e6 * const.e
    charge = 2 * const.e

    # Begin by finding min/max values for r and z based on separatrix.
    rho   = np.array([0,1])
    theta = np.array([1])
    rz = a5.get_rhotheta_rz(rho, theta, 0, 0 )
    rmin = rz[0][0]
    rmax = rz[0][0]
    zmin = rz[1][0]
    zmax = rz[1][0]
    for i in range(360):
        rz = a5.get_rhotheta_rz(rho, i*theta, 0, 0 )
        rmin = np.minimum(rmin, rz[0][1])
        rmax = np.maximum(rmax, rz[0][1])
        zmin = np.minimum(zmin, rz[1][1])
        zmax = np.maximum(zmax, rz[1][1])

    rgrid = np.linspace(rmin, rmax, 20)
    zgrid = np.linspace(zmin, zmax, 40)
    kgrid = np.linspace(-1, 1, 51)

    nP  = 100
    nmu = 100
    try:
        fig  = plt.figure()

        Pgrid, mugrid, = initgridPmu(a5, mass, charge, energy,
                                     nP, nmu, padding=1e-3)

        # Make uniform population on (R,z,ksi) grid.
        rvals, zvals, kvals = np.meshgrid(
            rgrid[:-1] + (rgrid[1] - rgrid[0]) / 2,
            zgrid[:-1] + (zgrid[1] - zgrid[0]) / 2,
            kgrid[:-1] + (kgrid[1] - kgrid[0]) / 2, indexing='ij')
        weights = np.ones(rvals.shape)

        # Transform it to (P,mu) space and plot
        weights, Pvals, muvals = maprzk2Pmu(a5, mass, charge, energy,
                                            rvals.ravel(), zvals.ravel(),
                                            kvals.ravel(), Pgrid, mugrid,
                                            weights=weights.ravel())

        axes = fig.add_subplot(2,3,1)
        h = axes.pcolormesh(Pgrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        # Now transform the population back to (R,z,ksi) grid and plot
        weights = np.ones(Pvals.shape)
        weights = mapPmu2rzk(a5, mass, charge, energy,
                             Pvals.ravel(), muvals.ravel(),
                             rgrid, zgrid, kgrid, weights=weights.ravel())

        # Rz
        axes = fig.add_subplot(2,3,2)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(weights, axis=2).transpose())
        plt.colorbar(h, ax=axes)

        # ksi
        axes = fig.add_subplot(2,3,3)
        h = axes.plot(kgrid[:-1] + (kgrid[1] - kgrid[0])/2,
                      np.sum(np.sum(weights, axis=0), axis=0))

        # Repeat the process other way around
        Pvals, muvals = np.meshgrid(
            Pgrid[:-1] + (Pgrid[1] - Pgrid[0]) / 2,
            mugrid[:-1] + (mugrid[1] - mugrid[0]) / 2,
            indexing='ij')
        weights = np.ones(Pvals.shape)
        weights = mapPmu2rzk(a5, mass, charge, energy,
                             Pvals.ravel(), muvals.ravel(),
                             rgrid, zgrid, kgrid, weights=weights.ravel())

        # Rz
        axes = fig.add_subplot(2,3,4)
        h = axes.pcolormesh(rgrid, zgrid, np.sum(weights, axis=2).transpose())
        plt.colorbar(h, ax=axes)

        # ksi
        axes = fig.add_subplot(2,3,5)
        h = axes.plot(kgrid[:-1] + (kgrid[1] - kgrid[0])/2,
                      np.sum(np.sum(weights, axis=0), axis=0))

        # Transform back to (P,mu) space and plot
        rvals, zvals, kvals = np.meshgrid(
            rgrid[:-1] + (rgrid[1] - rgrid[0]) / 2,
            zgrid[:-1] + (zgrid[1] - zgrid[0]) / 2,
            kgrid[:-1] + (kgrid[1] - kgrid[0]) / 2, indexing='ij')
        weights = maprzk2Pmu(a5, mass, charge, energy,
                             rvals.ravel(), zvals.ravel(), kvals.ravel(),
                             Pgrid, mugrid, weights=weights.ravel())[0]

        axes = fig.add_subplot(2,3,6)
        h = axes.pcolormesh(Pgrid, mugrid, weights.transpose())
        plt.colorbar(h, ax=axes)

        plt.show(block=False)

    except Exception as e:
        a5.free(bfield=True)
        raise e

    a5.free(bfield=True)


def maprhoksi2rzk(a5, mass, charge, energy, rhovals, ksivals, rgrid, zgrid,
                  ksigrid, weights=None):
    """
    Convert (rho_omp,ksi_omp) to (r,z,ksi) via (P,mu).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    Pgrid  = np.linspace(-1,1,5)
    mugrid = np.linspace(-1,1,5)
    _, Pvals, muvals = maprhoksi2Pmu(a5, mass, charge, energy, rhovals, ksivals,
                                     Pgrid, mugrid, weights=None)

    ordinate = mapPmu2rzk(a5, mass, charge, energy, Pvals, muvals,
                          rgrid, zgrid, ksigrid, weights=weights)

    if free:
        a5.free(bfield=True)

    return ordinate


def maprzk2rhoksi(a5, mass, charge, energy, rvals, zvals, ksivals, rhogrid,
                  ksigrid, weights=None):
    """
    Convert (r,z,ksi) to (rho_omp,ksi_omp) via (P,mu).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    Pgrid  = np.linspace(-1,1,5)
    mugrid = np.linspace(-1,1,5)
    #_, Pvals, muvals = maprzk2Pmu(a5, mass, charge, energy, rvals, zvals,
    #                              ksivals, Pgrid, mugrid, weights=None)
    Pvals, muvals = evalPmu(a5, mass, charge, energy, rvals, zvals, ksivals)

    ordinate, rhovals, ksivals = mapPmu2rhoksi(a5, mass, charge, energy, Pvals,
                                               muvals, rhogrid, ksigrid,
                                               weights=weights)

    if free:
        a5.free(bfield=True)

    return (ordinate, rhovals, ksivals)


def maprhoksi2rhomu(a5, mass, charge, energy, rhovals, ksivals, rhogrid,
                    mugrid, weights=None):
    """
    Convert (rho_omp,ksi_omp) to (rho,mu) via (P,mu).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    Pgrid  = np.linspace(-1,1,5)
    _, Pvals, muvals = maprhoksi2Pmu(a5, mass, charge, energy, rhovals, ksivals,
                                     Pgrid, mugrid, weights=None)

    ordinate, rhovals, muvals = mapPmu2rhomu(a5, mass, charge, energy, Pvals,
                                             muvals, rhogrid, mugrid,
                                             weights=weights)

    if free:
        a5.free(bfield=True)

    return (ordinate, rhovals, muvals)


def maprhomu2rhoksi(a5, mass, charge, energy, rhovals, muvals, rhogrid,
                    ksigrid, weights=None):
    """
    Convert (rho_omp,ksi_omp) to (rho,mu) via (P,mu).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    Pgrid  = np.linspace(-1,1,5)
    mugrid = np.linspace(-1,1,5)

    _, Pvals, muvals = maprhomu2Pmu(a5, mass, charge, energy, rhovals, muvals,
                                    Pgrid, mugrid, weights=None)

    ordinate, rhovals, ksivals = mapPmu2rhoksi(a5, mass, charge, energy, Pvals,
                                               muvals, rhogrid, ksigrid,
                                               weights=weights)

    if free:
        a5.free(bfield=True)

    return (ordinate, rhovals, ksivals)


def maprhomu2rzk(a5, mass, charge, energy, rhovals, muvals, rgrid, zgrid,
                 ksigrid, weights=None):
    """
    Convert (rho_omp,mu) to (r,z,ksi) via (P,mu).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    Pgrid  = np.linspace(-1,1,5)
    mugrid = np.linspace(-1,1,5)
    _, Pvals, muvals = maprhomu2Pmu(a5, mass, charge, energy, rhovals, muvals,
                                    Pgrid, mugrid, weights=None)

    ordinate = mapPmu2rzk(a5, mass, charge, energy, Pvals, muvals,
                          rgrid, zgrid, ksigrid, weights=weights)

    if free:
        a5.free(bfield=True)

    return ordinate


def maprzk2rhomu(a5, mass, charge, energy, rvals, zvals, ksivals,
                 rhogrid, mugrid, weights=None):
    """
    Convert (r,z,ksi) to (rho_omp,mu) via (P,mu).
    """
    free = False
    if isinstance(a5, str):
        from a5py.ascotpy.ascotpy import Ascotpy
        free = True
        a5 = Ascotpy(a5)
        a5.init(bfield=True)

    Pgrid  = np.linspace(-1,1,5)
    _, Pvals, muvals = maprzk2Pmu(a5, mass, charge, energy, rvals, zvals,
                                  ksivals, Pgrid, mugrid, weights=None)

    ordinate, rhovals, _ = mapPmu2rhomu(a5, mass, charge, energy, Pvals,
                                        muvals, rhogrid, mugrid,
                                        weights=weights)

    if free:
        a5.free(bfield=True)

    return (ordinate, rhovals, muvals)
