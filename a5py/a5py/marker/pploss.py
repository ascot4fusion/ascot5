"""
Study particle losses in phase space.

This module contains tools to plot how losses originate in phase space, to
evaluate transport coefficients, and to use the coefficients for estimating
losses.

pploss.py
"""
import numpy as np
from scipy.interpolate import interp1d

import a5py.marker.phasespace as phasespace
from a5py.ascotpy.ascotpy import Ascotpy

import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import ListedColormap


def plotpp(x, y, xgrid, ygrid, quantity=None, addlosscontour=None,
           addtrappedcontour=None, mask=None, axes=None, xlabel=None,
           xlim=None, ylabel=None, ylim=None, displaycolormap=True):
    """
    Plot marker distribution in phase space.

    This plot is similar to histogram in that it takes a pair of coordinates,
    and value associated with each coordinate, divides those among the mesh
    and plots the mean value on each cell.

    Args:
        x : array_like <br>
            Marker x coordinates.
        y : array_like <br>
            Marker y coordinates.
        xgrid : array_like <br>
            Grid edges for x coordinate.
        ygrid : array_like <br>
            Grid edges for y coordinate.
        quantity : array_like, optional <br>
            Marker values which are plotted. If None, plot shows marker density.
        addlosscontour : array_like, optional <br>
            Boolean array indicating which markers were lost. If given, adds
            contours that show limits for 10 % and 90 % losses.
        addtrappedcontour : array_like, optional <br>
            Boolean array indicating which markers are trapped. If given, adds
            contour showing trapped-passing limit.
        mask : array_like, optional <br>
            Boolean array indicating which values are taken into account when
            plotting the quantity (i.e. plot only lost markers). All values are
            still used to process addlosscontour and addtrappedcontour.
        displaycolormap : bool, optional <br>
            Show colormap.
    """
    if axes is None:
        fig = plt.figure()
        gs = GridSpec(1,1)
        ax = fig.add_subplot(gs[0,0])
    else:
        ax = axes

    if mask is None:
        mask = np.ones(x.shape) == 1

    # Find cell counts
    ncount = np.histogram2d(x[mask], y[mask], bins=[xgrid,ygrid])[0]

    # Default quantity is marker density
    if quantity is None:
        quantity = np.ones(x.shape)
        ncount   = np.ones(ncount.shape)

    density = np.histogram2d(x[mask], y[mask], bins=[xgrid,ygrid],
                             weights=quantity[mask])[0] / ncount

    alpha=None
    if addlosscontour is not None:
        ncount   = np.histogram2d(x, y, bins=[xgrid,ygrid])[0]
        lossfrac = np.histogram2d(x, y, bins=[xgrid,ygrid],
                                  weights=addlosscontour)[0] / ncount

    mesh = ax.pcolormesh(xgrid, ygrid, np.log10(density.transpose()),
                         alpha=alpha)
    if displaycolormap:
        plt.colorbar(mesh)


    if addtrappedcontour is not None:
        ncount  = np.histogram2d(x, y, bins=[xgrid,ygrid])[0]
        trapped = np.histogram2d(x, y, bins=[xgrid,ygrid],
                                 weights=addtrappedcontour)[0] / ncount
        ax.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                   ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                   (trapped.transpose() > 0)*1, [1], alpha=0.7,
                   colors='grey', linestyles='dotted', linewidths=2)

    if addlosscontour is not None:
        ncount   = np.histogram2d(x, y, bins=[xgrid,ygrid])[0]
        lossfrac = np.histogram2d(x, y, bins=[xgrid,ygrid],
                                  weights=addlosscontour)[0] / ncount
        ax.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                   ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                   lossfrac.transpose(), [0.1, 0.9], colors=['black', 'white'],
                   alpha=0.7, linewidths=2)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(xlim)


    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(ylim)


    if axes is None:
        plt.show(block=False)


def plotlossmap(a5, mass, charge, energy, r, z ,pitch, rhogrid, ksigrid,
                time, lost, weights, axis):
    xgrid = rhogrid
    ygrid = ksigrid
    _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                           pitch, xgrid, ygrid,
                                           weights=None)
    P,mu = phasespace.evalPmu(a5, mass, charge, energy, r, z, pitch)
    trapped = phasespace.istrapped(a5, mass, charge, energy,
                                   P, mu, rmin=4)

    lostxy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                            weights=weights[lost])[0]
    timexy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                            weights=(time*weights)[lost])[0] / lostxy

    timexy = np.log10(timexy)

    weightxy = np.histogram2d(x, y, bins=[xgrid,ygrid], weights=weights)[0]
    lossfrac = lostxy / weightxy

    trapxy = np.histogram2d(x, y, bins=[xgrid,ygrid],
                            weights=trapped*weights)[0] / weightxy

    cmap = plt.cm.get_cmap("viridis_r", 5)

    alphacmap = cmap(np.arange(cmap.N))
    alphacmap[:,:] = 1
    alphacmap[:,-1] = np.linspace(0, 1, cmap.N)
    alphacmap = ListedColormap(alphacmap)

    mesh = axis.pcolormesh(xgrid, ygrid, timexy.transpose(),
                           cmap=cmap, vmin=-5, vmax=-1)
    axis.figure.canvas.draw()
    colors = mesh.get_facecolor()

    def alpha_to_white(color):
        white = np.array([1,1,1])
        alpha = color[-1]
        color = color[:-1]
        return alpha*color + (1 - alpha)*white

    colors[:,3] = lossfrac.transpose().ravel()
    colors = np.array([alpha_to_white(color) for color in colors])
    mesh.set_facecolor(colors)

    axis.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                 ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                 (trapxy.transpose() > 0)*1, [1], alpha=0.7,
                 colors='red', linestyles='dotted', linewidths=1)

    axis.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                 ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                 lossfrac.transpose(), [0.1, 0.9], colors=['black', 'black'],
                 alpha=0.7, linewidths=[1, 3])


def evalandplotpp(a5, mass, charge, energy, r, z, pitch, xy="rhoksi",
                  edges=None, addtrappedcontour=False, **kwargs):

    if xy == "rhoksi":
        if edges is None:
            edges = [np.linspace(0, 1, 10), np.linspace(-1,1,10)]

        _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                           pitch, edges[0], edges[1],
                                           weights=None)
        if "xlabel" not in kwargs:
            kwargs["xlabel"] = "rho"
        if "ylabel" not in kwargs:
            kwargs["ylabel"] = "pitch"

    if xy == "rhomu":
        if edges is None:
            edges = [None, None]
            edges[1] = phasespace.initgridPmu(a5, mass, charge, energy, nP=10,
                                              nmu=20, padding=1e-1)[1].ravel()
            edges[0] = np.linspace(0, 1, 10)

        _, x, y = phasespace.maprzk2rhomu(a5, mass, charge, energy, r, z,
                                          pitch, edges[0], edges[1],
                                          weights=None)

        if "xlabel" not in kwargs:
            kwargs["xlabel"] = "rho"
        if "ylabel" not in kwargs:
            kwargs["ylabel"] = "mu"

    if xy == "Pmu":
        if edges is None:
            edges = [None, None]
            edges[0], edges[1] = phasespace.initgridPmu(a5, mass, charge,
                                                        energy, nP=10, nmu=20,
                                                        padding=1e-1)
            edges[0] = edges[0].ravel()
            edges[1] = edges[1].ravel()

        _, x, y = phasespace.maprzk2Pmu(a5, mass, charge, energy, r, z,
                                          pitch, edges[0], edges[1],
                                          weights=None)

        if "xlabel" not in kwargs:
            kwargs["xlabel"] = "P"
        if "ylabel" not in kwargs:
            kwargs["ylabel"] = "mu"

    if addtrappedcontour == False:
        addtrappedcontour = None
    else:
        P,mu = phasespace.evalPmu(a5, mass, charge, energy, r, z, pitch)
        addtrappedcontour = phasespace.istrapped(a5, mass, charge, energy,
                                                 P, mu, rmin=4)

    plotpp(x, y, edges[0], edges[1], addtrappedcontour=addtrappedcontour,
           **kwargs)


def evaltransportcoefficients(a5, mass, charge, energy, r0, z0, k0,
                              xgrid, ygrid, drho, dtime, endtime, lost,
                              rhoedge=1, lossfrac=0.3, xy="rhoksi", plot=False):
    if xy == "rhoksi":
        _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r0, z0,
                                           k0, xgrid, ygrid, weights=None)

        xi = np.searchsorted(xgrid, x, side='left')
        yi = np.searchsorted(ygrid, y, side='left')

    K = np.zeros((xgrid.size-1,ygrid.size-1))
    D = np.zeros((xgrid.size-1,ygrid.size-1))
    for i in range(xgrid.size-1):
        for j in range(ygrid.size-1):
            ids = np.logical_and.reduce([xi == i, yi == j])

            if np.sum(ids) == 0:
                continue

            mean   = np.nanmean( drho[ids] )
            K[i,j] = np.nanmean( drho[ids] / dtime[ids] )
            D[i,j] = 0.5 * np.nanmean( ( drho[ids] - mean )**2 / dtime[ids] )

            losses = np.logical_and.reduce([lost, ids])
            if np.sum(losses)/np.sum(ids) > lossfrac:
                c1 = np.mean(endtime[losses])
                c2 = c1*c1*c1 / np.var(endtime[losses])

                delta  = 1 - ( xgrid[i] + xgrid[i+1] ) / 2
                K[i,j] = delta / c1
                D[i,j] = 2 * delta * delta / c2


    if plot:
        fig = plt.figure()
        gs = GridSpec(1,2)

        ax = fig.add_subplot(gs[0,0])
        mesh = ax.pcolormesh(xgrid, ygrid, K.transpose(), vmax=100)
        plt.colorbar(mesh, ax=ax)

        ax = fig.add_subplot(gs[0,1])
        mesh = ax.pcolormesh(xgrid, ygrid, D.transpose(), vmax=10)
        plt.colorbar(mesh, ax=ax)

        plt.show(block=False)


    return (xgrid, ygrid, K, D)


def runtransportmodel(coefficients, markers, tmax, Nt,
                      xmin=0, xmax=1, Nmrk=1000, axis=None):
    xgrid = coefficients[0]
    ygrid = coefficients[1]
    lossfrac = np.zeros((xgrid.size-1,ygrid.size-1))
    losstime = np.zeros((xgrid.size-1,ygrid.size-1))
    endstate = np.zeros((xgrid.size-1,ygrid.size-1))

    for j in range(ygrid.size-1):
        dt   = tmax/Nt

        x    = markers[0] + np.random.rand(Nmrk,1) * (markers[-1] - markers[0])
        x0   = x.copy()
        time = np.zeros(x.shape)
        lost = np.zeros(x.shape) == 1

        K = interp1d(xgrid[:-1] + (xgrid[0]-xgrid[0])/2,
                     coefficients[2][:,j], kind='nearest',
                     fill_value='extrapolate')
        D = interp1d(xgrid[:-1] + (xgrid[0]-xgrid[0])/2,
                     coefficients[3][:,j], kind='nearest',
                     fill_value='extrapolate')

        for i in range(Nt):
            # Prepare
            w = 1 - 2*np.random.rand(Nmrk,1)
            k = K(x)
            d = D(x)

            # Step
            x += (lost == False) * (k*dt + w*np.sqrt(2*d*dt))
            time[lost == False] += dt

            # Enforce boundary conditions
            reflect = x < xmin
            lost    = x >= xmax

            x[reflect] = 2*xmin - x[reflect]
            x[lost] = xmax

        for i in range(xgrid.size-1):
            ids = np.logical_and.reduce([x0 >= xgrid[i],
                                         x0 <  xgrid[i+1]])
            lossfrac[i,j] = np.sum(lost[ids])/np.sum(ids)

            ids = np.logical_and.reduce([lost, ids])
            losstime[i,j] = np.mean(time[ids])

            ids = np.logical_and.reduce([x >= xgrid[i],
                                         x <  xgrid[i+1]])
            endstate[i,j] = np.sum(ids)

    if axis is not None:
        cmap = plt.cm.get_cmap("viridis_r", 5)

        mesh = axis.pcolormesh(xgrid, ygrid, np.log10(losstime.transpose()),
                               cmap=cmap, vmin=-5, vmax=-1)

        axis.figure.canvas.draw()
        colors = mesh.get_facecolor()

        def alpha_to_white(color):
            white = np.array([1,1,1])
            alpha = color[-1]
            color = color[:-1]
            return alpha*color + (1 - alpha)*white

        colors[:,3] = lossfrac.transpose().ravel()
        colors = np.array([alpha_to_white(color) for color in colors])
        mesh.set_facecolor(colors)

        axis.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                 ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                 lossfrac.transpose(), [0.1, 0.9], colors=['black', 'black'],
                 alpha=0.7, linewidths=[1, 3])



def evalcoefs(a5, mass, charge, energy, r, z, pitch, weights, xgrid, ygrid,
              drift, diff, lost, endtime, xedge):
    _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                       pitch, xgrid, ygrid, weights=None)

    if weights is None:
        weights = np.ones(drift.shape)

    xi = np.searchsorted(xgrid, x, side='left')
    yi = np.searchsorted(ygrid, y, side='left')

    K = np.zeros((xgrid.size-1,ygrid.size-1))
    D = np.zeros((xgrid.size-1,ygrid.size-1))
    for i in range(xgrid.size-1):
        for j in range(ygrid.size-1):
            idx = np.logical_and.reduce([xi == i, yi == j])

            if np.sum(idx) == 0:
                # No markers in this cell
                continue

            iii = np.logical_and.reduce([idx, lost == False])
            if np.sum(lost[idx]) < np.sum(idx):
                #K[i,j] = np.nanmean(drift[iii])
                #D[i,j] = np.nanmean(diff[iii])
                K[i,j] = np.average(drift[iii], weights=weights[iii])
                D[i,j] = np.average(diff[iii], weights=weights[iii])
            else:
                K[i,j] = 0
                D[i,j] = 0

            if np.sum(lost[idx]) > 0:
                jjj = np.logical_and.reduce([lost, idx])
                c1 = np.average(endtime[jjj], weights=weights[jjj])

                v = np.average((endtime[jjj] - c1)**2, weights=weights[jjj])
                #c2 = c1*c1*c1 / np.var(endtime[jjj])
                c2 = c1*c1*c1/v
                delta  = xedge - ( xgrid[i] + xgrid[i+1] ) / 2

                frac = np.sum(lost[idx]) / np.sum(idx)

                K[i,j] = (1-frac) * K[i,j] + frac * (delta / c1)
                D[i,j] = (1-frac) * D[i,j] + frac * (2 * delta * delta / c2)


    return (xgrid, ygrid, K, D)
