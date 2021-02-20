"""
Study particle losses in phase space.

This module contains tools to plot how losses originate in phase space, to
evaluate transport coefficients, and to use the coefficients for estimating
losses.

pploss.py
"""
import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline, interp2d

import a5py.marker.phasespace as phasespace
from a5py.ascotpy.ascotpy import Ascotpy

import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import ListedColormap


def plotlossmap(a5, mass, charge, energy, r, z, pitch, rhogrid, ksigrid,
                time, lost, weights, rmin, axis, muin=False):
    """
    Plot lost markers in (rho_omp,ksi_omp) or (rho_omp,mu) space.

    This plot constructs a loss map for tokamaks. Loss map shows the fraction of
    lost particles and their mean loss time. Loss time is shown in color while
    the fraction of particles lost is illustrated with shade and black contours
    (10 % losses = thin contour, 90 % losses = thick contour).  The plot also
    shows the passing-trapped boundary in red.

    Args:
        a5 : Ascotpy <br>
            Ascotpy object with magnetic field initialized.
        mass : float <br>
            Marker mass [kg].
        charge : float <br>
            Marker charge [C].
        energy : float <br>
            Marker energy [J].
        r : array_like (n,1) <br>
            Marker major radius coordinates [m].
        z : array_like (n,1) <br>
            Marker z coordinates [m].
        pitch : array_like (n,1) <br>
            Marker pitch coordinates (v_para/v_tot).
        rhogrid : array_like <br>
            OMP rho grid values.
        ksigrid : array_like <br>
            OMP pitch grid values.
        time : array_like (n,1) <br>
            Marker end state time [s].
        lost : array_like (n,1) <br>
            Boolean array indicating which markers are lost.
        weights : array_like (n,1) <br>
            Marker weights.
        rmin : float <br>
            R coordinate for the inner mid plane separatrix.
        axis : Axis <br>
            Axis on which the lossmap is plotted.
        muin : bool, optional <br>
            mu grid if the lossmap is plotted in (rho_omp,mu) instead.
    """

    # Find the phase-space coordinates for each marker
    if muin:
        xgrid = rhogrid
        ygrid = ksigrid
        _, x, y = phasespace.maprzk2rhomu(a5, mass, charge, energy, r, z,
                                          pitch, xgrid, ygrid,
                                          weights=None)
    else:
        xgrid = rhogrid
        ygrid = ksigrid
        _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                           pitch, xgrid, ygrid,
                                           weights=None)

    # Which particles are trapped
    P,mu = phasespace.evalPmu(a5, mass, charge, energy, r, z, pitch)
    trapped = phasespace.istrapped(a5, mass, charge, energy,
                                   P, mu, rmin=rmin)

    # Evaluate the plotted quantities
    lostxy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                            weights=weights[lost])[0]
    timexy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                            weights=(time*weights)[lost])[0] / lostxy

    timexy = np.log10(timexy)

    weightxy = np.histogram2d(x, y, bins=[xgrid,ygrid], weights=weights)[0]
    lossfrac = lostxy / weightxy

    trapxy = np.histogram2d(x, y, bins=[xgrid,ygrid],
                            weights=trapped*weights)[0] / weightxy

    # Plot loss time with 4 different colours
    cmap = plt.cm.get_cmap("viridis_r", 4)
    alphacmap = cmap(np.arange(cmap.N))
    alphacmap[:,:] = 1
    alphacmap[:,-1] = np.linspace(0, 1, cmap.N)
    alphacmap = ListedColormap(alphacmap)

    mesh = axis.pcolormesh(xgrid, ygrid, timexy.transpose(),
                           cmap=cmap, vmin=-9, vmax=-5)
    axis.figure.canvas.draw()

    # Get the colors and adjust shade so it corresponds to fraction of prts lost
    colors = mesh.get_facecolor()
    def alpha_to_white(color):
        white = np.array([1,1,1])
        alpha = color[-1]
        color = color[:-1]
        return alpha*color + (1 - alpha)*white

    colors[:,3] = lossfrac.transpose().ravel()
    colors = np.array([alpha_to_white(color) for color in colors])
    mesh.set_facecolors(colors)

    # Passing-trapped boundary
    axis.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                 ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                 (trapxy.transpose() > 0)*1, [1], alpha=0.7,
                 colors='red', linestyles='dotted', linewidths=1)

    # 10 % and 90 % contours
    axis.contour(xgrid[:-1] + (xgrid[1]-xgrid[0])/2,
                 ygrid[:-1] + (ygrid[1]-ygrid[0])/2,
                 lossfrac.transpose(), [0.1, 0.9], colors=['black', 'black'],
                 alpha=0.7, linewidths=[1, 3])

    axis.figure.canvas.draw()


def evallossfrac(a5, mass, charge, energy, r, z ,pitch, rhogrid, ksigrid,
                 time, lost, weights, endenergy=None):
    xgrid = rhogrid
    ygrid = ksigrid
    _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                       pitch, xgrid, ygrid,
                                       weights=None)
    P,mu = phasespace.evalPmu(a5, mass, charge, energy, r, z, pitch)

    lostxy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                            weights=weights[lost])[0]
    timexy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                            weights=(time*weights)[lost])[0] / lostxy

    weightxy = np.histogram2d(x, y, bins=[xgrid,ygrid], weights=weights)[0]
    lossfrac = lostxy / weightxy

    if endenergy is None:
        return (lossfrac, timexy)
    else:
        lostxy = np.histogram2d(x[lost], y[lost], bins=[xgrid,ygrid],
                                weights=(weights*endenergy)[lost])[0]
        weightxy = np.histogram2d(x, y, bins=[xgrid,ygrid],
                                  weights=weights*energy)[0]
        return (lossfrac, timexy, lostxy / weightxy)


def evalpdens(a5, mass, charge, energy, r, z ,pitch, rhogrid, ksigrid, weights):
    """
    """
    xgrid = rhogrid
    ygrid = ksigrid
    _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                       pitch, xgrid, ygrid,
                                       weights=None)
    P,mu = phasespace.evalPmu(a5, mass, charge, energy, r, z, pitch)


    weightxy = np.histogram2d(x, y, bins=[xgrid,ygrid],
                              weights=weights*energy)[0]

    return weightxy


def runtransportmodel(coefficients, markers, tmax, Nt,
                      xmin=0, xmax=1, Nmrk=1000, axis=None, ksigrid=None,
                      a5=None, mass=None, charge=None, energy=None):
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

    if ksigrid is not None:
        xvals, yvals = np.meshgrid(xgrid[:-1]+(xgrid[1]-xgrid[0])/2,
                                   ksigrid[:-1]+(ksigrid[1]-ksigrid[0])/2,
                                   indexing='ij')
        _, x, y = phasespace.maprhoksi2rhomu(a5, mass, charge, energy,
                                             xvals.ravel(), yvals.ravel(),
                                             xgrid, ygrid, weights=None)

        xvals, yvals = np.meshgrid(xgrid[:-1]+(xgrid[1]-xgrid[0])/2,
                                   ygrid[:-1]+(ygrid[1]-ygrid[0])/2,
                                   indexing='ij')
        f = RectBivariateSpline(xgrid[:-1]+(xgrid[1]-xgrid[0])/2,
                                ygrid[:-1]+(ygrid[1]-ygrid[0])/2,
                                lossfrac)
        lossfrac = f(x,y,grid=False).reshape(xgrid.size-1, ksigrid.size-1)
        #f = interp2d(xvals.ravel(),
        #             yvals.ravel(),
        #             lossfrac.ravel())
        #lossfrac = (np.diag(f(x,y, grid=False)).reshape(ksigrid.size-1, xgrid.size-1)).T
        #print(lossfrac.shape)
        #lossfrac = f(xgrid, y)


        lossfrac[lossfrac>1] = 1
        lossfrac[lossfrac<0] = 0

        losstime[np.isnan(losstime)] = 0
        f = RectBivariateSpline(xgrid[:-1]+(xgrid[1]-xgrid[0])/2,
                                ygrid[:-1]+(ygrid[1]-ygrid[0])/2,
                                losstime)
        losstime = f(x,y,grid=False).reshape(xgrid.size-1, ksigrid.size-1)

        #f = interp2d(xvals,
        #             yvals,
        #             losstime)
        #losstime = f(xgrid[:-1]+(xgrid[1]-xgrid[0])/2, ygrid[:-1]+(ygrid[1]-ygrid[0])/2).T#f(xgrid, ygrid)#(x, y).reshape(xgrid.size-1, ksigrid.size-1)
        #losstime = (np.diag(f(x,y)).reshape(ksigrid.size-1, xgrid.size-1)).T
        #losstime.setflags(write=1)

        ygrid = ksigrid

        #n = np.histogram2d(x, y, bins=[xgrid,ygrid])[0]
        #lossfrac = np.histogram2d(x, y, bins=[xgrid,ygrid],
        #                          weights=lossfrac.ravel())[0] / n
        #losstime = np.histogram2d(x, y, bins=[xgrid,ygrid],
        #                          weights=losstime.ravel())[0] / n
        #lossfrac[np.isnan(lossfrac)] = 0
        #losstime[np.isnan(losstime)] = 1


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


    return (lossfrac,losstime)



def evalcoefs(a5, mass, charge, energy, r, z, pitch, weights, xgrid, ygrid,
              drift, diff, lost, endtime, xedge, mutype=False):
    if mutype:
        _, x, y = phasespace.maprzk2rhomu(a5, mass, charge, energy, r, z,
                                          pitch, xgrid, ygrid, weights=None)
    else:
        _, x, y = phasespace.maprzk2rhoksi(a5, mass, charge, energy, r, z,
                                           pitch, xgrid, ygrid, weights=None)

    if weights is None:
        weights = np.ones(drift.shape)

    xi = np.searchsorted(xgrid, x, side='left')-1
    yi = np.searchsorted(ygrid, y, side='left')-1

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
