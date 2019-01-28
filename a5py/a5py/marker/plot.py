"""
Routines for plotting markers.

File: marker/plot.py
"""
import numpy as np
from numpy.random import randint
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

def plot_line(x, y=None, z=None, ids=None, equal=False,
              xlabel=None, ylabel=None, zlabel=None, axes=None):
    """
    Plot continuous line f(x).

    Args:
        x : array_like <br>
            x data.
        y : array_like, optional <br>
            y data.
        z : array_like, optional <br>
            z data, in which case plot will be 3D.
        ids : array_like, optional <br>
            Array with same shape as input data containing marker ID for each
            point. If this is provided, each marker is plotted individually.
        equal : bool, optional <br>
            Flag for making axes equal.
        xlabel : str, optional <br>
            Label on x axis.
        ylabel : str, optional <br>
            Label on y axis.
        zlabel : str, optional <br>
            Label on z axis.
        axes : Axes, optional <br>
            Axes where plot is plotted. If None, a new figure is created.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    if y is None:
        y = np.linspace(0, x.size, x.size)

    # If markers are separated by their ids, call this function recursively
    # for each marker.
    if ids is not None:
        uids = np.unique(ids)
        for i in uids:
            if z is None:
                plot_line(x[i==ids], y[i==ids], equal=equal, xlabel=xlabel,
                          ylabel=ylabel, axes=axes)
            else:
                plot_line(x[i==ids], y[i==ids], z[i==ids], equal=equal,
                          xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                          axes=axes)
    else:
        if z is None:
            axes.plot(x, y)
        else:
            axes.plot(x, y, z)

    if equal:
        axes.axis("scaled")

    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    if z is not None and zlabel is not None:
        axes.set_zlabel(zlabel)

    if newfig:
        plt.show(block=False)

    return axes

def plot_histogram(x, bins=None, weights=None, logy=False, xlabel=None,
                   axes=None):
    """
    Plot histogram.

    Args:
        x : array_like <br>
            Histogram quantity.
        bins : int or array_like, optional <br>
            Number of bins or array storing bin edges.
        weights : array_like, optional <br>
            Weights for x.
        logy : bool, optional <br>
            Make y-axis logarithmic.
        xlabel : str, optional <br>
            x axis label.
        axes : Axes, optional <br>
            Axes where plot is plotted. If None, a new figure is created.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    if bins is None:
        bins = 10
        #bins = np.amax(10, np.floor(x.size / 10))
        #bins = np.amin(100, bins)

    axes.hist(x, bins, normed=False, stacked=True, log=logy,
              weights=weights, rwidth=1)

    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if weights is not None:
        axes.set_ylabel("Particles per bin")
    else:
        axes.set_ylabel("Markers per bin")

    if not logy:
        axes.ticklabel_format(style="sci", axis="y", scilimits=(0,0))

    if newfig:
        plt.show(block=False)

    return axes

def plot_scatter(x, y=None, z=None, c=None, equal=False, ids=None,
                 xlabel=None, ylabel=None, zlabel=None, axes=None):
    """
    Plot a scatter plot.

    Args:
        x : array_like <br>
            x data.
        y : array_like, optional <br>
            y data.
        z : array_like, optional <br>
            z data. If None, scatter will be plotted in 2D.
        c : array_like, optional <br>
            Color data.
        kwargs : <br>
            Same arguments as in plot_line.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        axes = plt.figure()
        axes = plt.gca()

    cmap = plt.cm.get_cmap("hsv", 10)
    if z is None and c is None:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                k = randint(0, 10)
                h = axes.scatter(x[ids==i], y[ids==i], c=cmap(k), edgecolors="none")
        else:
            h = axes.scatter(x, y)

    elif z is not None and c is None:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                h = axes.scatter(x[ids==i], y[ids==i], zs=z[ids==i])

        else:
            h = axes.scatter(x, y, zs=z)

    elif z is None and c is not None:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                h = axes.scatter(x[ids==i], y[ids==i], c=c[ids==i])

        else:
            h = axes.scatter(x, y, c=c)
            plt.colorbar(h, ax=axes)
    else:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                h = axes.scatter(x[ids==i], y[ids==i], zs=z[ids==i],
                                 c=c[ids==i])

        else:
            h = axes.scatter(x, y, zs=z, c=c)
            plt.colorbar(h, ax=axes)

    if equal:
        axes.axis("scaled")

    axes.set_xlabel(xlabel);
    axes.set_ylabel(ylabel);
    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', direction='out')

    if z is not None:
        axes.set_zlabel(zlabel);
        axes.tick_params(axis='z', direction='out')

    if newfig:
        plt.show(block=False)

    return axes
