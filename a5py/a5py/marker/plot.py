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

def plot_histogram(x, xbins=None, y=None, ybins=None, weights=None,
                   logscale=False, xlabel=None, ylabel=None, axes=None):
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

    if xbins is None:
        xbins = np.linspace(np.amin(x), np.amax(x), 10)

    if y is None:
        axes.hist(x, xbins, density=False, stacked=True, log=logscale,
                  weights=weights, rwidth=1)

        if xlabel is not None:
            axes.set_xlabel(xlabel)
        if weights is not None:
            axes.set_ylabel("Particles per bin")
        else:
            axes.set_ylabel("Markers per bin")

        if not logscale:
            axes.ticklabel_format(style="sci", axis="y", scilimits=(0,0))

    else:
        if ybins is None:
            ybins = np.linspace(np.amin(y), np.amax(y), 10)

        axes.hist2d(x, y, bins=[xbins, ybins], weights=weights)

    if newfig:
        plt.show(block=False)

    return axes

def plot_scatter(x, y=None, z=None, c=None, prune=1, equal=False, ids=None,
                 xlabel=None, ylabel=None, zlabel=None, axes=None,
                 **kwargs):
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
        prune : int, optional <br>
            Plot only x[::prune] data points.
        kwargs : <br>
            Same arguments as in plot_line.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        axes = plt.figure()
        axes = plt.gca()

    cmap = plt.cm.get_cmap("viridis", 5)
    if z is None and c is None:
        if ids is not None:
            uids  = np.unique(ids)
            cpick = np.arange(5)
            np.random.shuffle(cpick)
            for i in range(uids.size):
                c = np.asarray( [cmap( cpick[np.mod(i, 5)] )]
                                * np.sum( ids==uids[i] ) )
                h = axes.scatter(x[ids==uids[i]][::prune],
                                 y[ids==uids[i]][::prune], c=c[::prune],
                                 edgecolors="none", **kwargs)
        else:
            h = axes.scatter(x[::prune], y[::prune], **kwargs)

    elif z is not None and c is None:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                h = axes.scatter(x[ids==i][::prune], y[ids==i][::prune],
                                 zs=z[ids==i][::prune])

        else:
            h = axes.scatter(x, y, zs=z)

    elif z is None and c is not None:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                h = axes.scatter(x[ids==i][::prune], y[ids==i][::prune],
                                 c=c[ids==i][::prune])

        else:
            h = axes.scatter(x, y, c=c)
            plt.colorbar(h, ax=axes)
    else:
        if ids is not None:
            uids = np.unique(ids)
            for i in uids:
                h = axes.scatter(x[ids==i][::prune], y[ids==i][::prune],
                                 zs=z[ids==i][::prune], c=c[ids==i][::prune])

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
