"""
Routines for plotting distributions.

File: plot.py
"""
import numpy as np
import importlib

plt = importlib.util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

def plot_dist_1D(dist, axes=None):
    """
    Plot distribution as a 1D plot.

    This function assumes the given distribution is squeezed so that a single
    dimension remains.

    Args:
        dist : dict_like <br>
            Distribution dictionary.
        axes : Axes, optional <br>
            Axes to which the distribution is plotted. If None, a new figure is
            created and displayed.
    """

    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    ordinate = None
    if "density" in dist:
        ordinate = dist["density"]
    elif "histogram" in dist:
        ordinate = dist["histogram"]

    x = dist["abscissae"][0]

    axes.plot(dist[x], ordinate)
    axes.set_xlabel(x);
    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', direction='out')

    if newfig:
        plt.show(block=False)

def plot_dist_2D(dist, *args, axes=None):
    """
    Plot distribution as a 2D plot (pcolormesh).

    This function assumes the given distribution is squeezed so that only two
    dimensions remain.

    Args:
        dist : dict_like <br>
            Distribution dictionary.
        args : str, str <br>
            Name of the x and y coordinates e.g. "R", "z".
        axes : Axes, optional <br>
            Axes to which the distribution is plotted. If None, a new figure is
            created and displayed.
    """

    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    ordinate = None
    if "density" in dist:
        ordinate = dist["density"]
    elif "histogram" in dist:
        ordinate = dist["histogram"]

    if len(args) == 0:
        x = dist["abscissae"][0]
        y = dist["abscissae"][1]
    else:
        x = args[0]
        y = args[1]

    if ordinate.shape[0] == dist[x].size:
        ordinate = np.transpose(ordinate)

    axes.pcolormesh(dist[x + "_edges"], dist[y + "_edges"], ordinate)
    axes.set_xlabel(x);
    axes.set_ylabel(y);
    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', direction='out')

    if newfig:
        plt.show(block=False)
