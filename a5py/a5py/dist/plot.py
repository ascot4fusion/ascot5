"""
Routines for plotting distributions.

File: dist/plot.py
"""
import numpy as np
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

def plot_dist_1D(dist, logscale=False, axes=None):
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
    if "distribution" in dist:
        ordinate = dist["distribution"]
    elif "histogram" in dist:
        ordinate = dist["histogram"]

    x = dist["abscissae"][0]

    if logscale:
        ordinate = np.log10(ordinate)

    axes.plot(dist[x], ordinate)
    axes.set_xlabel(x);
    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', direction='out')

    if newfig:
        plt.show(block=False)


def plot_dist_2D(dist, *args, logscale=False, equal=False, axes=None):
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
    if "distribution" in dist:
        ordinate = dist["distribution"]
    elif "histogram" in dist:
        ordinate = dist["histogram"]

    if len(args) == 0:
        x = dist["abscissae"][0]
        y = dist["abscissae"][1]
    else:
        x = args[0]
        y = args[1]

    if x == dist["abscissae"][0]:
        ordinate = np.transpose(ordinate)

    if logscale:
        ordinate = np.log10(ordinate)

    ordinate = np.ma.masked_invalid(ordinate)
    mesh = axes.pcolormesh(dist[x], dist[y], ordinate,
                           vmin=np.nanmin(ordinate), vmax=np.nanmax(ordinate))

    # https://stackoverflow.com/a/16125413/190597 (Joe Kington)
    # and https://stackoverflow.com/a/35905483
    axes.patch.set(hatch='x', edgecolor=[0.9, 0.9, 0.9])
    plt.colorbar(mesh, ax=axes)

    axes.set_xlabel(x);
    axes.set_ylabel(y);
    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', direction='out')

    if equal:
        axes.axis("image")
    else:
        axes.axis("tight")

    if newfig:
        plt.show(block=False)
