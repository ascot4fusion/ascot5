"""
Routines for plotting markers.

File: marker/plot.py
"""
import numpy as np
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

def plot_line(x, y=None, z=None, ids=None, mask=None, equal=False,
              xlabel=None, ylabel=None, axes=None):
    """
    Plot continuous line f(x).

    Args:
        x : array_like <br>
            
        y : array_like, optional <br>
            
        z : array_like, optional <br>
        
        ids : array_like, optional <br>

        mask : array_like, optional <br>
            
        equal : bool, optional <br>
            
        xlabel : str, optional <br>

        ylabel : str, optional <br>

        zlabel : str, optional <br>

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

    if mask is not None:
        mask = np.ones(x.shape) == 1

    if ids is not None:
        uids = np.unique(ids)
        for i in uids:
            idx = np.where( (i==ids) & (i==ids) )[0]
            print(idx)
            axes.plot(x, y)
    else:
        axes.plot(x[mask], y[mask])

    if equal:
        axes.axis("scaled")

    if newfig:
        plt.show(block=False)

    return axes

def plot_histogram(x, bins=None, weights=None, logy=None, xlabel=None,
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
        bins = np.amax(10, np.floor(x.size / 10))
        bins = np.amin(100, bins)

    axes.hist(x, bins, normed=True, stacked=True, log=logy,
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

def plot_scatter(x, y=None, z=None, equal=None, xlabel=None, ylabel=None,
                 axes=None):
    """
    Plot a scatter plot.

    Args:
        x : array_like <br>
            
        y : array_like, optional <br>
            
        z : array_like, optional <br>
        
        ids : array_like, optional <br>

        mask : array_like, optional <br>
            
        equal : bool, optional <br>
            
        xlabel : str, optional <br>

        ylabel : str, optional <br>

        zlabel : str, optional <br>
        axes : Axes, optional <br>
            Axes where plot is plotted. If None, a new figure is created.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    if equal:
        axes.axis("scaled")

    if newfig:
        plt.show(block=False)

    return axes
