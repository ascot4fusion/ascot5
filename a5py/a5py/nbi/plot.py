"""
Routines for plotting nbi inputs.

File: nbi/plot.py

"""
import numpy as np
from numpy.random import randint
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

def plot_arrow_3D(x, y, z, dx, dy, dz, equal=False,
                  xlabel=None, ylabel=None, zlabel=None, axes=None, **kwargs):
    """
    Plot arrows in 3D plot.
    
    Args:
        x : array_like <br>
            x data.
        y : array_like <br>
            y data.
        z : array_like <br>
            z data.
        dx : array_like <br>
            dx data.
        dy : array_like <br>
            dy data.
        dz : array_like <br>
            dz data.
        equal : bool, optional <br>
            Flag for making axes equal.
        xlabel : str, optional <br>
            Label on x axis.
        ylabel : str, optional <br>
            Label on y axis.
        zlabel : str, optional <br>
            Label on z axis.
        axes : Axes, optional <br>
            Axes where plot is plotted. Remember it must have 3D projections enabled.
            If None, a new figure is created.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        fig = plt.figure()
        axes = fig.add_subplot(1,1,1,projection="3d")

    axes.quiver3D(x,y,z,
                  dx,dy,dz, **kwargs)

    if equal:
        axes.axis("scaled")

    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    if zlabel is not None:
        axes.set_zlabel(zlabel)

    if newfig:
        plt.show(block=False)

    return axes

## Continue here

def plot_scatter_3D(x, y, z, equal=False, xlabel=None, ylabel=None, zlabel=None,
                    axes=None, **kwargs):
    """
    Plot a scatter plot.

    Args:
        x : array_like <br>
            x data.
        y : array_like, optional <br>
            y data.
        z : array_like, optional <br>
            z data. If None, scatter will be plotted in 2D.
        equal : bool, optional <br>
            Flag for making axes equal.
        xlabel : str, optional <br>
            Label on x axis.
        ylabel : str, optional <br>
            Label on y axis.
        zlabel : str, optional <br>
            Label on z axis.
        axes : Axes, optional <br>
            Axes where plot is plotted. Remember it must have 3D projections enabled.
            If None, a new figure is created.

    Returns:
        Axes where plot is plotted.
    """
    newfig = axes is None
    if newfig:
        fig = plt.figure()
        axes = fig.add_subplot(1,1,1,projection="3d")

    axes.scatter(x,y,z, **kwargs)

    if equal:
        axes.axis("scaled")

    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    if zlabel is not None:
        axes.set_zlabel(zlabel)

    if newfig:
        plt.show(block=False)

    return axes
