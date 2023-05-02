import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .helpers import openfigureifnoaxes

@openfigureifnoaxes(projection=None)
def scatter2d(x, y, c="C0", log=[False, False, False], nc=9, cmap="viridis",
              xlabel=None, ylabel=None, clabel=None, axesequal=False, axes=None,
              cax=None):
    """
    Make a scatter plot in 2D+1 where color can be one dimension.

    For better performance this function uses matplotlib's plot function instead
    of the scatter function. The only practical difference is that the marker
    size can't vary and the color can't be continuous.

    Args:
        x : array_like <br>
            Marker x-coordinates.
        y : array_like <br>
            Marker y-coordinates.
        c : {str, array_like}, optional <br>
            Color data or string indicating the color.
        log : [bool, bool, bool], optional <br>
            Make [x-axis, y-axis, color axis] logarithmic.
        nc : int, optional <br>
            Number of colors used if c contains data. Since we are using plot
            instead of data, the color scale can't be continuous.
        cmap : str, optional <br>
            Name of the colormap where nc colors are picked if c contains data.
        xlabel : str, optional <br>
            Label for the x-axis.
        ylabel : str, optional <br>
            Label for the y-axis.
        clabel : str, optional <br>
            Label for the color axis.
        axesequal : bool, optional <br>
            Flag to set aspect ratio of [x,y] axes equal.
        axes : Axes, optional <br>
            The Axes object to draw on.
        cax : Axes, optional <br>
            The Axes object for the color data (if c contains data), otherwise
            taken from axes.
    """

    cmap = plt.cm.get_cmap(cmap, nc)
    cbar = None

    if log[0]:
        x      = np.log10(np.absolute(x))
        xlabel = "log10(| " + xlabel + " |)"

    if log[1]:
        y      = np.log10(np.absolute(y))
        ylabel = "log10(| " + ylabel + " |)"

    if isinstance(c, str):
        # Simple plot with a single color
        axes.plot(x, y, color=c, linestyle="None", marker="o")
    else:
        if log[1]:
            c      = np.log10(np.absolute(c))
            clabel = "log10(| " + clabel + " |)"

        # Sort inputs by color values and then find the indices that divide the
        # color range in even intervals
        idx  = np.argsort(c)
        x    = x[idx]
        y    = y[idx]
        c    = c[idx]
        clim = np.linspace(c[0], c[-1], nc+1)
        idx  = np.searchsorted(c, clim) + 1

        # Now plot markers corresponding to each interval with a different color
        i1 = 0
        for i in range(nc):
            i2 = idx[i+1]
            axes.plot(x[i1:i2], y[i1:i2], color=cmap(i/nc), linestyle="None",
                      marker="o")
            i1 = i2

        # Make colorbar with where color scale has the intervals
        norm = mpl.colors.BoundaryNorm(clim, nc)
        smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = plt.colorbar(smap, ax=axes, cax=cax)
        cax  = cbar.ax

    # Apply decorations
    if axesequal:
        axes.set_aspect("equal", adjustable="box")

    axes.set_xlabel(xlabel);
    axes.set_ylabel(ylabel);
    axes.tick_params(axis="x", direction="out")
    axes.tick_params(axis="y", direction="out")

    if cbar is not None:
        cbar.set_label(clabel)

    return (axes, cax)


@openfigureifnoaxes(projection="3d")
def scatter3d(x, y, z, c="C0", nc=9, cmap="viridis",
              log=[False, False, False, False], xlabel=None, ylabel=None,
              zlabel=None, clabel=None, axesequal=False, axes=None, cax=None):
    """
    Make a scatter plot in 3D+1 where color can be one dimension.

    For better performance this function uses matplotlib's plot function instead
    of the scatter function. The only practical difference is that the marker
    size can't vary and the color can't be continuous.

    Args:
        x : array_like <br>
            Marker x-coordinates.
        y : array_like <br>
            Marker y-coordinates.
        z : array_like <br>
            Marker z-coordinates.
        c : {str, array_like}, optional <br>
            Color data or string indicating the color.
        log : [bool, bool, bool, bool], optional <br>
            Make [x-axis, y-axis, z-axis, color axis] logarithmic.
        nc : int, optional <br>
            Number of colors used if c contains data. Since we are using plot
            instead of data, the color scale can't be continuous.
        cmap : str, optional <br>
            Name of the colormap where nc colors are picked if c contains data.
        xlabel : str, optional <br>
            Label for the x-axis.
        ylabel : str, optional <br>
            Label for the y-axis.
        zlabel : str, optional <br>
            Label for the z-axis.
        clabel : str, optional <br>
            Label for the color axis.
        axesequal : bool, optional <br>
            Flag to set aspect ratio of [x,y,z] axes equal.
        axes : Axes, optional <br>
            The Axes object to draw on. If None, a new figure is displayed.
        cax : Axes, optional <br>
            The Axes object for the color data (if c contains data), otherwise
            taken from axes.
    """

    cmap = plt.cm.get_cmap(cmap, nc)
    cbar = None

    if log[1]:
        x      = np.log10(np.absolute(x))
        xlabel = "log10(| " + xlabel + " |)"
    if log[1]:
        y      = np.log10(np.absolute(y))
        ylabel = "log10(| " + ylabel + " |)"
    if log[1]:
        z      = np.log10(np.absolute(z))
        zlabel = "log10(| " + zlabel + " |)"

    if isinstance(c, str):
        # Simple plot with a single color
        axes.plot(x, y, z, color=c, linestyle="None", marker="o")
    else:
        if log[1]:
            c      = np.log10(np.absolute(c))
            clabel = "log10(| " + clabel + " |)"

        # Sort inputs by color values and then find the indices that divide the
        # color range in even intervals
        idx  = np.argsort(c)
        x    = x[idx]
        y    = y[idx]
        z    = z[idx]
        c    = c[idx]
        clim = np.linspace(c[0], c[-1], nc+1)
        idx  = np.searchsorted(c, clim) + 1

        # Now plot markers corresponding to each interval with a different color
        i1 = 0
        for i in range(nc):
            i2 = idx[i+1]
            axes.plot(x[i1:i2], y[i1:i2], z[i1:i2], color=cmap(i/nc),
                      linestyle="None", marker="o")
            i1 = i2

        # Make colorbar with where color scale has the intervals
        norm = mpl.colors.BoundaryNorm(clim, nc)
        smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = plt.colorbar(smap, ax=axes, cax=cax)
        cax  = cbar.ax

    # Apply decorations
    if axesequal:
        axes.set_aspect("equal", adjustable="box")

    axes.set_xlabel(xlabel);
    axes.set_ylabel(ylabel);
    axes.set_zlabel(ylabel);
    axes.tick_params(axis="x", direction="out")
    axes.tick_params(axis="y", direction="out")
    axes.tick_params(axis="z", direction="out")

    if cbar is not None:
        cbar.set_label(clabel)
