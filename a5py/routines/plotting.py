"""Module for generating common plots with ASCOT5.

This module should be imported everywhere where plotting is done instead of
using matplotlib directly. The reason is that on some platforms matplotlib
is not available and even there we want to able to use all functionality
that doesn't require plotting.

So either import this module or use try-except when importing matplotlib.
"""
import numpy as np
import warnings
import unyt
import a5py.physlib as physlib

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import mpl_toolkits
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    warnings.warn("Could not import matplotlib. Plotting disabled.")
    plt = None
    mpl = None

try:
    import pyvista as pv
except ImportError:
    warnings.warn("Could not import pyvista. 3D wall plotting disabled.")
    pv = None

try:
    from scipy.spatial import cKDTree
    import alphashape
except ImportError:
    warnings.warn(
        "Could not import cKDTree and alphashape. "
        "3D surface plotting disabled."
    )

from functools import wraps

def setpaperstyle(latex=True):
    """Set default figure settings (label sizes etc.) so that the figure is
    suitable for publications (looks nice on A4).

    This function modifies the matplotlib style settings so one call is changes
    the style for the entire session.

    Parameters
    ----------
    latex : bool, optional
        Use LaTex interpreter.
    """
    mpl.style.use({
        "figure.autolayout":False,
        "font.family":"serif",
        "pdf.fonttype":42,
        "ps.fonttype":42,
        "axes.labelsize":14,
        "axes.titlesize":16,
        "axes.titlepad":6,
        "xtick.labelsize":12,
        "ytick.labelsize":12,
        "axes.labelpad":6,
        "legend.fontsize":12,
        "legend.numpoints":1,
        "legend.scatterpoints":1,
        "grid.linewidth":0.8,
        "lines.linewidth":1.4,
        "patch.linewidth":0.24,
        "lines.markersize":5.6,
        "lines.markeredgewidth":0,
        "xtick.major.width":0.8,
        "ytick.major.width":0.8,
        "xtick.minor.width":0.4,
        "ytick.minor.width":0.4,
        "xtick.major.pad":5.6,
        "ytick.major.pad":5.6,
        "savefig.dpi":300,
        "axes.formatter.limits":[-2,2]
    })
    if latex:
        mpl.style.use({
            "font.serif":"ComputerModern",
            "text.usetex":True,
        })

def setguistyle(latex=True):
    """Set default figure settings (label sizes etc.) so that the figure is
    suitable for GUI and presentations (large labels).

    This function modifies the matplotlib style settings so one call is changes
    the style for the entire session.

    Parameters
    ----------
    latex : bool, optional
        Use LaTex interpreter.
    """
    mpl.style.use({
        "figure.autolayout":False,
        "font.family":"serif",
        "axes.labelsize":18,
        "axes.titlesize":18,
        "axes.titlepad":12,
        "xtick.labelsize":16,
        "ytick.labelsize":16,
        "axes.labelpad":6,
        "legend.fontsize":16,
        "legend.numpoints":1,
        "legend.scatterpoints":1,
        "grid.linewidth":0.8,
        "lines.linewidth":1.4,
        "patch.linewidth":0.24,
        "lines.markersize":5.6,
        "lines.markeredgewidth":0,
        "xtick.major.width":0.8,
        "ytick.major.width":0.8,
        "xtick.minor.width":0.4,
        "ytick.minor.width":0.4,
        "xtick.major.pad":5.6,
        "ytick.major.pad":5.6,
        "savefig.dpi":300,
        "axes.formatter.limits":[-2,2]
    })
    if latex:
        mpl.style.use({
            "font.serif":"ComputerModern",
            "text.usetex":True,
        })

def figuresinglecolumn(aspectratio=3/2):
    """Return figure that has a size suitable for printing in A4 single-column
    width (when the paper has a double-column format).

    Parameters
    ----------
    aspectratio : float
        Width / height ratio of the returned figure.

    Returns
    -------
    """
    return plt.figure(figsize=(3.504, 3.504/aspectratio))

def figuredoublecolumn(aspectratio=3/2):
    """Return figure that has a size suitable for printing in A4 double-column
    width (when the paper has a double-column format).
    """
    return plt.figure(figsize=(7.205, 7.205/aspectratio))

def openfigureifnoaxes(projection="rectilinear"):
    """Decorator for creating and displaying a new figure if axes are not
    provided.

    Parameters
    ----------
    projection : str, {None, '3d', 'aitoff', 'hammer', 'lambert', 'mollweide',
    'polar', 'rectilinear'}, optional
        The axes projection type, see Matplotlib documentation for details.
    """
    def actualdecorator(plotfun):
        """Decorator that takes the plotting routine.
        """
        @wraps(plotfun)
        def wrapper(*args, axes=None, **kwargs):
            """Create a new figure if axes is None and pass *args and **kwargs
            for the plotter.
            """

            # Sometimes you don't want to give axes, but you want to call
            # several plotting functions for the same axes. Hence, you need to
            # skip plt.show() initially.
            skipshow = kwargs.pop("skipshow", False)
            if axes is None:
                fig  = plt.figure()
                axes = fig.add_subplot(111, projection=projection)
                result = plotfun(*args, axes=axes, **kwargs)
                if not skipshow:
                    plt.show()
                return result
            else:
                if projection != None and axes.name != projection:
                    raise ValueError(
                        "Invalid projection \"%s\" on axes: expected \"%s\"" %
                        (axes.name, projection))
                return plotfun(*args, axes=axes, **kwargs)

        return wrapper

    return actualdecorator

def getmathtextsciformatter(format):
    """Returns a label tick formatter that shows numbers in format "a x 10^b".

    Credit: https://stackoverflow.com/a/49330649

    Examples
    --------
    >>> plt.gca().yaxis.set_major_formatter(getmathtextsciformatter("%1.2e"))
    """
    class MathTextSciFormatter(mpl.ticker.Formatter):

        def __init__(self, format="%1.2e"):
            self.fmt = format

        def __call__(self, x, pos=None):
            s = self.fmt % x
            decimal_point = '.'
            positive_sign = '+'
            tup = s.split('e')
            significand = tup[0].rstrip(decimal_point)
            sign = tup[1][0].replace(positive_sign, '')
            exponent = tup[1][1:].lstrip('0')
            if exponent:
                exponent = '10^{%s%s}' % (sign, exponent)
            if significand and exponent:
                s =  r'%s{\times}%s' % (significand, exponent)
            else:
                s =  r'%s%s' % (significand, exponent)
            return "${}$".format(s)

    return MathTextSciFormatter(format)

def show():
    return plt.show()

def legend():
    return plt.legend()

def tight_layout(axes):
    fig = axes.get_figure()
    return fig.tight_layout()

def subplots():
    return plt.subplots()

@openfigureifnoaxes(projection=None)
def scatter2d(x, y, c=None, xlog="linear", ylog="linear", clog="linear",
              xlabel=None, ylabel=None, clabel=None, cint=9, cmap=None,
              axesequal=False, axes=None, cax=None, marker="o", markersize=6.0,
              markerfacecolor=None, alpha=1.0,
              title=None, label=None, skipshow=False):
    """Make a scatter plot in 2D+1 where color can be one dimension.

    For better performance this function uses matplotlib's plot function instead
    of the scatter function. The only practical difference is that the marker
    size can't vary and the color can't be continuous.

    Parameters
    ----------
    x : array_like
        Marker x-coordinates.
    y : array_like
        Marker y-coordinates.
    c : str or array_like, optional
        Color data or string indicating the color.
    xlog : {"linear", "log", "symlog"}, optional
        x-axis scaling.
    ylog : {"linear", "log", "symlog"}, optional
        y-axis scaling.
    clog : {"linear", "log", "symlog"}, optional
        color-axis scaling.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    cint : int, optional
        Number of colors used if c contains data. Since we are using plot
        instead of data, the color scale can't be continuous.
    cmap : str, optional
        Name of the colormap where nc colors are picked if c contains data.
    axesequal : bool, optional
        Flag to set aspect ratio of [x,y] axes equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    markersize : float, optional
        size of the drawn markers
    markerfacecolor : str, optional
        color with which the drawn markers are filled
    alpha : float, optional
        Transparency (0: transparent, 1:opaque)
    title : str, optional
        title for the axes
    label : str, optional
        label for the legend
    skipshow : bool, optional
        Skip plt.show() if True. Only used by the decorator when no axes given.
    """
    axes.set_xscale(xlog)
    axes.set_yscale(ylog)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if axesequal: axes.set_aspect("equal", adjustable="box")
    if title is not None: axes.set_title(title)

    if c is None or isinstance(c, str):
        # Simple plot with a single color
        axes.plot(x, y, color=c, linestyle="None", marker=marker,
                  markersize=markersize, alpha=alpha,
                  markerfacecolor=markerfacecolor, label=label)
        return axes

    # Sort inputs by color values and then find the indices that divide the
    # color range in given intervals
    idx  = np.argsort(c)
    x    = x[idx]
    y    = y[idx]
    c    = c[idx]
    if isinstance(cint, int):
        cint = np.linspace(c[0], c[-1], cint+1)
    idx  = np.searchsorted(c, cint) + 1
    nc = cint.size - 1

    # Choose default colormap depending on whether values change sign, and pick
    # nc colors
    if cmap is None:
        cmap = "viridis"
        if any(c < 0) and any(c > 0): cmap = "bwr"
    cmap = plt.get_cmap(cmap, nc)

    # Now plot markers corresponding to each interval with a different color
    i1 = 0
    for i in range(nc):
        i2 = idx[i+1]
        axes.plot(x[i1:i2], y[i1:i2], color=cmap(i/nc), linestyle="None",
                  marker=marker, markersize=markersize, alpha=alpha,
                  markerfacecolor=markerfacecolor, label=label)
        i1 = i2

    # Make colorbar with where color scale has the intervals
    norm = mpl.colors.BoundaryNorm(cint, nc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    ticks = []
    for b in cint:
        log = np.floor(np.log10(np.abs(b)))
        mul = b / 10**log
        ticks.append(r"$%.2f\times 10^{%.0f}$" % (mul, log))
    cbar.ax.set_yticklabels(ticks)
    cbar.set_label(clabel)
    return axes

@openfigureifnoaxes(projection="3d")
def scatter3d(x, y, z, c=None, xlog="linear", ylog="linear", zlog="linear",
              clog="linear", xlabel=None, ylabel=None, zlabel=None, clabel=None,
              cint=None, cmap=None, axesequal=False, axes=None, cax=None,
              markersize=6.0, alpha=1.0, title=None):
    """Make a scatter plot in 3D+1 where color can be one dimension.

    For better performance this function uses matplotlib's plot function instead
    of the scatter function. The only practical difference is that the marker
    size can't vary and the color can't be continuous.

    Parameters
    ----------
    x : array_like
        Marker x-coordinates.
    y : array_like
        Marker y-coordinates.
    z : array_like
        Marker z-coordinates.
    c : str or array_like, optional
        Color data or string indicating the color.
    xlog : {"linear", "log", "symlog"}, optional
        x-axis scaling.
    ylog : {"linear", "log", "symlog"}, optional
        y-axis scaling.
    zlog : {"linear", "log", "symlog"}, optional
        z-axis scaling.
    clog : {"linear", "log", "symlog"}, optional
        color-axis scaling.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    zlabel : str, optional
        Label for the z-axis.
    clabel : str, optional
        Label for the color axis.
    cint : int, optional
        Number of colors used if c contains data. Since we are using plot
        instead of data, the color scale can't be continuous.
    cmap : str, optional
        Name of the colormap where nc colors are picked if c contains data.
    axesequal : bool, optional
        Flag to set aspect ratio of [x,y] axes equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    markersize: float, optional
        Sets the size of the markers. Bigger number, bigger marker.
    alpha : float, optional
        Transparency (0: transparent, 1:opaque)
    title : str, optional
        title for the axes
    """
    axes.set_xscale(xlog)
    axes.set_yscale(ylog)
    axes.set_zscale(zlog)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_zlabel(zlabel)
    if axesequal:
        # axesequal=True requires manually setting symmetric axis limits;
        # set_aspect("equal") alone does not enforce equal scaling in 3D.
        all_points = np.array([x, y, z])
        mins = np.min(all_points, axis=1)
        maxs = np.max(all_points, axis=1)
        centers = (mins + maxs) / 2
        max_range = (maxs - mins).max() / 2

        axes.set_xlim(centers[0] - max_range, centers[0] + max_range)
        axes.set_ylim(centers[1] - max_range, centers[1] + max_range)
        axes.set_zlim(centers[2] - max_range, centers[2] + max_range)
        axes.set_aspect("equal", adjustable="box")

    if title is not None: axes.set_title(title)

    if c is None or isinstance(c, str):
        # Simple plot with a single color
        axes.plot(x, y, z, color=c, linestyle="None", marker="o",
                  markersize=markersize, alpha=alpha)
        return

    # Sort inputs by color values and then find the indices that divide the
    # color range in even intervals
    idx  = np.argsort(c)
    x    = x[idx]
    y    = y[idx]
    z    = z[idx]
    c    = c[idx]
    if cint is None:
        cint = np.linspace(c[0], c[-1], 10)
    elif isinstance(cint, int):
        cint = np.linspace(c[0], c[-1], cint+1)
    idx  = np.searchsorted(c, cint) + 1
    nc = cint.size - 1

    # Choose default colormap depending on whether values change sign, and pick
    # nc colors
    if cmap is None:
        cmap = "viridis"
        if any(c < 0) and any(c > 0): cmap = "bwr"
    cmap = plt.get_cmap(cmap, nc)

    # Now plot markers corresponding to each interval with a different color
    i1 = 0
    for i in range(nc):
        i2 = idx[i+1]
        axes.plot(x[i1:i2], y[i1:i2], z[i1:i2], color=cmap(i/nc),
                  linestyle="None", marker="o", markersize=markersize,
                  alpha=alpha)
        i1 = i2

    # Make colorbar with where color scale has the intervals
    norm = mpl.colors.BoundaryNorm(cint, nc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)

    ticks = []
    for b in cint:
        log = np.floor(np.log10(np.abs(b)))
        mul = b / 10**log
        ticks.append(r"$%.2f\times 10^{%.0f}$" % (mul, log))
    cbar.ax.set_yticklabels(ticks)
    cbar.set_label(clabel)

@openfigureifnoaxes(projection=None)
def hist1d(x, xbins=None, weights=None, xlog="linear", logscale=False,
           xlabel=None, ylabel=None, title=None, legend=None, axes=None,skipshow=False, histtype="bar"):
    """Plot (stacked) marker histogram in 1D.

    Parameters
    ----------
    x : array_like or list [array_like]
        Array or a list of arrays to be binned and plotted.

        If list is given, the resulting histogram is stacked with each
        array in the list corresponding to one layer in the stacked histogram.
    xbins : int or array_like, optional
        Number of bins or array storing bin edges for the x coordinate.
    weights : array_like or list [array_like], optional
        Values the datapoints are weighted with.

        Same format as for x.
    xlog : {"linear", "log"}, optional
        x-axis scaling.
    logscale : bool, optional
        Show histogram in logarithmic scale.
    xlabel : str, optional
        Label for the x-axis.
    legend : str, array_like
        List of strings to label the data with.

        The length of the list must be same as the number of data arrays
        provided.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    """
    axes.set_xlabel(xlabel)
    axes.set_xscale(xlog)
    axes.set_title(title)
    if ylabel is None:
        ylabel = "Markers per bin" if weights is None else "Particles per bin"
    axes.set_ylabel(ylabel)
    if not logscale:
        axes.ticklabel_format(style="sci", axis="y", scilimits=(0,0))

    if isinstance(xbins, int) and xlog == "log":
        xmin = x[0][0]; xmax = x[0][0]
        for i in range(len(x)):
            xmin = np.amin([np.amin(x[i]), xmin])
            xmax = np.amax([np.amax(x[i]), xmax])
        xbins = np.logspace(np.log10(xmin), np.log10(xmax), xbins)

    # Plot and legend
    axes.hist(x, xbins, density=False, stacked=True, log=logscale,
              weights=weights, rwidth=2, histtype=histtype)
    if legend is not None:
        axes.legend(legend, frameon=False)
    return axes

@openfigureifnoaxes(projection=None)
def hist2d(x, y, xbins=None, ybins=None, weights=None, xlog="linear",
           ylog="linear", logscale=False, xlabel=None, ylabel=None,
           axesequal=False, axes=None, cax=None):
    """Plot marker histogram in 2D.

    Parameters
    ----------
    x : array_like, (n,)
        x-coordinates of the data to be binned and plotted.
    y : array_like, (n,)
        y-coordinates of the data to be binned and plotted.
    xbins : int or array_like, optional
        Number of bins or array storing bin edges for the x coordinate.
    ybins : int or array_like, optional
        Number of bins or array storing bin edges for the y coordinate.
    weights : array_like, (n,), optional
        Values the datapoints are weighted with.

        If weights are included, the colorbar label changes from "Markers"
        to "Particles".
    xlog : {"linear", "log"}, optional
        x-axis scaling.
    ylog : {"linear", "log"}, optional
        y-axis scaling.
    logscale : bool, optional
        Show histogram in logarithmic scale.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    axesequal : bool, optional
        Flag to set the aspect ratio equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    axes.set_xlabel(xlabel)
    axes.set_xscale(xlog)
    axes.set_ylabel(ylabel)
    axes.set_yscale(ylog)
    clabel = "Markers per bin" if weights is None else "Particles per bin"
    if axesequal: axes.set_aspect("equal", adjustable="box")

    if isinstance(xbins, int) and xlog == "log":
        xmin = np.amin(x)
        xmax = np.amax(x)
        xbins = np.logspace(np.log10(xmin), np.log10(xmax), xbins)

    if isinstance(ybins, int) and ylog == "log":
        ymin = np.amin(y)
        ymax = np.amax(y)
        ybins = np.logspace(np.log10(ymin), np.log10(ymax), ybins)

    norm = None
    if logscale: norm = mpl.colors.LogNorm()

    if weights is not None:
        try:
            weights = weights.v # Cannot have units in weights yet in 2D histogram
        except AttributeError:
            pass
    h,_,_,m = axes.hist2d(x, y, bins=[xbins, ybins], weights=weights, norm=norm)

    cbar = plt.colorbar(m, ax=axes, cax=cax)
    cbar.set_label(clabel)

@openfigureifnoaxes(projection=None)
def mesh1d(x, y, xlabel=None, ylabel=None, axes=None,
           logscale=False, label=None):
    """Plot 1D distribution.

    Parameters
    ----------
    x : array_like (nx,)
        Abscissa edges for the x-axis.
    y : array_like (nx-1,)
        Data to be plotted.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    logscale: bool, optional
        Whether the plot is in logarithmic scale.
    label : str, optional
        Label if you are using a legend.
    """
    xc = np.zeros((y.size*2,))
    xc[1:-1:2] = x[1:-1]
    xc[2:-1:2] = x[1:-1]
    xc[0]      = x[0]
    xc[-1]     = x[-1]
    yc = np.zeros((y.size*2,))
    yc[::2]    = y
    yc[1::2]   = y

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_xlim(x[0], x[-1])
    if logscale:
        axes.set_yscale('log')
    axes.plot(xc, yc,label=label)

@openfigureifnoaxes(projection=None)
def mesh2d(x, y, z, logscale=False, diverging=False, xlabel=None, ylabel=None,
           clabel=None, clim=None, cmap=None, axesequal=False,
           axes=None, cax=None):
    """Make a mesh (surface) plot in 2D.

    Parameters
    ----------
    x : array_like (nx,) or (nx+1,)
        Abscissa or abscissa edges for the x-axis.
    y : array_like (ny,) or (ny+1,)
        Abscissa or abscissa edges for the y-axis.
    z : array_like (nx,ny)
        Data to be plotted.
    logscale : bool, optional
        Make color axis logarithmic.
    diverging : bool, optional
        Use diverging colormapping which is centered at zero.

        Works with logarithmic scale as well.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color-axis.
    clim : [float, float], optional
        Color [min, max] limits.
    cmap : str, optional
        Colormap.
    axesequal : bool, optional
        Flag to set the aspect ratio equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    if clim is None: clim = [None, None]
    if clim[0] is None:
        clim[0] = np.nanmin(z)
    if clim[1] is None:
        clim[1] = np.nanmax(z)
    z = np.ma.masked_invalid(z)

    if logscale:
        if diverging:
            if cmap == None: cmap = "bwr"
            norm = mpl.colors.SymLogNorm(linthresh=10, linscale=1.0,
                                 vmin=clim[0], vmax=clim[1], base=10)
        else:
            if cmap == None: cmap = "viridis"
            if clim[0] <=0: clim[0] = np.min(z[z>0])
            norm = mpl.colors.LogNorm(vmin=clim[0], vmax=clim[1])
    else:
        if diverging:
            if cmap == None: cmap = "bwr"
            norm = mpl.colors.CenteredNorm(halfrange=np.amax(np.abs(clim)))
        else:
            if cmap == None: cmap = "viridis"
            norm = mpl.colors.Normalize(vmin=clim[0], vmax=clim[1])

    shading = "nearest" if z.shape == (x.size,y.size) else "flat"
    pcm = axes.pcolormesh(x, y, z.T, norm=norm, cmap=cmap, shading=shading)
    axes.patch.set(hatch='x', edgecolor=[0.9, 0.9, 0.9])

    axes.set_xlim(x[0], x[-1])
    axes.set_ylim(y[0], y[-1])

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    if axesequal:
        axes.set_aspect("equal", adjustable="box")

    cbar = plt.colorbar(pcm, ax=axes, cax=cax)
    cbar.set_label(clabel)

@openfigureifnoaxes(projection=None)
def contour2d(x, y, z, contours=None, xlabel=None, ylabel=None, axesequal=False,
              colors=None, linestyles=None, linewidths=None, axes=None,
              contourlabels=False, labelfontsize=10):
    """Plot contour on 2D plane.
    """
    if contours==None:
        # The level argument needs to be left out if the contour levels are not
        # specified
        cs = axes.contour(x, y, z.T, colors=colors, linestyles=linestyles,
                 linewidths=linewidths)
    else:
        cs = axes.contour(x, y, z.T, contours, colors=colors, linestyles=linestyles,
                 linewidths=linewidths)
    if contourlabels: axes.clabel(cs, inline=True, fontsize=labelfontsize)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if axesequal:
        axes.set_aspect("equal", adjustable="box")

@openfigureifnoaxes(projection=None)
def arrows2d():
    """Plot vector field on 2D plane.
    """
    pass

@openfigureifnoaxes(projection=None)
def line2d(x, y, c=None, xlog="linear", ylog="linear", clog="linear",
           xlabel=None, ylabel=None, clabel=None, bbox=None,
           cmap=None, axesequal=False, axes=None, cax=None, marker=None,
           markerfacecolor=None, title=None,label=None, skipshow=False):
    """Plot line segments in 2D.

    x : array_like
        Marker x-coordinates.
    y : array_like
        Marker y-coordinates.
    c : str or array_like, optional
        Color data or string indicating the color.
    xlog : {"linear", "log", "symlog"}, optional
        x-axis scaling.
    ylog : {"linear", "log", "symlog"}, optional
        y-axis scaling.
    clog : {"linear", "log", "symlog"}, optional
        color-axis scaling.
    cmap : str, optional
        Name of the colormap.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    axesequal : bool, optional
        Flag to set aspect ratio of [x,y] axes equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    marker : str, optional
        Marker type.
    markerfacecolor : str, optional
        Color of the marker face.
    title : str, optional
        Title of the figure.
    label : str, optional
        Label for the legend.
    skipshow : bool, optional
        Flag to skip the show function. Only used by the openfigureifnoaxes decorator.
    """
    axes.set_xscale(xlog)
    axes.set_yscale(ylog)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title)
    if axesequal: axes.set_aspect("equal", adjustable="box")
    if c is None or isinstance(c, str):
        # Simple plot with a single color
        for i in range(len(x)):
            if i == 0:
                axes.plot(x[i], y[i], color=c, marker=marker,
                          markerfacecolor=markerfacecolor,label=label)
            else:
                axes.plot(x[i], y[i], color=c, marker=marker,
                          markerfacecolor=markerfacecolor)
        return axes

    if bbox is None:
        bbox = "[xmin, xmax, ymin, ymax, cmin, cmax]"
        raise ValueError("Specify bounding box: bbox="+bbox)
    axes.set_xlim(bbox[0], bbox[1])
    axes.set_ylim(bbox[2], bbox[3])
    if clog == "linear":
        if bbox[-2]*bbox[-1] < 0:
            if cmap == None: cmap = "bwr"
            norm = mpl.colors.CenteredNorm(halfrange=np.amax(np.abs(bbox[-2:])))
        else:
            if cmap == None: cmap = "viridis"
            norm = plt.Normalize(bbox[-2], bbox[-1])
    elif clog == "log":
        if cmap == None: cmap = "viridis"
        norm = mpl.colors.LogNorm(vmin=bbox[-2], vmax=bbox[-1])
    elif clog == "symlog":
        if cmap == None: cmap = "bwr"
        norm = mpl.colors.SymLogNorm(linthresh=10, linscale=1.0,
                                     vmin=bbox[-2], vmax=bbox[-1], base=10)

    for i in range(len(x)):
        points = np.array([x[i], y[i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = mpl.collections.LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(c[i][1:])
        line = axes.add_collection(lc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = plt.colorbar(smap, ax=axes, cax=cax)
    cax.ax.set_ylabel(clabel, rotation=90, labelpad=10)
    return axes

@openfigureifnoaxes(projection="3d")
def line3d(x, y, z, c=None, xlog="linear", ylog="linear", zlog="linear",
           clog="linear", xlabel=None, ylabel=None, zlabel=None, clabel=None,
           bbox=None, cmap=None, axesequal=False, axes=None, cax=None):
    """Plot line segments in 3D.

    Parameters
    ----------
    x : array_like
        Marker x-coordinates.
    y : array_like
        Marker y-coordinates.
    z : array_like
        Marker z-coordinates.
    c : str or array_like, optional
        Color data or string indicating the color.
    xlog : {"linear", "log", "symlog"}, optional
        x-axis scaling.
    ylog : {"linear", "log", "symlog"}, optional
        y-axis scaling.
    zlog : {"linear", "log", "symlog"}, optional
        z-axis scaling.
    clog : {"linear", "log", "symlog"}, optional
        color-axis scaling.
    cmap : str, optional
        Name of the colormap.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    zlabel : str, optional
        Label for the z-axis.
    clabel : str, optional
        Label for the color axis.
    axesequal : bool, optional
        Flag to set aspect ratio of [x,y,z] axes equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    axes.set_xscale(xlog)
    axes.set_yscale(ylog)
    axes.set_zscale(zlog)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_zlabel(zlabel)
    if axesequal: axes.set_aspect("equal", adjustable="box")
    if c is None or isinstance(c, str):
        # Simple plot with a single color
        for i in range(len(x)):
            axes.plot(x[i], y[i], z[i], color=c)
        return

    if bbox is None:
        bbox = "[xmin, xmax, ymin, ymax, cmin, cmax]"
        raise ValueError("Specify bounding box: bbox="+bbox)
    axes.set_xlim(bbox[0], bbox[1])
    axes.set_ylim(bbox[2], bbox[3])
    axes.set_zlim(bbox[4], bbox[5])
    if clog == "linear":
        if bbox[-2]*bbox[-1] < 0:
            if cmap == None: cmap = "bwr"
            norm = mpl.colors.CenteredNorm(halfrange=np.amax(np.abs(bbox[-2:])))
        else:
            if cmap == None: cmap = "viridis"
            norm = plt.Normalize(bbox[-2], bbox[-1])
    elif clog == "log":
        if cmap == None: cmap = "viridis"
        norm = mpl.colors.LogNorm(vmin=bbox[-2], vmax=bbox[-1])
    elif clog == "symlog":
        if cmap == None: cmap = "bwr"
        norm = mpl.colors.SymLogNorm(linthresh=10, linscale=1.0,
                                     vmin=bbox[-2], vmax=bbox[-1], base=10)
    for i in range(len(x)):
        points = np.array([x[i], y[i], z[i]]).T.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = mpl_toolkits.mplot3d.art3d.Line3DCollection(
            segments, cmap=cmap, norm=norm)
        lc.set_array(c[i][1:])
        line = axes.add_collection(lc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    if clabel is not None:
        cbar.set_label(clabel)

@openfigureifnoaxes(projection="3d")
def surface3d(x_grid1d, y_grid1d, z_grid1d, qnt_grid1d=None, diverging=False,
              logscale=False, axesequal=False, xlabel=None, ylabel=None,
              zlabel=None, clabel=None, clim=None, cmap=None, axes=None,
              cax=None, meshalpha=0.6, tri_lc="gray", tri_lw=0.3, opacity=0.8,
              surfacecolor="purple"):
    """
    Make a 3D surface plot from 1D coordinate arrays and optional quantity values.

    You may encounter "WARNING:root:Singular matrix. Likely caused by all points
    lying in an N-1 space." or two but, as the name suggests, it is just a
    warning.

    Parameters
    ----------
    x_grid1d : array_like (n,)
        X-coordinates of the surface points.
    y_grid1d : array_like (n,)
        Y-coordinates of the surface points.
    z_grid1d : array_like (n,)
        Z-coordinates of the surface points.
    qnt_grid1d : array_like (n,), optional
        Quantity values associated with each point, used for coloring the surface.
    diverging : bool, optional
        Use a diverging colormap centered at zero.
    logscale : bool, optional
        Apply logarithmic scaling to the colormap.
    axesequal : bool, optional
        Set 3D axes to have equal scaling.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    zlabel : str, optional
        Label for the z-axis.
    clabel : str, optional
        Label for the colorbar.
    clim : [float, float], optional
        Limits for the colormap [min, max].
    cmap : str or Colormap, optional
        The colormap to use for surface coloring.
    axes : :obj:`~matplotlib.axes._subplots.Axes3D`, optional
        The 3D axes to plot on, or None to create a new figure and axes.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        Axes to draw the colorbar on, or None to place it next to the main axes.
    meshalpha : float, optional
        A parameter used when generating the triangles (google: alpha shape).
    tri_lc : str or color, optional
        Color of triangle edges.
    tri_lw : float, optional
        Line width of triangle edges.
    opacity : float, optional
        Overall opacity of the surface (0: transparent, 1: opaque).
    surfacecolor : str, optional
        If no qnt is given, this color is used. By default purple because
        magneticfield itself is purple, everyone knows that.
    """

    # Put co-ordinate triplets into 2d array and throw away nans
    points = np.vstack((x_grid1d, y_grid1d, z_grid1d)).T
    mask = ~np.isnan(points).any(axis=1)
    points = points[mask]
    if qnt_grid1d is not None: quantity = qnt_grid1d[mask]

    # Create the triangles
    alpha_shape = alphashape.alphashape(points, meshalpha)
    triangles = np.array(list(alpha_shape.faces))
    vertices = np.array(alpha_shape.vertices)

    # Get the surface color from the desired quantity
    if qnt_grid1d is not None:
        tree = cKDTree(points)
        _, idx = tree.query(vertices, k=1)
        vertex_quantity = quantity[idx]
        triangle_colors = vertex_quantity[triangles].mean(axis=1) #average of vertex

        if clim is None: clim = [None, None]
        if clim[0] is None:
            clim[0] = np.nanmin(triangle_colors)
        if clim[1] is None:
            clim[1] = np.nanmax(triangle_colors)

        if logscale:
            if diverging:
                if cmap == None: cmap = "bwr"
                norm = mpl.colors.SymLogNorm(linthresh=10, linscale=1.0,
                                    vmin=clim[0], vmax=clim[1], base=10)
            else:
                if cmap == None: cmap = "viridis"
                if clim[0] <=0: clim[0] = np.min(quantity[quantity>0])
                norm = mpl.colors.LogNorm(vmin=clim[0], vmax=clim[1])
        else:
            if diverging:
                if cmap == None: cmap = "bwr"
                norm = mpl.colors.CenteredNorm(halfrange=np.amax(np.abs(clim)))
            else:
                if cmap == None: cmap = "viridis"
                norm = mpl.colors.Normalize(vmin=clim[0], vmax=clim[1])

        colors = mpl.cm.viridis(norm(triangle_colors))
    else:
        # No qnt given for coloring => plot just the flux surface
        colors = surfacecolor

    # Plotting
    mesh = mpl_toolkits.mplot3d.art3d.Poly3DCollection(
        vertices[triangles], facecolors=colors, edgecolor=tri_lc, lw=tri_lw,
        alpha=opacity)

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111, projection='3d')
    axes.add_collection3d(mesh)

    axes.set_xlim(vertices[:,0].min(), vertices[:,0].max())
    axes.set_ylim(vertices[:,1].min(), vertices[:,1].max())
    axes.set_zlim(vertices[:,2].min(), vertices[:,2].max())
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_zlabel(zlabel)

    if axesequal:
        axes.set_aspect("equal", adjustable="box")

    if qnt_grid1d is not None:
        plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                 cax=cax, shrink=0.5, label=clabel)

@openfigureifnoaxes(projection=None)
def poincare(x, y, ids, connlen=None, xlim=None, ylim=None, xlabel=None,
             ylabel=None, clabel=None, axesequal=False, markersize=2,
             axes=None, cax=None):
    """Poincaré plot where color separates markers or shows the connection
    length.

    For better performance this function uses matplotlib's plot function instead
    of the scatter function. The only practical difference is that the marker
    color can't be continuous.

    Parameters
    ----------
    x : array_like
        Orbit x-coordinates.
    y : array_like
        Orbit y-coordinates.
    ids : array_like
        Array of marker IDs showing to which marker the points in x and y
        arrays correspond to.
    connlen : array_like, optional
        Connection length at the position (x,y).

        Negative if the marker is confined. If given, the color scale shows
        connection length instead of marker ID. The confined markers are still
        shown with shades of red.
    xlim : tuple(float, float), optional
        Min and max values for the x-axis.
    ylim : tuple(float, float), optional
        Min and max values for the y-axis.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    axesequal : bool, optional
        If True, x and y axis have equal aspect ratio.
    markersize : int, optional
        Marker size on plot.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    nc = 5 # How many colors are used

    if connlen is not None:
        # If we don't have lost markers then ignore the connection length.
        if np.argwhere(connlen > 0).size == 0:
            connlen = None

    if connlen is None:
        # Shuffle colors and plot
        cmap = plt.cm.get_cmap("Reds", nc)
        uids  = np.unique(ids)
        cpick = np.arange(nc)
        np.random.shuffle(cpick)
        for i in range(nc):
            idx = np.in1d(ids, uids[i::nc])
            axes.plot(x[idx], y[idx], color=cmap(cpick[i]/nc),
                      linestyle="None", marker=".", markersize=markersize)

    else:
        # Sort by connection length (confined markers are indicated with a
        # negative connlen).
        idx = np.argsort(connlen)
        x   = x[idx]
        y   = y[idx]
        ids = ids[idx]
        connlen = connlen[idx]

        # Find where the line between confined and lost markers is.
        # Set connlen positive again for confined markers and rearrange again
        # by the connection length.
        idx    = np.argwhere(connlen > 0).ravel()[0]
        x      = np.append(x[idx:],      np.flip(x[:idx]))
        y      = np.append(y[idx:],      np.flip(y[:idx]))
        ids    = np.append(ids[idx:],    np.flip(ids[:idx]))
        connlen = np.append(connlen[idx:], np.flip(-connlen[:idx]))
        idx = (connlen.size - 1) - idx

        # The color has meaning only for lost markers so find the scale
        cmin = connlen[0]
        cmax = connlen[idx]

        logscale = False
        if cmin / cmax < 0.1:
            logscale = True
            cmin   = np.log10(cmin)
            cmax   = np.log10(cmax)
            connlen = np.log10(connlen)
            cunits = 1
        else:
            cunits = cmax.units

        cmin = np.floor(cmin)
        cmax = np.ceil(cmax)
        nc_b = int(cmax - cmin)
        clim = np.linspace(cmin, cmax, nc_b+1)

        # Confined markers are plotted separately.
        connlen[idx+1:] = cmax + (1/nc) * cunits
        clim = np.append(clim, cmax + np.linspace(0, 1/nc, nc) * cunits)

        # Create colormap and colorbar
        cmapred  = plt.cm.get_cmap("Reds_r")
        cmapblue = plt.cm.get_cmap("Blues")
        colours  = [None] * (nc_b + nc)
        for i in range(nc_b):
            colours[i]  = cmapblue(i/nc_b)
        for i in range(nc):
            colours[i+nc_b] = cmapred(i/nc)

        cmap = mpl.colors.ListedColormap(colours)
        norm = mpl.colors.BoundaryNorm(boundaries=clim, ncolors=nc_b+nc )
        ticks = clim[:nc_b+1]
        if logscale:
            ticklabels = list(clim[:nc_b]) + [r"$\inf$"]
        else:
            ticklabels = list(clim[:nc_b]) + [r"$\inf$"]
        deltac = clim[1] - clim[2]
        cb = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
            spacing='proportional', extend='min', ticks=ticks,
            boundaries=[clim[0]-deltac/2] + clim + [clim[-1]+deltac/2])
        cb.ax.set_yticklabels(ticklabels)
        if clabel is not None:
            clabel = "log10( "+clabel+" )" if logscale else clabel
            cb.ax.set_ylabel(clabel)

        # Find the indices that divide the color range in even intervals
        idx  = np.searchsorted(connlen, clim) + 1

        # Now plot markers on to each interval with a different color
        i1 = 0
        for i in range(nc_b):
            i2 = idx[i+1]
            axes.plot(x[i1:i2], y[i1:i2], color=colours[i],
                      linestyle="None", marker="o", markersize=markersize)
            i1 = i2

        # Plot confined
        uids  = np.unique(ids[i2:])
        cpick = np.arange(nc)
        np.random.shuffle(cpick)
        for i in range(nc):
            idx = np.in1d(ids, uids[i::nc])
            axes.plot(x[idx], y[idx], color=colours[nc_b+cpick[i]],
                      linestyle="None", marker=".", markersize=markersize)

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    if axesequal:
        axes.set_aspect("equal", adjustable="box")

@openfigureifnoaxes(projection=None)
def still(wallmesh, points=None, orbit=None, data=None, log=False, clim=None,
          cpos=None, cfoc=None, cang=None, axes=None, cax=None, **kwargs):
    """Take a still shot of the mesh and display it using matplotlib backend.

    The rendering is done using vtk but the vtk (interactive) window is not
    displayed. It is recommended to use the interactive plot to find desired
    camera position and produce the actual plot using this method. The plot
    is shown using imshow and acts as a regular matplotlib plot.

    Parameters
    ----------
    wallmesh : Polydata
        Mesh representing the wall.
    points : array_like, optional
        Array Npoint x 3 defining points (markers) to be shown. For
        each point [x, y, z] coordinates are given.
    orbit : array_like, (n,3), optional
        Cartesian coordinates for an orbit to be plotted.
    data : str, optional
        Name of the cell data in the wall mesh that is shown in color.
    log : bool, optional
        Color range is logarithmic if True.
    clim : [float, float], optional
        Color [min, max] limits.
    cpos : array_like, optional
        Camera position coordinates [x, y, z].
    cfoc : array_like, optional
        Camera focal point coordinates [x, y, z].
    cang : array_like, optional
        Camera angle [azimuth, elevation, roll].
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    **kwargs
        Keyword arguments passed to :obj:`~pyvista.Plotter`.
    """
    p = pv.Plotter(off_screen=True, **kwargs)
    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9])
        clim = None
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9,0.9,0.9])
        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])
        if clim is None: clim = [minval, maxval]
        if clim[0] is None: clim[0] = minval
        if clim[1] is None: clim[1] = maxval

        p.add_mesh(wallmesh, scalars=data, cmap=cmap, clim=clim, log_scale=log)
        p.remove_scalar_bar()

    if points is not None:
        p.theme.color = 'black'
        p.add_points(points, render_points_as_spheres=True, point_size=10)

    if orbit is not None:
        orbit = pv.lines_from_points(orbit)
        p.add_mesh(orbit, color="red")

    # Set camera
    if cpos is not None:
        p.camera.position = cpos
    if cfoc is not None:
        p.camera.focal_point = cfoc
    if cang is not None:
        p.camera.azimuth   = cang[0]
        p.camera.elevation = cang[1]
        p.camera.roll      = cang[2]

    p.show()
    axes.imshow(p.image)
    axes.set_xticks([])
    axes.set_yticks([])

    if data is not None:
        if log:
            norm = mpl.colors.LogNorm(vmin=clim[0], vmax=clim[1])
        else:
            norm = mpl.colors.Normalize(vmin=clim[0], vmax=clim[1])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, ax=axes, cax=cax)

        if data == "eload":
            cbar.set_label(r"Power load W/m$^2$")

def interactive(wallmesh,
                *args,
                points=None,
                orbit=None,
                data=None,
                log=False,
                clim=None,
                cmap=None,
                cbar_title=None,
                cpos=None,
                cfoc=None,
                cang=None,
                p=None,
                phi_lines=None,
                const_phi_planes=None,
                theta_lines=None,
                a5=None,
                skipshow=False,
                **kwargs):
    """Open VTK window to display interactive view of the wall mesh.

    Parameters
    ----------
    wallmesh : Polydata
        Mesh representing the wall.
    *args : tuple (str, method), optional
        Key (str) method pairs. When key is pressed when the plot is
        displayed, the associated method is called. The method should take
        Plotter instance as an argument.
    points : array_like, optional
        Array Npoint x 3 defining points (markers) to be shown. For
        each point [x, y, z] coordinates are given.
    orbit : array_like, (n,3), optional
        Cartesian coordinates for an orbit to be plotted.
    data : str, optional
        Name of the cell data in the wall mesh that is shown in color.
    log : bool, optional
        Color range is logarithmic if True.
    clim : [float, float], optional
        Color [min, max] limits.
    cmap : str, optional
        Colormap name.
    cbar_title : str, optional
        Color bar title.
    cpos : array_like, optional
        Camera position coordinates [x, y, z].
    cfoc : array_like, optional
        Camera focal point coordinates [x, y, z].
    cang : array_like, optional
        Camera angle [azimuth, elevation, roll].
    p : :obj:`~pyvista.Plotter`, optional
        PyVista plotter instance.
    phi_lines : array_like, optional
        Array of phi values in degrees.
    const_phi_planes : array_like, optional
        Array of phi values in degrees.
    theta_lines : array_like, optional
        Array of theta values in degrees.
    a5 : a5py.ASCOT, optional
        ASCOT instance.
    skipshow : bool, optional
        If True, do not show the plot. Useful if you want to draw something else
        like a flux surface in the same Plotter.
    **kwargs
        Keyword arguments passed to :obj:`~pyvista.Plotter`.
    """
    if p is None:
        p = pv.Plotter(**kwargs)
    p.disable()
    cameracontrols(p)

    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9], log_scale=log)
    else:
        if cmap is None:
            cmap = mpl.colormaps["Reds"].copy()
            cmap.set_bad(color=[0.9,0.9,0.9])
        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])
        if clim is None: clim = [minval, maxval]
        if clim[0] is None: clim[0] = minval
        if clim[1] is None: clim[1] = maxval
        if cbar_title is None: cbar_title = data
        p.add_mesh(wallmesh,
                   scalars=data,
                   cmap=cmap,
                   clim=clim,
                   log_scale=log,
                   scalar_bar_args={
                        "title": cbar_title,
                        "title_font_size": 22,
                        "label_font_size": 18,
                    },
        )

    if points is not None:
        p.theme.color = 'black'
        p.add_points(points, render_points_as_spheres=True, point_size=10)

    if orbit is not None:
        orbit = pv.lines_from_points(orbit)
        p.add_mesh(orbit, color="red")

    # Set events
    for i in range(len(args)):
        p.add_key_event(args[i][0], lambda : args[i][1](p))

    # Set camera
    if cpos is not None:
        p.camera.position = cpos
    if cfoc is not None:
        p.camera.focal_point = cfoc
    if cang is not None:
        p.camera.azimuth   = cang[0]
        p.camera.elevation = cang[1]
        p.camera.roll      = cang[2]

    # Constant phi lines with z=0, extending radially outward
    if phi_lines is not None:
        add_phi_lines(wallmesh, p, phi_lines, a5)

    # Constant phi planes with constant theta lines on them
    if const_phi_planes is not None:
        add_phi_planes(wallmesh=wallmesh,
                    plotter=p,
                    phi_planes=const_phi_planes,
                    )
        if theta_lines is not None:
            add_theta_lines(p,
                            const_phi_planes,
                            theta_lines,
                            a5,
                            rminor_wall=4*unyt.m)

    p.renderer.ResetCameraClippingRange()
    if not skipshow:
        p.show()

def add_highlighted_edges(plotter,
                          wallmesh_highlight,
                          color="yellow",
                          line_width=3):
    """
    Highlight the edges of a wall mesh on an existing PyVista plotter.

    Parameters
    ----------
    plotter : pv.Plotter
        The PyVista plotter instance where the mesh is already displayed.
    wallmesh_highlight : pv.PolyData
        The wall mesh whose edges you want to highlight.
    color : str, optional
        Color of the edges (default "yellow").
    line_width : int, optional
        Thickness of the edge lines (default 3).
    """
    plotter.add_mesh(
        wallmesh_highlight,
        style="wireframe",
        color=color,
        line_width=line_width
    )

@physlib.parseunits(phi_lines="deg")
def add_phi_lines(wallmesh,
                  plotter,
                  phi_lines,
                  a5,
                  ):
    """
    Add phi lines to a PyVista plotter.

    Parameters
    ----------
    wallmesh : pv.PolyData
        The wall mesh to which the phi lines will be added.
    plotter : pv.Plotter
        The PyVista plotter instance where the mesh is already displayed.
    phi_lines : np.ndarray
        Array of phi values.
    a5 : a5py.ASCOT
        An instance of the ASCOT class.
    """
    bounds = wallmesh.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
    xmax = max(abs(bounds[0]), abs(bounds[1]))
    ymax = max(abs(bounds[2]), abs(bounds[3]))
    max_r = np.sqrt(xmax**2 + ymax**2)
    x = max_r * np.cos(phi_lines)
    y = max_r * np.sin(phi_lines)
    out = a5._eval_bfield(1*unyt.m, phi_lines.to("rad"), 1*unyt.m, 0, evalaxis=True)
    x_axis = out["axisr"]*np.cos(phi_lines)
    y_axis = out["axisr"]*np.sin(phi_lines)

    for i in range(len(phi_lines)):
        # line from origin to outer radius in XY plane, z=0
        line = pv.Line([0.05*x[i], 0.05*y[i], 0], [x[i], y[i], 0])
        plotter.add_mesh(line,
                         color="green",
                         opacity=1.0,
                         line_width=2,
                         pickable=False,
                         reset_camera=False,
                         )
        labels = plotter.add_point_labels([[x_axis[i], y_axis[i], 0.05]],
                                          [f"phi = {phi_lines[i].to('deg').to_value():.0f}°"],
                                          point_size=10,
                                          font_size=22,
                                          text_color="green",
                                          shape_opacity=0.0,
                                          show_points=False,
                                          )
        labels.GetProperty().SetDisplayLocationToForeground()

@physlib.parseunits(phi_planes="deg")
def add_phi_planes(wallmesh,
                   plotter,
                   phi_planes,
                   ):
    """
    Add phi planes to a PyVista plotter.

    Parameters
    ----------
    wallmesh : pv.PolyData
        The wall mesh to which the phi planes will be added.
    plotter : pv.Plotter
        The PyVista plotter instance where the mesh is already displayed.
    phi_planes : np.ndarray
        Array of phi values.
    """
    bounds = wallmesh.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
    xmax = max(abs(bounds[0]), abs(bounds[1]))
    ymax = max(abs(bounds[2]), abs(bounds[3]))
    zmin, zmax = bounds[4], bounds[5]
    max_r = np.sqrt(xmax**2 + ymax**2)
    phi = 45*unyt.deg
    for phi in phi_planes:
        # direction vector in XY plane
        dx, dy = np.cos(phi), np.sin(phi)
        # define the four corners of the rectangle (plane)
        corners = np.array([
            [0, 0, zmin],
            [0, 0, zmax],
            [max_r*dx, max_r*dy, zmax],
            [max_r*dx, max_r*dy, zmin],
        ])
        # make a quad surface (polygon)
        faces = np.hstack([[4, 0, 1, 2, 3]])  # 4-point face
        plane = pv.PolyData(corners, faces)
        # add with some transparency
        actor = plotter.add_mesh(plane, color="blue", opacity=0.2)
        #actor.SetForceTranslucent(True)

@physlib.parseunits(theta_lines="deg")
def add_theta_lines(plotter,
                    phi_array,
                    theta_lines,
                    a5,
                    rminor_wall,
                    ):
    """
    Add theta lines to a PyVista plotter.

    Parameters
    ----------
    plotter : pv.Plotter
        The PyVista plotter instance where the mesh is already displayed.
    phi_array : np.ndarray
        Array of phi values.
    theta_lines : np.ndarray
        Array of theta values.
    a5 : a5py.ASCOT
        An instance of the ASCOT class.
    rminor_wall : float
        Minor radius of the wall. Used to determine the theta label positions.
    """
    out = a5._eval_bfield(1*unyt.m, phi_array.to("rad"), 1*unyt.m, 0, evalaxis=True)
    axis_x = out["axisr"]*np.cos(phi_array)
    axis_y = out["axisr"]*np.sin(phi_array)
    axis_z = out["axisz"]
    for i in range(len(phi_array)):
        phi = phi_array[i]
        for theta in theta_lines:
            d = rminor_wall

            # line beginning co-ordinates close to magnetic axis
            x0 = axis_x[i] + 0.05*d*np.cos(theta)*np.cos(phi)
            y0 = axis_y[i] + 0.05*d*np.cos(theta)*np.sin(phi)
            z0 = axis_z[i] + 0.05*d*np.sin(theta)

            dr1 = d * np.cos(theta)
            dz1 = d * np.sin(theta)

            # line end co-ordinates outside the wall
            R1 = out["axisr"][i]+dr1
            x1 = R1 * np.cos(phi)
            y1 = R1 * np.sin(phi)
            z1 = out["axisz"][i]+dz1

            dr2 = dr1/2
            dz2 = dz1/2

            # Co-ordinates for labels, half way between rminor=0 and wall
            R2 = out["axisr"][i]+dr2
            x2 = R2 * np.cos(phi)
            y2 = R2 * np.sin(phi)
            z2 = out["axisz"][i]+dz2

            line = pv.Line([x0, y0, z0], [x1, y1, z1])
            plotter.add_mesh(line,
                             color="blue",
                             opacity=1.0,
                             line_width=2,
                             pickable=False,
                             reset_camera=False,
                             )
            labels = plotter.add_point_labels([[x2, y2, z2]],
                                              [f"{theta:.0f}°"],
                                              point_size=10,
                                              font_size=22,
                                              text_color="blue",
                                              shape_opacity=0.0,
                                              )
            labels.GetProperty().SetDisplayLocationToForeground()



@openfigureifnoaxes(projection=None)
def loadvsarea(wetted, loads, axes=None):
    """Plot histogram showing minimum load vs area.
    """
    idx = np.argsort(-loads)
    wetted = np.cumsum(wetted[idx])
    loads  = loads[idx]

    axes.set_xscale('log')
    axes.set_yscale('linear')
    axes.plot(loads, wetted)
    axes.set_xlabel(r"Load above [" + str(loads.units) + "]")
    axes.set_ylabel(r"Wetted area [" + str(wetted.units) + "]")

@openfigureifnoaxes(projection=None)
def triangularpatch(
        patches, color, log=False, xlim=None, ylim=None, clim=None, xlabel=None,
        ylabel=None, clabel=None, cmap=None, axes=None, cax=None):
    """Plot triangular patches.
    """
    if clim    is None: clim = [None, None]
    if clim[0] is None: clim[0] = np.nanmin(color)
    if clim[1] is None: clim[1] = np.nanmax(color)

    if log:
        norm = mpl.colors.LogNorm(clim[0], clim[1])
    else:
        norm = mpl.colors.Normalize(clim[0], clim[1])

    coll = mpl.collections.PolyCollection(patches, array=color, norm=norm,
                                          cmap=cmap)
    axes.add_collection(coll)

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    plt.colorbar(coll, norm=norm, ax=axes, cax=cax, label=clabel)

@openfigureifnoaxes(projection="polar")
def momentumpolarplot(pnorm_edges, pitch_edges, dist, axes=None, cax=None):
    """Plot momentum space distribution in polar coordinates.

    Parameters
    ----------
    pnorm_edges : array_like
        Momentum abscissa edges.
    pitch_edges : array_like
        Pitch abscissa edges.
    dist : array_like
        Values of the distribution.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    theta = np.arccos(pitch_edges)
    h = axes.pcolormesh(theta, pnorm_edges, dist)
    plt.colorbar(h, ax=axes, cax=cax)
    axes.set_thetamin(0)
    axes.set_thetamax(180)
    axes.set_xticks(np.array([0, 45, 90, 135, 180])*np.pi/180)
    axes.set_xticklabels([1.0, 0.5, 0.0, -0.5, -1.0])

@openfigureifnoaxes(projection=None)
def momentumpolargrid(pnorm_edges, pitch_edges, axes=None):
    """Plot momentum space polar coordinate grid in Cartesian basis.

    Parameters
    ----------
    pnorm_edges : array_like
        Momentum abscissa edges.
    pitch_edges : array_like
        Pitch abscissa edges.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    """
    ang = np.linspace(0, np.pi, 60)
    for v in pnorm_edges:
        axes.plot(v*np.cos(ang), v*np.sin(ang), color="black")

    p = pnorm_edges[-1]
    for v in pitch_edges:
        axes.plot([0, v*p], [0, np.sqrt(1.0 - v**2)*p], color="black")

@openfigureifnoaxes(projection=None)
def radialprofile(x, y1, y2=None, xlim=None, y1lim=None, y2lim=None,
                  xlabel=None, y1label=None, y2label=None, y1legends=None,
                  y2legends=None, axes=None, tightlayout=True, title=None):
    """Plot 1D profiles on axes that can have two y-axes and the y-axis combines
    both linear and logarithmic scale.

    Parameters
    ----------
    x : array_like, (n,)
        The x grid where ``y1`` (and ``y2``) values are provided.
    y1 : array_like or [array_like]
        The values (or a list of values in which case they are separated by
        colour) plotted on the left y-axis.
    y2 : array_like or [array_like], optional
        The values (or a list of values) plotted on the right y-axis.
    xlim : [float, float], optional
        Limits on x-axis.
    y1lim : [float, float, float], optional
        Limits on the first y axis where the middle value is when the scale
        changes from logarithmic to linear.
    y2lim : [float, float, float], optional
        Limits on the second y axis where the middle value is when the scale
        changes from logarithmic to linear.
    xlabel : str, optional
        Label on the x-axis.
    y1label : str, optional
        Label on the first y-axis.
    y2label : str, optional
        Label on the second y-axis.
    y1legends: [str], optional
        Legends for the values plotted on the first y-axis.

        Number of legend values must be the same as the number of ``y1``.
    y2legends: [str], optional
        Legends for the values plotted on the second y-axis.

        Number of legend values must be the same as the number of ``y2``.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    tightlayout : bool, optional
        Whether to use tight layout.
    title : str, optional
        Title of the figure.
    """
    if y1lim is None: raise ValueError("y1lim must be provided")
    if xlim is not None: axes.set_xlim(xlim)
    axes.set_title(title)

    # Create linear left axis
    axleftlin = axes
    axleftlin.set_yscale('linear')
    axleftlin.spines['right'].set_visible(False)
    axleftlin.spines['bottom'].set_visible(False)
    # Create log left axis
    divider = make_axes_locatable(axes)
    axleftlog = divider.append_axes('bottom', size=1.0, pad=0, sharex=axes)
    axleftlog.set_yscale('log')
    axleftlog.spines['right'].set_visible(False)
    axleftlog.spines['top'].set_visible(False)
    axleftlog.yaxis.set_ticks_position('left')

    axleftlog.tick_params(axis='y', which='minor', left=False)
    axleftlin.yaxis.set_major_formatter(getmathtextsciformatter("%1.1e"))

    axleftlog.set_xlabel(xlabel)
    plt.setp(axleftlin.get_xticklabels(), visible=False)

    if y2 is not None:
        # Create linear right axis
        axrightlin = axes.twinx()
        axrightlin.set_yscale('linear')
        axrightlin.spines['left'].set_visible(False)
        axrightlin.spines['bottom'].set_visible(False)

        # Create log right axis
        divider = make_axes_locatable(axrightlin)
        axrightlog = divider.append_axes("bottom", size=1.0, pad=0, sharex=axes)
        axrightlog.set_yscale('log')
        axrightlog.spines['left'].set_visible(False)
        axrightlog.spines['top'].set_visible(False)
        axrightlog.yaxis.set_ticks_position('right')
        axrightlog.xaxis.set_visible(False)
        axrightlog.set_facecolor('none')
        axrightlog.yaxis.set_label_position("right")

        axrightlog.tick_params(axis='y', which='minor', right=False)
        axrightlin.yaxis.set_major_formatter(getmathtextsciformatter("%1.1e"))

        axlin      = axrightlin
        axrightlin = axleftlin
        axleftlin  = axlin
        axlog      = axrightlog
        axrightlog = axleftlog
        axleftlog  = axlog

        axrightlin.set_ylim((y2lim[1], y2lim[2]))
        axrightlog.set_ylim((y2lim[0], y2lim[1]))
        axrightlin.set_ylabel(y2label, loc='bottom')

        axrightlin.spines['right'].set_color('C3')
        axrightlog.spines['right'].set_color('C3')
        axrightlin.tick_params(axis='y', colors='C3')
        axrightlog.tick_params(axis='y', colors='C3')
        axrightlin.yaxis.label.set_color('C3')

        if not isinstance(y2, list): y2 = [y2]
        handles2 = []
        for i, y in enumerate(y2):
            ls = '-' if i == 0 else '--'
            h, = axrightlin.plot(x, y, color='C3', ls=ls)
            axrightlog.plot(x, y, color='C3', ls=ls)
            handles2.append(h)
        axrightlog.set_xlim(xlim)

        legend2 = plt.legend(handles2, y2legends, loc='upper left',
                             bbox_to_anchor=(0.0,1.0), frameon=False)
        axleftlog.add_artist(legend2)

    axleftlin.set_ylabel(y1label, loc='bottom')
    axleftlin.set_ylim((y1lim[1], y1lim[2]))
    axleftlog.set_ylim((y1lim[0], y1lim[1]))

    if not isinstance(y1, list): y1 = [y1]
    handles1 = []
    for i, y in enumerate(y1):
        ls = '-' if i == 0 else '--'
        c = 'C'+str(i) if i < 3 else 'C'+str(i+1)
        axleftlin.plot(x, y, ls=ls, color=c)
        h, = axleftlog.plot(x, y, ls=ls, color=c)
        handles1.append(h)
    legend1 = plt.legend(handles1, y1legends, loc='upper left',
                         bbox_to_anchor=(0.6,3.1), frameon=False)
    axleftlog.add_artist(legend1)
    fig = axes.figure
    if tightlayout:
        fig.tight_layout()

    return axleftlin, axrightlin, axleftlog, axrightlog

@openfigureifnoaxes(projection=None)
def radialprofile_single_scale_axes(x,
                                    y1,
                                    y2=None,
                                    xlim=None,
                                    y1lim=None,
                                    y2lim=None,
                                    xlabel=None,
                                    y1label=None,
                                    y2label=None,
                                    y1legends=None,
                                    y2legends=None,
                                    axes=None,
                                    tightlayout=True, yscale="linear",
                                    ):
    """Plot 1D profiles in the same figure with two axes. Each y-axis has only
    linear or logarithmic scale. For contrast, see: radialprofile()

    Parameters
    ----------
    x : array_like
        Marker x-coordinates.
    y1 : array_like
        Data to be plotted.
    y2 : array_like, optional
        Data to be plotted.
    xlim : array_like, optional
        x-axis limits.
    y1lim : array_like, optional
        y1-axis limits.
    y2lim : array_like, optional
        y2-axis limits.
    xlabel : str, optional
        Label for the x-axis.
    y1label : str, optional
        Label for the y1-axis.
    y2label : str, optional
        Label for the y2-axis.
    y1legends : list of str, optional
        Legends for the y1-axis.
    y2legends : list of str, optional
        Legends for the y2-axis.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    tightlayout : bool, optional
        If True, adjust the subplot parameters to provide certain padding.
    yscale : str, optional
        Either "linear" or "log".
    """
    if xlim is not None: axes.set_xlim(xlim)
    axleft=axes
    axleft.set_yscale(yscale)
    axleft.spines['right'].set_visible(False)
    axleft.spines['bottom'].set_visible(False)
    axleft.set_xlabel(xlabel)

    if y2 is not None:
        # Create linear right axis
        axright = axes.twinx()
        axright.set_yscale(yscale)
        axright.spines['left'].set_visible(False)
        axright.spines['bottom'].set_visible(False)

        if y2lim is not None: axright.set_ylim((y2lim[0], y2lim[1]))
        axright.set_ylabel(y2label)

        axright.spines['right'].set_color('C3')
        axright.tick_params(axis='y', colors='C3')
        axright.yaxis.label.set_color('C3')

        if not isinstance(y2, list): y2 = [y2]
        handles2 = []
        for i, y in enumerate(y2):
            ls = '-' if i == 0 else '--'
            h, = axright.plot(x, y, color='C3', ls=ls)
            handles2.append(h)
        axright.set_xlim(xlim)

    axleft.set_ylabel(y1label)
    if y1lim is not None: axleft.set_ylim((y1lim[0], y1lim[1]))

    if not isinstance(y1, list): y1 = [y1]
    handles1 = []
    for i, y in enumerate(y1):
        ls = '-' if i == 0 else '--'
        c = 'C'+str(i) if i < 3 else 'C'+str(i+1)
        h, = axleft.plot(x, y, ls=ls, color=c)
        handles1.append(h)

    all_handles = handles1 + handles2
    # Here we assume that all or none of the legends should be shown
    if y1legends is not None and y2legends is not None:
        all_legend_lagels = y1legends + y2legends
        axleft.legend(all_handles,
                      all_legend_lagels,
                      loc='upper right',
                      frameon=False)

    fig = axes.figure
    if tightlayout:
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.25)

    return axleft, axright

def defaultcamera(wallmesh):
    """Get default camera (helper function for the 3D plots).

    Default camera is located at R = (Rmax + Rmin) / 2, phi = 0 deg,
    z = (zmax, zmin) / 2, where the min/max values are taken from the bounding
    box of the wall mesh. The focal point is at same (R, z) coordinates but
    phi = 10 deg. The camera angle parameters are all set to zero.

    Parameters
    ----------
    wallmesh : :obj:`~pyvista.Polydata`
        The wall mesh.

    Returns
    -------
    cpos : array_like
        Camera position coordinates [x, y, z].
    cfoc : array_like
        Camera focal point coordinates [x, y, z].
    cang : array_like
        Camera angle [azimuth, elevation, roll].
    """
    # Find min/max values
    Rmax = np.sqrt( np.amax(wallmesh.points[:,0]**2 + wallmesh.points[:,1]**2) )
    Rmin = np.sqrt( np.amin(wallmesh.points[:,0]**2 + wallmesh.points[:,1]**2) )
    zmax = np.amax( wallmesh.points[:,2] )
    zmin = np.amin( wallmesh.points[:,2] )

    # Camera position in (R,phi,z)
    cpos = np.array([( Rmax + Rmin ) / 2,  0, ( zmax + zmin ) / 2])
    cfoc = np.array([( Rmax + Rmin ) / 2, 10, ( zmax + zmin ) / 2])
    cang = np.array([0, 0, -90])

    # Coordinate transformation (R,phi,z) -> (x,y,z)
    cpos = np.array([cpos[0] * np.cos( np.pi * cpos[1] / 180 ),
                     cpos[0] * np.sin( np.pi * cpos[1] / 180 ),
                     cpos[2]])
    cfoc = np.array([cfoc[0] * np.cos( np.pi * cfoc[1] / 180 ),
                     cfoc[0] * np.sin( np.pi * cfoc[1] / 180 ),
                     cfoc[2]])

    return (cpos, cfoc, cang)

def cameracontrols(plotter):
    """Set FPS camera controls vor interactive plotting.

    Parameters
    ----------
    plotter : :obj:`~pyvista.Plotter`
        The active plotter.
    """
    def control_camera(action):
        cpos = np.array(plotter.camera.position)
        cfoc = np.array(plotter.camera.focal_point)
        cup  = np.array(plotter.camera.up)
        dir = np.array(plotter.camera.direction)
        if action == 'move_forward':
            cpos += dir * 0.2
            cfoc += dir * 0.2
        elif action == 'move_backward':
            cpos -= dir * 0.2
            cfoc -= dir * 0.2
        elif action == 'move_left':
            cpos -= np.cross(dir, cup) * 0.2
            cfoc -= np.cross(dir, cup) * 0.2
        elif action == 'move_right':
            cpos += np.cross(dir, cup) * 0.2
            cfoc += np.cross(dir, cup) * 0.2
        elif action == 'move_up':
            cpos += cup * 0.2
            cfoc += cup * 0.2
        elif action == 'move_down':
            cpos -= cup * 0.2
            cfoc -= cup * 0.2
        elif action == 'rotate_cw':
            cup += np.cross(dir, cup) * 0.05
            cup /= np.sqrt(np.sum(cup**2))
        elif action == 'rotate_ccw':
            cup -= np.cross(dir, cup) * 0.05
            cup /= np.sqrt(np.sum(cup**2))
        elif action == 'look_up':
            vec   = np.cross(dir, cup)
            cfoc += cup * 0.05
            cup   = np.cross(vec, cfoc - cpos)
            cup  /= np.sqrt(np.sum(cup**2))
        elif action == 'look_down':
            vec   = np.cross(dir, cup)
            cfoc -= cup * 0.05
            cup   = np.cross(vec, cfoc - cpos)
            cup  /= np.sqrt(np.sum(cup**2))
        elif action == 'look_left':
            vec   = np.cross(dir, cup)
            cfoc -= vec * 0.05
        elif action == 'look_right':
            vec   = np.cross(dir, cup)
            cfoc += vec * 0.05

        plotter.camera.position    = cpos
        plotter.camera.focal_point = cfoc
        plotter.camera.up = cup
        plotter.update()
        plotter.disable() # This disables some default keys

    # Many of the keys are already taken, hence these unorthodox keybindings
    plotter.add_key_event('d', lambda : control_camera('move_forward'))
    plotter.add_key_event('x', lambda : control_camera('move_backward'))
    plotter.add_key_event('z', lambda : control_camera('move_left'))
    plotter.add_key_event('c', lambda : control_camera('move_right'))
    plotter.add_key_event('n', lambda : control_camera('move_up'))
    plotter.add_key_event('m', lambda : control_camera('move_down'))
    plotter.add_key_event('o', lambda : control_camera('rotate_cw'))
    plotter.add_key_event('u', lambda : control_camera('rotate_ccw'))
    plotter.add_key_event('i', lambda : control_camera('look_up'))
    plotter.add_key_event('k', lambda : control_camera('look_down'))
    plotter.add_key_event('j', lambda : control_camera('look_left'))
    plotter.add_key_event('l', lambda : control_camera('look_right'))
