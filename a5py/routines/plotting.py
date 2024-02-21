"""Module for generating common plots with ASCOT5.

This module should be imported everywhere where plotting is done instead of
using matplotlib directly. The reason is that on some platforms matplotlib
is not available and even there we want to able to use all functionality
that doesn't require plotting.

So either import this module or use try-except when importing matplotlib.
"""
import numpy as np
import warnings

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
            if axes is None:
                fig  = plt.figure()
                axes = fig.add_subplot(111, projection=projection)
                plotfun(*args, axes=axes, **kwargs)
                plt.show()
            else:
                if projection != None and axes.name != projection:
                    raise ValueError(
                        "Invalid projection \"%s\" on axes: expected \"%s\"" %
                        (axes.name, projection))
                plotfun(*args, axes=axes, **kwargs)

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

@openfigureifnoaxes(projection=None)
def scatter2d(x, y, c=None, xlog="linear", ylog="linear", clog="linear",
              xlabel=None, ylabel=None, clabel=None, cint=9, cmap=None,
              axesequal=False, axes=None, cax=None):
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
    """
    axes.set_xscale(xlog)
    axes.set_yscale(ylog)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if axesequal: axes.set_aspect("equal", adjustable="box")

    if c is None or isinstance(c, str):
        # Simple plot with a single color
        axes.plot(x, y, color=c, linestyle="None", marker="o")
        return

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
                  marker="o")
        i1 = i2

    # Make colorbar with where color scale has the intervals
    norm = mpl.colors.BoundaryNorm(cint, nc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    ticks = []
    for b in cint.v:
        log = np.floor(np.log10(np.abs(b)))
        mul = b / 10**log
        ticks.append(r"$%.2f\times 10^{%.0f}$" % (mul, log))
    cbar.ax.set_yticklabels(ticks)
    cbar.set_label(clabel)

@openfigureifnoaxes(projection="3d")
def scatter3d(x, y, z, c=None, xlog="linear", ylog="linear", zlog="linear",
              clog="linear", xlabel=None, ylabel=None, zlabel=None, clabel=None,
              cint=None, cmap=None, axesequal=False, axes=None, cax=None):
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
        axes.plot(x, y, z, color=c, linestyle="None", marker="o")
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
                  linestyle="None", marker="o")
        i1 = i2

    # Make colorbar with where color scale has the intervals
    norm = mpl.colors.BoundaryNorm(cint, nc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)

    ticks = []
    for b in cint.v:
        log = np.floor(np.log10(np.abs(b)))
        mul = b / 10**log
        ticks.append(r"$%.2f\times 10^{%.0f}$" % (mul, log))
    cbar.ax.set_yticklabels(ticks)
    cbar.set_label(clabel)

@openfigureifnoaxes(projection=None)
def hist1d(x, xbins=None, weights=None, xlog="linear", logscale=False,
           xlabel=None, legend=None, axes=None):
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
              weights=weights, rwidth=2)
    axes.legend(legend, frameon=False)

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

    h,_,_,m = axes.hist2d(x, y, bins=[xbins, ybins], weights=weights, norm=norm)

    cbar = plt.colorbar(m, ax=axes, cax=cax)
    cbar.set_label(clabel)

@openfigureifnoaxes(projection=None)
def mesh1d(x, y, log=False, xlabel=None, ylabel=None, axes=False):
    """Plot 1D distribution.

    Parameters
    ----------
    x : array_like (nx,)
        Abscissa edges for the x-axis.
    z : array_like (nx-1,)
        Data to be plotted.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
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
    axes.plot(xc, yc)

@openfigureifnoaxes(projection=None)
def mesh2d(x, y, z, log=False, diverging=False, xlabel=None, ylabel=None,
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
    log : bool, optional
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

    if log:
        if diverging:
            if cmap == None: cmap = "bwr"
            norm = mpl.colors.SymLogNorm(linthresh=10, linscale=1.0,
                                 vmin=clim[0], vmax=clim[1], base=10)
        else:
            if cmap == None: cmap = "viridis"
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
def contour2d(x, y, z, contours, xlabel=None, ylabel=None, axesequal=False,
              colors=None, linestyles=None, linewidths=None, axes=None):
    """Plot contour on 2D plane.
    """
    axes.contour(x, y, z.T, contours, colors=colors, linestyles=linestyles,
                 linewidths=linewidths)
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
           cmap=None, axesequal=False, axes=None, cax=None):
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
    """
    axes.set_xscale(xlog)
    axes.set_yscale(ylog)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if axesequal: axes.set_aspect("equal", adjustable="box")
    if c is None or isinstance(c, str):
        # Simple plot with a single color
        for i in range(len(x)):
            axes.plot(x[i], y[i], color=c)
        return

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
    plt.colorbar(smap, ax=axes, cax=cax)

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
    plt.colorbar(smap, ax=axes, cax=cax)

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

        cmin = np.floor(cmin)
        cmax = np.ceil(cmax)
        nc_b = int(cmax - cmin)
        clim = np.linspace(cmin, cmax, nc_b+1)

        # Confined markers are plotted separately.
        connlen[idx+1:] = cmax + 1/nc
        clim = np.append(clim, cmax + np.linspace(0, 1/nc, nc))

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
        cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                          spacing='proportional', extend='min', ticks=ticks,
                          boundaries=[clim[0]-0.5] + clim + [clim[-1]+0.5])
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


def interactive(wallmesh, *args, points=None, orbit=None, data=None, log=False,
                clim=None, cpos=None, cfoc=None, cang=None, **kwargs):
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
    cpos : array_like, optional
        Camera position coordinates [x, y, z].
    cfoc : array_like, optional
        Camera focal point coordinates [x, y, z].
    cang : array_like, optional
        Camera angle [azimuth, elevation, roll].
    **kwargs
        Keyword arguments passed to :obj:`~pyvista.Plotter`.
    """
    p = pv.Plotter(**kwargs)
    p.disable()
    cameracontrols(p)

    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9], log_scale=log)
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9,0.9,0.9])
        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])
        if clim is None: clim = [minval, maxval]
        if clim[0] is None: clim[0] = minval
        if clim[1] is None: clim[1] = maxval
        p.add_mesh(wallmesh, scalars=data, cmap=cmap, clim=clim, log_scale=log)

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

    p.show()

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
                  y2legends=None, axes=None):
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
    """
    if y1lim is None: raise ValueError("y1lim must be provided")
    if xlim is not None: axes.set_xlim(xlim)

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
    axleftlin.yaxis.set_major_formatter(getmathtextsciformatter("%1.0e"))

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
        axrightlin.yaxis.set_major_formatter(getmathtextsciformatter("%1.0e"))

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

    # Not all keys are available for us so we make do
    plotter.add_key_event('w', lambda : control_camera('move_forward'))
    plotter.add_key_event('s', lambda : control_camera('move_backward'))
    plotter.add_key_event('a', lambda : control_camera('move_left'))
    plotter.add_key_event('d', lambda : control_camera('move_right'))
    plotter.add_key_event('n', lambda : control_camera('move_up'))
    plotter.add_key_event('m', lambda : control_camera('move_down'))
    plotter.add_key_event('r', lambda : control_camera('rotate_cw'))
    plotter.add_key_event('y', lambda : control_camera('rotate_ccw'))
    plotter.add_key_event('t', lambda : control_camera('look_up'))
    plotter.add_key_event('g', lambda : control_camera('look_down'))
    plotter.add_key_event('f', lambda : control_camera('look_left'))
    plotter.add_key_event('h', lambda : control_camera('look_right'))
