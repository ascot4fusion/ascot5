"""Module for generating common plots with ASCOT5.
"""
import numpy as np
import warnings

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import mpl_toolkits
except ImportError:
    warnings.warn("Could not import matplotlib. Plotting disabled.")

try:
    import pyvista as pv
except ImportError:
    warnings.warn("Could not import pyvista. 3D wall plotting disabled.")

from mpl_toolkits.axes_grid1 import make_axes_locatable
from functools import wraps

def setpaperstyle(height=5, halfpage=False):
    """TODO make this style suitable for publications
    """
    mpl.style.use(
        {"figure.figsize":(8., 6.),
        "figure.autolayout":True,
        "font.family":"serif",
        "font.serif":"ComputerModern",
        "text.usetex":True,
        "pdf.fonttype":42,
        "ps.fonttype":42,
        "axes.labelsize":28,
        "axes.titlesize":28,
        "axes.titlepad":12,
        "xtick.labelsize":24,
        "ytick.labelsize":24,
        "axes.labelpad":6,
        "legend.fontsize":22,
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
        "ytick.major.pad":5.6,}
    )

def setguistyle():
    """Syle used in GUI.
    """
    mpl.style.use({
        "figure.autolayout":False,
        "font.family":"serif",
        "font.serif":"ComputerModern",
        "text.usetex":True,
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
    })

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
    ticks = [];
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

    ticks = [];
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

    h,_,_,m = axes.hist2d(x, y, bins=[xbins, ybins], weights=weights)

    norm = None
    if logscale: norm = mpl.colors.LogNorm(np.amin(h), np.amax(h))
    cbar = plt.colorbar(m, norm=norm, ax=axes, cax=cax)
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
           clabel=None, clim=[None, None], cmap=None, axesequal=False,
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
    z = np.ma.masked_invalid(z)
    if clim[0] is None:
        clim[0] = np.nanmin(z)
    if clim[1] is None:
        clim[1] = np.nanmax(z)

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

def contour2d():
    """Plot contour on 2D plane.
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
             ylabel=None, clabel=None, axesequal=False, axes=None, cax=None):
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
                      linestyle="None", marker=".", markersize=1)

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
                      linestyle="None", marker="o", markersize=1)
            i1 = i2

        # Plot confined
        uids  = np.unique(ids[i2:])
        cpick = np.arange(nc)
        np.random.shuffle(cpick)
        for i in range(nc):
            idx = np.in1d(ids, uids[i::nc])
            axes.plot(x[idx], y[idx], color=colours[nc_b+cpick[i]],
                      linestyle="None", marker=".", markersize=1)

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    if axesequal:
        axes.set_aspect("equal", adjustable="box")

@openfigureifnoaxes(projection=None)
def still(wallmesh, points=None, data=None, log=False, cpos=None, cfoc=None,
          cang=None, axes=None, cax=None):
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
    """
    p = pv.Plotter(off_screen=True)
    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9])
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9,0.9,0.9])
        p.add_mesh(wallmesh, scalars=data, cmap=cmap, log_scale=log)
        p.remove_scalar_bar()

        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])

    if points is not None:
        p.theme.color = 'black'
        p.add_points(points, render_points_as_spheres=True, point_size=10)

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
            norm = mpl.colors.LogNorm(vmin=minval, vmax=maxval)
        else:
            norm = mpl.colors.Normalize(vmin=0, vmax=maxval)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, ax=axes, cax=cax)
        #cbar.set_label(r"Energy load J/m$^2$")


def interactive(wallmesh, *args, points=None, data=None, log=False, cpos=None,
                cfoc=None, cang=None):
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
    cpos : array_like, optional
        Camera position coordinates [x, y, z].
    cfoc : array_like, optional
        Camera focal point coordinates [x, y, z].
    cang : array_like, optional
        Camera angle [azimuth, elevation, roll].
    """
    p = pv.Plotter()

    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9], log_scale=log)
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9,0.9,0.9])
        p.add_mesh(wallmesh, scalars=data, cmap=cmap, log_scale=log)

    if points is not None:
        p.theme.color = 'black'
        p.add_points(points, render_points_as_spheres=True, point_size=10)

    # Set events
    for i in range(len(args)):
        def wrapper(*wargs):
            args[i][1](p)

        p.add_key_event(args[i][0], wrapper)

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
    """
    """
    idx = np.argsort(-loads)
    wetted = np.cumsum(wetted[idx])
    loads  = loads[idx]

    axes.set_xscale('log')
    axes.set_yscale('linear')
    #axes.set_ylim((1e4, 2e8))
    #ax.set_xlim((1,14))
    #ax.spines['left'].set_visible(False)
    #ax.yaxis.set_ticks_position('right')
    #ax.yaxis.set_visible(False)

    #divider = make_axes_locatable(axes)
    #axes1 = divider.append_axes("left", size=1.8, pad=0, sharex=axes)
    #axes1.set_yscale('log')
    #ax1.set_xscale('linear')
    #ax1.set_xlim((1e-4, 9.99e-1))
    #ax1.spines['right'].set_visible(False)
    #ax1.yaxis.set_ticks_position('left')
    #plt.setp(ax1.get_xticklabels(), visible=True)
    axes.plot(loads, wetted)
    axes.set_xlabel(r"Wet area [m$^2$]")
    axes.set_ylabel(r"Wet area [m$^2$]")

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
    cang = np.array([0, 0, 0])

    # Coordinate transformation (R,phi,z) -> (x,y,z)
    cpos = np.array([cpos[0] * np.cos( np.pi * cpos[1] / 180 ),
                     cpos[0] * np.sin( np.pi * cpos[1] / 180 ),
                     cpos[2]])
    cfoc = np.array([cfoc[0] * np.cos( np.pi * cfoc[1] / 180 ),
                     cfoc[0] * np.sin( np.pi * cfoc[1] / 180 ),
                     cfoc[2]])

    return (cpos, cfoc, cang)
