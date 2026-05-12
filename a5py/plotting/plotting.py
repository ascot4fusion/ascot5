"""Plotting routines customed for ASCOT5."""

from typing import Optional, TYPE_CHECKING

import numpy as np
import warnings

from .misc import (
    open_figure_if_no_axes,
    get_ax10b_formatter,
)

from a5py.utils import OptionalDependency

if TYPE_CHECKING:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pyvista
    import mpl_toolkits
else:
    mpl = OptionalDependency("matplotlib")
    plt = OptionalDependency("matplotlib.pyplot")
    pyvista = OptionalDependency("pyvista")
    mpl_toolkits = OptionalDependency("mpl_toolkits")

try:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    warnings.warn("Could not import matplotlib. Plotting disabled.")

def _set_bbox(bbox, xmin, xmax, xscale, ymin, ymax, yscale, zmin=None, zmax=None, zscale=None):
    dx = (xmax - xmin) / 10
    dy = (ymax - ymin) / 10
    if zmin is not None and zmax is not None:
        dz = (zmax - zmin) / 10
    if bbox is None:
        bbox = [None] * 4 if zmin is None else [None] * 6
    else:
        bbox = list(bbox)

    def find_lower_limit(min, d, scale):
        if scale == "log":
            return min * 2 if min < 0 else min / 2
        if min >= 0 and min - d < 0:
            return 0
        return min - d

    def find_upper_limit(max, d, scale):
        if scale == "log":
            return max / 2 if max < 0 else max * 2
        return max + d


    for i in range(len(bbox)):
        if bbox[i] is None:
            match i:
                case 0:
                    bbox[i] = find_lower_limit(xmin, dx, xscale)
                case 1:
                    bbox[i] = find_lower_limit(ymin, dy, yscale)
                case 2:
                    bbox[i] = (
                        find_lower_limit(zmin, dz, zscale) if zmin is not None
                        else find_upper_limit(xmax, dx, xscale)
                    )
                case 3:
                    bbox[i] = (
                        find_upper_limit(xmax, dx, xscale) if zmin is not None
                        else find_upper_limit(ymax, dy, yscale)
                    )
                case 4:
                    bbox[i] = find_upper_limit(ymax, dy, yscale)
                case 5:
                    bbox[i] = find_upper_limit(zmax, dz, zscale)
    return tuple(bbox)

def _set_clim(clim, cmin, cmax):
    if clim is None:
        return (cmin, cmax)
    clim = list(clim)
    if clim[0] is None:
        clim[0] = cmin
    if clim[1] is None:
        clim[1] = cmax
    return tuple(clim)

def _set_axes(axes, equal, bbox, *params):
    """
    Set axes properties (scale, label, aspect).

    Parameters
    ----------
    axes : :obj:`~matplotlib.axes.Axes`
        The axes where figure is plotted or otherwise new figure is created.
    equal : bool
        Whether the plot is in equal scale.
    bbox : tuple[float, float, float, float] or
           tuple[float, float, float, float, float, float], optional
        Bounding box for the plot.
    *params : tuple[str, str, str]
        Name, scale, and label of the axis to be set.
    """
    divergent = {"x": False, "y": False, "z": False}
    if len(bbox) == 6:
        xlim = (bbox[0], bbox[3])
        ylim = (bbox[1], bbox[4])
        zlim = (bbox[2], bbox[5])
        divergent["x"] = bbox[0] < 0 and bbox[3] > 0
        divergent["y"] = bbox[1] < 0 and bbox[4] > 0
        divergent["z"] = bbox[2] < 0 and bbox[5] > 0
    else:
        xlim = (bbox[0], bbox[2])
        ylim = (bbox[1], bbox[3])
        zlim = None
        divergent["x"] = bbox[0] < 0 and bbox[2] > 0
        divergent["y"] = bbox[1] < 0 and bbox[3] > 0

    for coord, scale, label in params:
        getattr(axes, f"set_{coord}label")(label)
        getattr(axes, f"set_{coord}scale")(scale)
        if scale == "log" and divergent[coord]:
            getattr(axes, f"set_{coord}scale")("symlog")

        if scale == "linear":
            formatter = mpl.ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((0, 0))
            getattr(axes, f"{coord}axis").set_major_formatter(formatter)

    if equal:
        axes.set_aspect("equal", adjustable="box")

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    if zlim is not None:
        axes.set_zlim(zlim)


def _set_color(clim, cscale, cmap):
    diverging = clim[0] < 0 and clim[1] > 0
    if cmap is None:
        cmap = "vanimo" if diverging else "viridis"
    if cscale == "log":
        if diverging:
            norm = mpl.colors.SymLogNorm(
                linthresh=10, linscale=1.0, vmin=clim[0], vmax=clim[1], base=10
            )
        else:
            norm = mpl.colors.LogNorm(vmin=clim[0], vmax=clim[1])
    elif cscale == "linear":
        if diverging:
            norm = mpl.colors.CenteredNorm(halfrange=np.amax(np.abs(clim)))
        else:
            norm = mpl.colors.Normalize(vmin=clim[0], vmax=clim[1])
    return norm, cmap

@open_figure_if_no_axes()
def scatter2d(
    x: np.ndarray,
    y: np.ndarray,
    c: Optional[np.ndarray | str]=None,
    color_intervals: int | np.ndarray=9,
    xscale: str="linear",
    yscale: str="linear",
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    equal: bool=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
    """Make a scatter plot in 2D+1 where color can be one dimension.

    For better performance this function uses matplotlib.pyplot.plot function
    instead of matplotlib.pyplot.scatter. From practical point of view, the
    marker size won't vary and the color is not continuous.

    Parameters
    ----------
    x : array_like
        Marker x-coordinates.
    y : array_like
        Marker y-coordinates.
    c : str or array_like, optional
        Color data or string indicating the color.
    xscale : {"linear", "log"}, optional
        Scaling of the x-axis.
    yscale : {"linear", "log"}, optional
        Scaling of the y-axis.
    cscale : {"linear", "log"}, optional
        Scaling of the color axis.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    color_intervals : int or array_like, optional
        Color intervals used if ``c`` contains data.

        Since we are using plot
        instead of data, the color scale can't be continuous.
    cmap : str, optional
        Name of the colormap where nc colors are picked if c contains data.
    equal : bool, optional
        Make the plot in equal scale.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    bbox = _set_bbox(
        bbox, np.nanmin(x), np.nanmax(x), xscale,
        np.nanmin(y), np.nanmax(y), yscale,
        )
    _set_axes(
        axes, equal, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel)
        )

    if c is None or isinstance(c, str):
        axes.plot(x, y, color=c, linestyle="None", marker="o")
        return

    sort_by_color = np.argsort(c)
    x_s, y_s, c_s = (
        x[sort_by_color], y[sort_by_color], c[sort_by_color]
    )
    clim = _set_clim(clim, c_s[0], c_s[-1])
    if isinstance(color_intervals, int):
        if cscale == "log":
            if clim[0] > 0:
                color_intervals = np.logspace(
                    np.log10(clim[0]), np.log10(clim[1]), color_intervals + 1,
                    )
            else:
                transform = mpl.scale.SymmetricalLogTransform(
                    base=10, linthresh=1., linscale=1.,
                    )
                t = np.linspace(
                    transform.transform(clim[0]),
                    transform.transform(clim[1]),
                    color_intervals + 1,
                    )
                color_intervals = transform.inverted().transform(t)
        else:
            color_intervals = np.linspace(clim[0], clim[1], color_intervals)
    n_colors = len(color_intervals) - 1

    if cmap is None:
        cmap = "vanimo" if (clim[0] < 0 < clim[1]) else "viridis"
    cmap = plt.get_cmap(cmap, n_colors)

    idx = np.searchsorted(c_s, color_intervals)
    for i, (start, stop) in enumerate(zip(idx[:-1], idx[1:])):
        axes.plot(
            x_s[start:stop],
            y_s[start:stop],
            linestyle="None",
            marker="o",
            color=cmap(i),
        )

    norm = mpl.colors.BoundaryNorm(color_intervals, n_colors)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    ticks = []
    for limit in color_intervals:
        if limit == 0:
            ticks.append(r"$0$")
        else:
            log = int(np.floor(np.log10(abs(limit))))
            mul = limit / 10**log
            ticks.append(rf"${mul:.2f}\times 10^{{{log}}}$")
    cbar.ax.set_yticklabels(ticks)
    cbar.set_label(clabel)


@open_figure_if_no_axes(projection="3d")
def scatter3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    c: Optional[np.ndarray]=None,
    color_intervals: int | np.ndarray=9,
    xscale: str="linear",
    yscale: str="linear",
    zscale: str="linear",
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    zlabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float, float, float]]=None,
    equal: Optional[plt.Axes]=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
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
        Flag to set aspect ratio of [x: np.ndarray,y] axes equal.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    bbox = _set_bbox(
        bbox, np.nanmin(x), np.nanmax(x), xscale,
        np.nanmin(y), np.nanmax(y), yscale,
        np.nanmin(z), np.nanmax(z), zscale,
        )
    _set_axes(
        axes, equal, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel),
        ("z", zscale, zlabel),
        )

    if c is None or isinstance(c, str):
        axes.plot(x, y, z, color=c, linestyle="None", marker="o")
        return

    sort_by_color = np.argsort(c)
    x_s, y_s, z_s, c_s = (
        x[sort_by_color], y[sort_by_color], z[sort_by_color],c[sort_by_color]
    )
    clim = _set_clim(clim, c_s[0], c_s[-1])
    if isinstance(color_intervals, int):
        if cscale == "log":
            if clim[0] > 0:
                color_intervals = np.logspace(
                    np.log10(clim[0]), np.log10(clim[1]), color_intervals + 1,
                    )
            else:
                transform = mpl.scale.SymmetricalLogTransform(
                    base=10, linthresh=1., linscale=1.,
                    )
                t = np.linspace(
                    transform.transform(clim[0]),
                    transform.transform(clim[1]),
                    color_intervals + 1,
                    )
                color_intervals = transform.inverted().transform(t)
        else:
            color_intervals = np.linspace(clim[0], clim[1], color_intervals)
    n_colors = len(color_intervals) - 1

    if cmap is None:
        cmap = "vanimo" if (clim[0] < 0 < clim[1]) else "viridis"
    cmap = plt.get_cmap(cmap, n_colors)

    idx = np.searchsorted(c_s, color_intervals)
    for i, (start, stop) in enumerate(zip(idx[:-1], idx[1:])):
        axes.plot(
            x_s[start:stop],
            y_s[start:stop],
            z_s[start:stop],
            linestyle="None",
            marker="o",
            color=cmap(i),
        )

    norm = mpl.colors.BoundaryNorm(color_intervals, n_colors)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    ticks = []
    for limit in color_intervals:
        if limit == 0:
            ticks.append(r"$0$")
        else:
            log = int(np.floor(np.log10(abs(limit))))
            mul = limit / 10**log
            ticks.append(rf"${mul:.2f}\times 10^{{{log}}}$")
    cbar.ax.set_yticklabels(ticks)
    cbar.set_label(clabel)


@open_figure_if_no_axes()
def hist1d(
    x: np.ndarray,
    xbins: int | np.ndarray=10,
    weights: Optional[np.ndarray]=None,
    legend: Optional[list[str]]=None,
    xscale: str="linear",
    yscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    axes: Optional[plt.Axes]=None,
):
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
    if isinstance(xbins, int) and xscale == "log":
        xmin, xmax = x[0][0], x[0][0]
        for i in range(len(x)):
            xmin = np.amin([np.amin(x[i]), xmin])
            xmax = np.amax([np.amax(x[i]), xmax])
        xbins = np.logspace(
            np.sign(xmin) * np.log10(np.abs(xmin)),
            np.sign(xmax) * np.log10(np.abs(xmax)),
            xbins,
            )

    bars, xedge, _ = axes.hist(
        x, xbins, density=False, stacked=True, log=yscale == "log", weights=weights, rwidth=2
    )
    bars = bars[-1, :]
    bbox = _set_bbox(bbox, xedge[0], xedge[-1], xscale, np.min(bars), np.max(bars), yscale)
    _set_axes(
        axes, False, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel),
        )
    if legend is not None:
        axes.legend(legend, frameon=False)


@open_figure_if_no_axes()
def hist2d(
    x: np.ndarray,
    y: np.ndarray,
    xbins: np.ndarray | int=10,
    ybins: np.ndarray | int=10,
    weights: Optional[np.ndarray]=None,
    xscale: str="linear",
    yscale: str="linear",
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    equal: bool=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
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
    if isinstance(xbins, int) and xscale == "log":
        xmin, xmax = np.amin(x), np.amax(x)
        xbins = np.logspace(
            np.sign(xmin) * np.log10(np.abs(xmin)),
            np.sign(xmax) * np.log10(np.abs(xmax)),
            xbins,
            )
    if isinstance(ybins, int) and yscale == "log":
        ymin, ymax = np.amin(y), np.amax(y)
        ybins = np.logspace(
            np.sign(ymin) * np.log10(np.abs(ymin)),
            np.sign(ymax) * np.log10(np.abs(ymax)),
            ybins,
            )

    c, xedge, yedge = np.histogram2d(x, y, bins=[xbins, ybins], weights=weights)
    if cscale == "log" and not (np.amin(c) < 0 and np.amax(c) > 0):
        c[c == 0] = np.nan
    mesh2d(
        xedge, yedge, c, xscale=xscale, yscale=yscale, cscale=cscale,
        xlabel=xlabel, ylabel=ylabel, clabel=clabel, bbox=bbox, equal=equal,
        clim=clim, cmap=cmap, axes=axes, cax=cax,
    )


@open_figure_if_no_axes()
def mesh1d(
    x: np.ndarray,
    y: np.ndarray,
    xscale: str="linear",
    yscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    axes: Optional[plt.Axes]=None,
    ):
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
    if (np.diff(x) <= 0).any():
        raise ValueError("x must be monotonically increasing.")
    if x.size - 1 != y.size:
        raise ValueError("Size of x must be one more than size of y.")
    if bbox is None:
        bbox = [None] * 4
    bbox = list(bbox)
    if bbox[0] is None:
        bbox[0] = x[0]
    if bbox[2] is None:
        bbox[2] = x[-1]
    bbox = _set_bbox(bbox, x[0], x[-1], xscale, np.nanmin(y), np.nanmax(y), yscale)
    _set_axes(
        axes, False, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel),
        )
    xc = np.r_[x[0], np.repeat(x[1:-1], 2), x[-1]]
    yc = np.repeat(y, 2)

    axes.plot(xc, yc)


@open_figure_if_no_axes()
def mesh2d(
    x: np.ndarray,
    y: np.ndarray,
    c: np.ndarray,
    xscale: str="linear",
    yscale: str="linear",
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    equal: bool=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
    """Make a mesh (surface) plot in 2D.

    Parameters
    ----------
    x : array_like (nx,) or (nx+1,)
        Abscissa or abscissa edges for the x-axis.
    y : array_like (ny,) or (ny+1,)
        Abscissa or abscissa edges for the y-axis.
    z : array_like (nx,ny)
        Data to be plotted.
    xscale : {"linear", "log"}
        Make x-axis logarithmic.
    yscale : {"linear", "log"}
        Make y-axis logarithmic.
    cscale : {"linear", "log"}
        Make color logarithmic.
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
    if c.shape == (x.size, y.size):
        shading = "nearest"
    elif c.shape == (x.size - 1, y.size - 1):
        shading = "flat"
    else:
        raise ValueError(
            f"Shape of 'c' was {c.shape}. Expected {(x.size, y.size)} or "
            f"{(x.size - 1, y.size - 1)}.",
            )
    c = np.ma.masked_invalid(c)
    clim = _set_clim(clim, np.nanmin(c), np.nanmax(c))
    norm, cmap = _set_color(clim, cscale, cmap)

    pcm = axes.pcolormesh(x, y, c.T, norm=norm, cmap=cmap, shading=shading)
    axes.patch.set(hatch="x", edgecolor=[0.9, 0.9, 0.9])

    if bbox is None:
        bbox = [None] * 4
    bbox = list(bbox)
    if bbox[0] is None:
        bbox[0] = x[0] if shading == "flat" else x[0] - 0.5 * np.diff(x)[0]
    if bbox[1] is None:
        bbox[1] = y[0] if shading == "flat" else y[0] - 0.5 * np.diff(y)[0]
    if bbox[2] is None:
            bbox[2] = x[-1] if shading == "flat" else x[-1] + 0.5 * np.diff(x)[-1]
    if bbox[3] is None:
            bbox[3] = y[-1] if shading == "flat" else y[-1] + 0.5 * np.diff(y)[-1]
    _set_axes(
        axes, equal, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel),
        )

    cbar = plt.colorbar(pcm, ax=axes, cax=cax)
    cbar.set_label(clabel)


@open_figure_if_no_axes()
def contour2d(
    x: np.ndarray,
    y: np.ndarray,
    c: np.ndarray,
    contours: np.ndarray,
    fill: bool=False,
    extend_min: Optional[float]=None,
    extend_max: Optional[float]=None,
    colors: str | list[str]="black",
    linestyles: str | list[str]="solid",
    linewidths: float | np.ndarray=1.,
    labels: bool | list[str]=False,
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    xscale: str="linear",
    yscale: str="linear",
    cscale: str="linear",
    bbox: Optional[tuple[float, float, float, float]]=None,
    equal: bool=False,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
    """
    Plot contour on 2D plane.

    Parameters
    ----------
    x : array_like (nx,)
        Abscissa for the x-axis.
    y : array_like (ny,)
        Abscissa for the y-axis.
    c : array_like (nx,ny)
        Data to be plotted.
    contours : array_like
        Contour levels.
    fill : bool, optional
        If true, fill the contour with color.
    extend_min : float, optional
        Color for the minimum value.
    extend_max : float, optional
        Color for the maximum value.
    colors : str | list[str], optional
        Colors for the contour lines.
    linestyles : str | list[str], optional
        Line styles for the contour lines.
    linewidths : float | array_like, optional
        Line widths for the contour lines.
    labels : bool | list[str], optional
        Labels for the contour lines.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color-axis.
    xscale : {"linear", "log"}, optional
        x-axis scaling.
    yscale : {"linear", "log"}, optional
        y-axis scaling.
    cscale : {"linear", "log"}, optional
        Color scaling.
    bbox : tuple[float, float, float, float], optional
        Bounding box for the plot.
    equal : bool, optional
        Flag to set the aspect ratio equal.
    cmap : str, optional
        Colormap.
    axes : plt.Axes, optional
        Axes for the plot.
    cax : plt.Axes, optional
        Axes for the color bar.
    """
    if isinstance(contours, int):
        if cscale == "linear":
            contours = np.linspace(c.min(), c.max(), contours)
        else:
            if np.amin(c) > 0:
                contours = np.logspace(
                    np.log10(np.amin(c)), np.log10(np.amax(c)), contours + 1,
                    )
            else:
                transform = mpl.scale.SymmetricalLogTransform(
                    base=10, linthresh=1., linscale=1.,
                    )
                t = np.linspace(
                    transform.transform(np.amin(c)),
                    transform.transform(np.amax(c)),
                    contours,
                    )
                contours = transform.inverted().transform(t)
    if fill:
        norm, cmap = _set_color((contours[0], contours[-1]), cscale, cmap)
        cmap = plt.get_cmap(cmap)
        extend = "neither"
        if extend_max is not None:
            cmap.set_over(extend_max)
            extend = "max"
        if extend_min is not None:
            cmap.set_under(extend_min)
            extend = "min" if extend == "neither" else "both"
        cf = axes.contourf(
            x, y, c, contours, norm=norm, cmap=cmap, extend=extend,
        )
        cbar = plt.colorbar(cf, cax=cax)
        if cscale == "log":
            ticks = []
            for limit in contours:
                if limit == 0:
                    ticks.append(r"$0$")
                else:
                    log = int(np.floor(np.log10(abs(limit))))
                    mul = limit / 10**log
                    ticks.append(rf"${mul:.2f}\times 10^{{{log}}}$")
            cbar.ax.set_yticklabels(ticks)
        cbar.set_label(clabel)
    cs = axes.contour(
        x, y, c, contours, colors=colors, linestyles=linestyles,
        linewidths=linewidths,
    )
    if labels:
        fmt = None
        if isinstance(labels, list):
            fmt = {}
            for l, s in zip(cs.levels, labels):
                fmt[l] = s
        axes.clabel(cs, cs.levels, fmt=fmt)
    if bbox is None:
        bbox = [None] * 4
    if bbox[0] is None:
        bbox[0] = np.amin(x)
    if bbox[1] is None:
        bbox[1] = np.amin(y)
    if bbox[2] is None:
        bbox[2] = np.amax(x)
    if bbox[3] is None:
        bbox[3] = np.amax(y)
    _set_axes(
        axes, equal, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel),
        )


@open_figure_if_no_axes()
def arrows2d(
    x: np.ndarray,
    y: np.ndarray,
    dx: np.ndarray,
    dy: np.ndarray,
    c: Optional[str | np.ndarray]=None,
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    equal: bool=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
    ):
    """
    Plot 2D arrows.

    Parameters
    ----------
    start : (2,) or (N, 2) array-like
        Arrow starting point(s).
    vec_or_end : (N, 2) array-like
        Either direction vectors or arrow endpoints.
    colors : None, color spec, or sequence of color specs
        Optional arrow color(s).
    endpoints : bool
        If True, vec_or_end is interpreted as endpoints.
    cscale : {"linear", "log"}
        Make color logarithmic.
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
    x = np.asarray(x)
    y = np.asarray(y)
    dx = np.asarray(dx)
    dy = np.asarray(dy)
    n = x.size
    if dx.size != n or y.size != n or dy.size != n:
        raise ValueError(
            "Input arrays 'x', 'y', 'dx', and 'dy' must have the same size"
            f" (was {x.size}, {y.size}, {dx.size}, and {dy.size})."
            )

    xmin = np.amin((np.amin(x), np.amin(x + dx)))
    xmax = np.amax((np.amax(x), np.amax(x + dx)))
    ymin = np.amin((np.amin(y), np.amin(y + dy)))
    ymax = np.amax((np.amax(y), np.amax(y + dy)))
    bbox = _set_bbox(
        bbox, xmin, xmax, "linear", ymin, ymax, "linear",
        )
    _set_axes(
        axes, equal, bbox, ("x", "linear", xlabel), ("y", "linear", ylabel)
        )

    params = [x, y, dx, dy]
    colors, norm = None, None
    if c is not None:
        if isinstance(c, str):
            colors = c
        else:
            if len(c) != n:
                raise ValueError(
                    f"'c' has {len(c)} elements (expected {n})."
                    )
            params.append(c)
            clim = _set_clim(clim, np.nanmin(c), np.nanmax(c))
            norm, cmap = _set_color(clim, cscale, cmap)

    q = axes.quiver(
        *params,
        angles="xy",
        scale_units="xy",
        scale=1,
        color=colors,
        norm=norm,
        cmap=cmap,
        )

    if c is not None and colors is None:
        cbar = plt.colorbar(q, ax=axes, cax=cax)
        cbar.set_label(clabel)



@open_figure_if_no_axes()
def line2d(
    x: list[np.ndarray],
    y: list[np.ndarray],
    c: Optional[list[np.ndarray] | str]=None,
    xscale: str="linear",
    yscale: str="linear",
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float]]=None,
    equal: bool=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
    """Plot line segments in 2D.

    Parameters
    ----------
    x : list[array_like]
        x-coordinates of the line segments for each line.
    y : list[array_like]
        y-coordinates of the line segments for each line.
    c : list[array_like] or str, optional
        Color data for each line or string indicating the color.
    xscale : {"linear", "log", "symlog"}, optional
        x-axis scaling.
    yscale : {"linear", "log", "symlog"}, optional
        y-axis scaling.
    cscale : {"linear", "log", "symlog"}, optional
        color-axis scaling.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    bbox : tuple[float, float, float, float], optional
        Bounding box for the plot limits.
    equal : bool, optional
        Flag to set aspect ratio of [x,y] axes equal.
    clim : tuple[float, float], optional
        Color limits.
    cmap : str, optional
        Name of the colormap.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same number of lines.")

    nlines = len(x)
    ismulticolored = not (c is None or isinstance(c, str))
    if ismulticolored and nlines != len(c):
        raise ValueError("The number of colors must match the number of lines.")

    if not ismulticolored:
        xmin, xmax = np.amin(x[0]), np.amax(x[0])
        ymin, ymax = np.amin(y[0]), np.amax(y[0])
        for i in range(nlines):
            xmin = min(xmin, np.amin(x[i]))
            xmax = max(xmax, np.amax(x[i]))
            ymin = min(ymin, np.amin(y[i]))
            ymax = max(ymax, np.amax(y[i]))
            axes.plot(x[i], y[i], color=c)
    else:
        xmin, xmax = np.amin(x[0]), np.amax(x[0])
        ymin, ymax = np.amin(y[0]), np.amax(y[0])
        cmin, cmax = np.amin(c[0]), np.amax(c[0])
        for i in range(nlines):
            xmin = min(xmin, np.amin(x[i]))
            xmax = max(xmax, np.amax(x[i]))
            ymin = min(ymin, np.amin(y[i]))
            ymax = max(ymax, np.amax(y[i]))
            cmin = min(cmin, np.amin(c[i]))
            cmax = max(cmax, np.amax(c[i]))

    bbox = _set_bbox(bbox, xmin, xmax, xscale, ymin, ymax, yscale)
    _set_axes(
        axes, equal, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel)
        )
    if not ismulticolored:
        return

    clim = _set_clim(clim, cmin, cmax)
    norm, cmap = _set_color(clim, cscale, cmap)
    for i in range(nlines):
        points = np.array([x[i], y[i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = mpl.collections.LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(c[i][1:])
        axes.add_collection(lc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    cbar.set_label(clabel)


@open_figure_if_no_axes(projection="3d")
def line3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    c: Optional[np.ndarray]=None,
    xscale: str="linear",
    yscale: str="linear",
    zscale: str="linear",
    cscale: str="linear",
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    zlabel: Optional[str]=None,
    clabel: Optional[str]=None,
    bbox: Optional[tuple[float, float, float, float, float, float]]=None,
    equal: bool=False,
    clim: Optional[tuple[float, float]]=None,
    cmap: Optional[str]=None,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
    """Plot line segments in 3D.

    Parameters
    ----------
    x : list[array_like]
        x-coordinates of the line segments for each line.
    y : list[array_like]
        y-coordinates of the line segments for each line.
    z : list[array_like]
        z-coordinates of the line segments for each line.
    c : list[array_like] or str, optional
        Color data for each line or string indicating the color.
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
    if len(x) != len(y) or len(x) != len(z):
        raise ValueError("x, y, and z must have the same number of lines.")

    nlines = len(x)
    ismulticolored = not (c is None or isinstance(c, str))
    if ismulticolored and nlines != len(c):
        raise ValueError("The number of colors must match the number of lines.")

    if not ismulticolored:
        xmin, xmax = np.amin(x[0]), np.amax(x[0])
        ymin, ymax = np.amin(y[0]), np.amax(y[0])
        zmin, zmax = np.amin(z[0]), np.amax(z[0])
        for i in range(nlines):
            xmin = min(xmin, np.amin(x[i]))
            xmax = max(xmax, np.amax(x[i]))
            ymin = min(ymin, np.amin(y[i]))
            ymax = max(ymax, np.amax(y[i]))
            zmin = min(zmin, np.amin(z[i]))
            zmax = max(zmax, np.amax(z[i]))
            axes.plot(x[i], y[i], z[i], color=c)
    else:
        xmin, xmax = np.amin(x[0]), np.amax(x[0])
        ymin, ymax = np.amin(y[0]), np.amax(y[0])
        zmin, zmax = np.amin(z[0]), np.amax(z[0])
        cmin, cmax = np.amin(c[0]), np.amax(c[0])
        for i in range(nlines):
            xmin = min(xmin, np.amin(x[i]))
            xmax = max(xmax, np.amax(x[i]))
            ymin = min(ymin, np.amin(y[i]))
            ymax = max(ymax, np.amax(y[i]))
            zmin = min(zmin, np.amin(z[i]))
            zmax = max(zmax, np.amax(z[i]))
            cmin = min(cmin, np.amin(c[i]))
            cmax = max(cmax, np.amax(c[i]))

    bbox = _set_bbox(bbox, xmin, xmax, xscale, ymin, ymax, yscale, zmin, zmax, zscale)
    _set_axes(
        axes, equal, bbox, ("x", xscale, xlabel), ("y", yscale, ylabel), ("z", zscale, zlabel)
        )
    if not ismulticolored:
        return

    clim = _set_clim(clim, cmin, cmax)
    norm, cmap = _set_color(clim, cscale, cmap)
    for i in range(len(x)):
        points = np.array([x[i], y[i], z[i]]).T.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = mpl_toolkits.mplot3d.art3d.Line3DCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(c[i][1:])
        axes.add_collection(lc)
    smap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(smap, ax=axes, cax=cax)
    cbar.set_label(clabel)


@open_figure_if_no_axes()
def poincare(
    x: np.ndarray,
    y: np.ndarray,
    ids: np.ndarray,
    connlen: Optional[np.ndarray]=None,
    xlim=None,
    ylim=None,
    xlabel: Optional[str]=None,
    ylabel: Optional[str]=None,
    clabel: Optional[str]=None,
    equal: bool=False,
    markersize: int=2,
    axes: Optional[plt.Axes]=None,
    cax: Optional[plt.Axes]=None,
):
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
    nc = 5  # How many colors are used

    if connlen is not None:
        # If we don't have lost markers then ignore the connection length.
        if np.argwhere(connlen > 0).size == 0:
            connlen = None

    if connlen is None:
        # Shuffle colors and plot
        cmap = plt.cm.get_cmap("Reds", nc)
        uids = np.unique(ids)
        cpick = np.arange(nc)
        np.random.shuffle(cpick)
        for i in range(nc):
            idx = np.in1d(ids, uids[i::nc])
            axes.plot(
                x[idx],
                y[idx],
                color=cmap(cpick[i] / nc),
                linestyle="None",
                marker=".",
                markersize=markersize,
            )

    else:
        # Sort by connection length (confined markers are indicated with a
        # negative connlen).
        idx = np.argsort(connlen)
        x = x[idx]
        y = y[idx]
        ids = ids[idx]
        connlen = connlen[idx]

        # Find where the line between confined and lost markers is.
        # Set connlen positive again for confined markers and rearrange again
        # by the connection length.
        idx = np.argwhere(connlen > 0).ravel()[0]
        x = np.append(x[idx:], np.flip(x[:idx]))
        y = np.append(y[idx:], np.flip(y[:idx]))
        ids = np.append(ids[idx:], np.flip(ids[:idx]))
        connlen = np.append(connlen[idx:], np.flip(-connlen[:idx]))
        idx = (connlen.size - 1) - idx

        # The color has meaning only for lost markers so find the scale
        cmin = connlen[0]
        cmax = connlen[idx]

        logscale = False
        if cmin / cmax < 0.1:
            logscale = True
            cmin = np.log10(cmin)
            cmax = np.log10(cmax)
            connlen = np.log10(connlen)
            cunits = 1
        else:
            cunits = cmax.units

        cmin = np.floor(cmin)
        cmax = np.ceil(cmax)
        nc_b = int(cmax - cmin)
        clim = np.linspace(cmin, cmax, nc_b + 1)

        # Confined markers are plotted separately.
        connlen[idx + 1 :] = cmax + (1 / nc) * cunits
        clim = np.append(clim, cmax + np.linspace(0, 1 / nc, nc) * cunits)

        # Create colormap and colorbar
        cmapred = plt.cm.get_cmap("Reds_r")
        cmapblue = plt.cm.get_cmap("Blues")
        colours = [None] * (nc_b + nc)
        for i in range(nc_b):
            colours[i] = cmapblue(i / nc_b)
        for i in range(nc):
            colours[i + nc_b] = cmapred(i / nc)

        cmap = mpl.colors.ListedColormap(colours)
        norm = mpl.colors.BoundaryNorm(boundaries=clim, ncolors=nc_b + nc)
        ticks = clim[: nc_b + 1]
        if logscale:
            ticklabels = list(clim[:nc_b]) + [r"$\inf$"]
        else:
            ticklabels = list(clim[:nc_b]) + [r"$\inf$"]
        deltac = clim[1] - clim[2]
        cb = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            ax=axes,
            spacing="proportional",
            extend="min",
            ticks=ticks,
            boundaries=[clim[0] - deltac / 2] + clim + [clim[-1] + deltac / 2],
        )
        cb.ax.set_yticklabels(ticklabels)
        if clabel is not None:
            clabel = "log10( " + clabel + " )" if logscale else clabel
            cb.ax.set_ylabel(clabel)

        # Find the indices that divide the color range in even intervals
        idx = np.searchsorted(connlen, clim) + 1

        # Now plot markers on to each interval with a different color
        i1 = 0
        for i in range(nc_b):
            i2 = idx[i + 1]
            axes.plot(
                x[i1:i2],
                y[i1:i2],
                color=colours[i],
                linestyle="None",
                marker="o",
                markersize=markersize,
            )
            i1 = i2

        # Plot confined
        uids = np.unique(ids[i2:])
        cpick = np.arange(nc)
        np.random.shuffle(cpick)
        for i in range(nc):
            idx = np.in1d(ids, uids[i::nc])
            axes.plot(
                x[idx],
                y[idx],
                color=colours[nc_b + cpick[i]],
                linestyle="None",
                marker=".",
                markersize=markersize,
            )

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    if axesequal:
        axes.set_aspect("equal", adjustable="box")


@open_figure_if_no_axes()
def still(
    wallmesh,
    points=None,
    orbit=None,
    data=None,
    log=False,
    clim=None,
    cpos=None,
    cfoc=None,
    cang=None,
    axes=None,
    cax=None,
    **kwargs,
):
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
        p.add_mesh(wallmesh, color=[0.9, 0.9, 0.9])
        clim = None
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9, 0.9, 0.9])
        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])
        if clim is None:
            clim = [minval, maxval]
        if clim[0] is None:
            clim[0] = minval
        if clim[1] is None:
            clim[1] = maxval

        p.add_mesh(wallmesh, scalars=data, cmap=cmap, clim=clim, log_scale=log)
        p.remove_scalar_bar()

    if points is not None:
        p.theme.color = "black"
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
        p.camera.azimuth = cang[0]
        p.camera.elevation = cang[1]
        p.camera.roll = cang[2]

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


def interactive(
    wallmesh,
    *args,
    points=None,
    orbit=None,
    data=None,
    log=False,
    clim=None,
    cpos=None,
    cfoc=None,
    cang=None,
    **kwargs,
):
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
        p.add_mesh(wallmesh, color=[0.9, 0.9, 0.9], log_scale=log)
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9, 0.9, 0.9])
        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])
        if clim is None:
            clim = [minval, maxval]
        if clim[0] is None:
            clim[0] = minval
        if clim[1] is None:
            clim[1] = maxval
        p.add_mesh(wallmesh, scalars=data, cmap=cmap, clim=clim, log_scale=log)

    if points is not None:
        p.theme.color = "black"
        p.add_points(points, render_points_as_spheres=True, point_size=10)

    if orbit is not None:
        orbit = pv.lines_from_points(orbit)
        p.add_mesh(orbit, color="red")

    # Set events
    for i in range(len(args)):
        p.add_key_event(args[i][0], lambda: args[i][1](p))

    # Set camera
    if cpos is not None:
        p.camera.position = cpos
    if cfoc is not None:
        p.camera.focal_point = cfoc
    if cang is not None:
        p.camera.azimuth = cang[0]
        p.camera.elevation = cang[1]
        p.camera.roll = cang[2]

    p.show()


@open_figure_if_no_axes()
def loadvsarea(wetted, loads, axes=None):
    """Plot histogram showing minimum load vs area."""
    idx = np.argsort(-loads)
    wetted = np.cumsum(wetted[idx])
    loads = loads[idx]

    axes.set_xscale("log")
    axes.set_yscale("linear")
    axes.plot(loads, wetted)
    axes.set_xlabel(r"Load above [" + str(loads.units) + "]")
    axes.set_ylabel(r"Wetted area [" + str(wetted.units) + "]")


@open_figure_if_no_axes()
def triangularpatch(
    patches,
    color,
    log=False,
    xlim=None,
    ylim=None,
    clim=None,
    xlabel: str=None,
    ylabel: str=None,
    clabel: str=None,
    cmap=None,
    axes=None,
    cax=None,
):
    """Plot triangular patches."""
    if clim is None:
        clim = [None, None]
    if clim[0] is None:
        clim[0] = np.nanmin(color)
    if clim[1] is None:
        clim[1] = np.nanmax(color)

    if log:
        norm = mpl.colors.LogNorm(clim[0], clim[1])
    else:
        norm = mpl.colors.Normalize(clim[0], clim[1])

    coll = mpl.collections.PolyCollection(patches, array=color, norm=norm, cmap=cmap)
    axes.add_collection(coll)

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    plt.colorbar(coll, norm=norm, ax=axes, cax=cax, label=clabel)


@open_figure_if_no_axes(projection="polar")
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
    axes.set_xticks(np.array([0, 45, 90, 135, 180]) * np.pi / 180)
    axes.set_xticklabels([1.0, 0.5, 0.0, -0.5, -1.0])


@open_figure_if_no_axes()
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
        axes.plot(v * np.cos(ang), v * np.sin(ang), color="black")

    p = pnorm_edges[-1]
    for v in pitch_edges:
        axes.plot([0, v * p], [0, np.sqrt(1.0 - v**2) * p], color="black")


@open_figure_if_no_axes()
def radialprofile(
    x: np.ndarray,
    y1,
    y2=None,
    xlim=None,
    y1lim=None,
    y2lim=None,
    xlabel: str=None,
    y1label: str=None,
    y2label: str=None,
    y1legends=None,
    y2legends=None,
    axes=None,
):
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
    if y1lim is None:
        raise ValueError("y1lim must be provided")
    if xlim is not None:
        axes.set_xlim(xlim)

    # Create linear left axis
    axleftlin = axes
    axleftlin.set_yscale("linear")
    axleftlin.spines["right"].set_visible(False)
    axleftlin.spines["bottom"].set_visible(False)
    # Create log left axis
    divider = make_axes_locatable(axes)
    axleftlog = divider.append_axes("bottom", size=1.0, pad=0, sharex=axes)
    axleftlog.set_yscale("log")
    axleftlog.spines["right"].set_visible(False)
    axleftlog.spines["top"].set_visible(False)
    axleftlog.yaxis.set_ticks_position("left")

    axleftlog.tick_params(axis="y", which="minor", left=False)
    axleftlin.yaxis.set_major_formatter(getmathtextsciformatter("%1.0e"))

    axleftlog.set_xlabel(xlabel)
    plt.setp(axleftlin.get_xticklabels(), visible=False)

    if y2 is not None:
        # Create linear right axis
        axrightlin = axes.twinx()
        axrightlin.set_yscale("linear")
        axrightlin.spines["left"].set_visible(False)
        axrightlin.spines["bottom"].set_visible(False)

        # Create log right axis
        divider = make_axes_locatable(axrightlin)
        axrightlog = divider.append_axes("bottom", size=1.0, pad=0, sharex=axes)
        axrightlog.set_yscale("log")
        axrightlog.spines["left"].set_visible(False)
        axrightlog.spines["top"].set_visible(False)
        axrightlog.yaxis.set_ticks_position("right")
        axrightlog.xaxis.set_visible(False)
        axrightlog.set_facecolor("none")
        axrightlog.yaxis.set_label_position("right")

        axrightlog.tick_params(axis="y", which="minor", right=False)
        axrightlin.yaxis.set_major_formatter(getmathtextsciformatter("%1.1e"))

        axlin = axrightlin
        axrightlin = axleftlin
        axleftlin = axlin
        axlog = axrightlog
        axrightlog = axleftlog
        axleftlog = axlog

        axrightlin.set_ylim((y2lim[1], y2lim[2]))
        axrightlog.set_ylim((y2lim[0], y2lim[1]))
        axrightlin.set_ylabel(y2label, loc="bottom")

        axrightlin.spines["right"].set_color("C3")
        axrightlog.spines["right"].set_color("C3")
        axrightlin.tick_params(axis="y", colors="C3")
        axrightlog.tick_params(axis="y", colors="C3")
        axrightlin.yaxis.label.set_color("C3")

        if not isinstance(y2, list):
            y2 = [y2]
        handles2 = []
        for i, y in enumerate(y2):
            ls = "-" if i == 0 else "--"
            (h,) = axrightlin.plot(x, y, color="C3", ls=ls)
            axrightlog.plot(x, y, color="C3", ls=ls)
            handles2.append(h)
        axrightlog.set_xlim(xlim)

        legend2 = plt.legend(
            handles2,
            y2legends,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            frameon=False,
        )
        axleftlog.add_artist(legend2)

    axleftlin.set_ylabel(y1label, loc="bottom")
    axleftlin.set_ylim((y1lim[1], y1lim[2]))
    axleftlog.set_ylim((y1lim[0], y1lim[1]))

    if not isinstance(y1, list):
        y1 = [y1]
    handles1 = []
    for i, y in enumerate(y1):
        ls = "-" if i == 0 else "--"
        c = "C" + str(i) if i < 3 else "C" + str(i + 1)
        axleftlin.plot(x, y, ls=ls, color=c)
        (h,) = axleftlog.plot(x, y, ls=ls, color=c)
        handles1.append(h)
    legend1 = plt.legend(
        handles1, y1legends, loc="upper left", bbox_to_anchor=(0.6, 3.1), frameon=False
    )
    axleftlog.add_artist(legend1)

    return axleftlin, axrightlin, axleftlog, axrightlog


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
    Rmax = np.sqrt(np.amax(wallmesh.points[:, 0] ** 2 + wallmesh.points[:, 1] ** 2))
    Rmin = np.sqrt(np.amin(wallmesh.points[:, 0] ** 2 + wallmesh.points[:, 1] ** 2))
    zmax = np.amax(wallmesh.points[:, 2])
    zmin = np.amin(wallmesh.points[:, 2])

    # Camera position in (R,phi,z)
    cpos = np.array([(Rmax + Rmin) / 2, 0, (zmax + zmin) / 2])
    cfoc = np.array([(Rmax + Rmin) / 2, 10, (zmax + zmin) / 2])
    cang = np.array([0, 0, -90])

    # Coordinate transformation (R,phi,z) -> (x,y,z)
    cpos = np.array(
        [
            cpos[0] * np.cos(np.pi * cpos[1] / 180),
            cpos[0] * np.sin(np.pi * cpos[1] / 180),
            cpos[2],
        ]
    )
    cfoc = np.array(
        [
            cfoc[0] * np.cos(np.pi * cfoc[1] / 180),
            cfoc[0] * np.sin(np.pi * cfoc[1] / 180),
            cfoc[2],
        ]
    )

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
        cup = np.array(plotter.camera.up)
        dir = np.array(plotter.camera.direction)
        if action == "move_forward":
            cpos += dir * 0.2
            cfoc += dir * 0.2
        elif action == "move_backward":
            cpos -= dir * 0.2
            cfoc -= dir * 0.2
        elif action == "move_left":
            cpos -= np.cross(dir, cup) * 0.2
            cfoc -= np.cross(dir, cup) * 0.2
        elif action == "move_right":
            cpos += np.cross(dir, cup) * 0.2
            cfoc += np.cross(dir, cup) * 0.2
        elif action == "move_up":
            cpos += cup * 0.2
            cfoc += cup * 0.2
        elif action == "move_down":
            cpos -= cup * 0.2
            cfoc -= cup * 0.2
        elif action == "rotate_cw":
            cup += np.cross(dir, cup) * 0.05
            cup /= np.sqrt(np.sum(cup**2))
        elif action == "rotate_ccw":
            cup -= np.cross(dir, cup) * 0.05
            cup /= np.sqrt(np.sum(cup**2))
        elif action == "look_up":
            vec = np.cross(dir, cup)
            cfoc += cup * 0.05
            cup = np.cross(vec, cfoc - cpos)
            cup /= np.sqrt(np.sum(cup**2))
        elif action == "look_down":
            vec = np.cross(dir, cup)
            cfoc -= cup * 0.05
            cup = np.cross(vec, cfoc - cpos)
            cup /= np.sqrt(np.sum(cup**2))
        elif action == "look_left":
            vec = np.cross(dir, cup)
            cfoc -= vec * 0.05
        elif action == "look_right":
            vec = np.cross(dir, cup)
            cfoc += vec * 0.05

        plotter.camera.position = cpos
        plotter.camera.focal_point = cfoc
        plotter.camera.up = cup
        plotter.update()
        plotter.disable()  # This disables some default keys

    # Not all keys are available for us so we make do
    plotter.add_key_event("w", lambda: control_camera("move_forward"))
    plotter.add_key_event("s", lambda: control_camera("move_backward"))
    plotter.add_key_event("a", lambda: control_camera("move_left"))
    plotter.add_key_event("d", lambda: control_camera("move_right"))
    plotter.add_key_event("n", lambda: control_camera("move_up"))
    plotter.add_key_event("m", lambda: control_camera("move_down"))
    plotter.add_key_event("r", lambda: control_camera("rotate_cw"))
    plotter.add_key_event("y", lambda: control_camera("rotate_ccw"))
    plotter.add_key_event("t", lambda: control_camera("look_up"))
    plotter.add_key_event("g", lambda: control_camera("look_down"))
    plotter.add_key_event("f", lambda: control_camera("look_left"))
    plotter.add_key_event("h", lambda: control_camera("look_right"))
