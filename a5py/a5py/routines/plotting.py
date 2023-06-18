"""Module for generating common plots with ASCOT5.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyvista as pv

from mpl_toolkits.axes_grid1 import make_axes_locatable
from functools import wraps

def openfigureifnoaxes(projection=None):
    """Decorator for creating and displaying a new figure if axes are not
    provided.

    Parameters
    ----------
    projection : str, {None, 'aitoff', 'hammer', 'lambert', 'mollweide',
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
            if not axes:
                fig  = plt.figure()
                axes = fig.add_subplot(111, projection=projection)
                plotfun(*args, axes=axes, **kwargs)
                plt.show()
            else:
                plotfun(*args, axes=axes, **kwargs)

        return wrapper

    return actualdecorator

@openfigureifnoaxes(projection=None)
def scatter2d(x, y, c="C0", log=[False, False, False], nc=9, cmap="viridis",
              xlabel=None, ylabel=None, clabel=None, axesequal=False, axes=None,
              cax=None):
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
    c : {str, array_like}, optional
        Color data or string indicating the color.
    log : [bool, bool, bool], optional
        Make [x-axis, y-axis, color axis] logarithmic.
    nc : int, optional
        Number of colors used if c contains data. Since we are using plot
        instead of data, the color scale can't be continuous.
    cmap : str, optional
        Name of the colormap where nc colors are picked if c contains data.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    axesequal : bool, optional
        Flag to set aspect ratio of [x,y] axes equal.
    axes : Axes, optional
        The Axes object to draw on.
    cax : Axes, optional
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
    c : {str, array_like}, optional
        Color data or string indicating the color.
    log : [bool, bool, bool, bool], optional
        Make [x-axis, y-axis, z-axis, color axis] logarithmic.
    nc : int, optional
        Number of colors used if c contains data. Since we are using plot
        instead of data, the color scale can't be continuous.
    cmap : str, optional
        Name of the colormap where nc colors are picked if c contains data.
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
    axes : Axes, optional
        The Axes object to draw on. If None, a new figure is displayed.
    cax : Axes, optional
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

@openfigureifnoaxes(projection=None)
def hist1d(x, xbins=None, weights=None, log=[False, False], xlabel=None,
           legend=None, axes=None):
    """Plot (stacked) marker histogram in 1D.

    Parameters
    ----------
    x : array_like
        Array or a list of arrays to be binned and plotted.

        If list is given, the resulting histogram is stacked with each
        array in the list corresponding to one layer in the stacked histogram.
    xbins : int or array_like, optional
        Number of bins or array storing bin edges for the x coordinate.
    weights : array_like, optional
        Values the datapoints are weighted with.

        Same format as for x.
    log : [bool, bool], optional
        Make [x, y] axes logarithmic.
    xlabel : str, optional
        Label for the x-axis.
    legend : str, array_like
        List of strings to label the data with.

        The length of the list must be same as the number of data arrays
        provided.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    """

    # Set log scales
    if log[0]: x = [np.log10(x) for x in x]
    if log[0] and ( xbins is not None and not isinstance(xbins, int) ):
        xbins = np.log10(xbins)

    if not log[1]: axes.ticklabel_format(style="sci", axis="y", scilimits=(0,0))

    # Plot and legend
    axes.hist(x, xbins, density=False, stacked=True, log=log[1],
              weights=weights, rwidth=1)
    axes.legend(legend, frameon=False)

    # Set labels
    if log[0] and xlabel is not None: xlabel = "log10( "+xlabel+" )"
    axes.set_xlabel(xlabel)

    ylabel = "Markers per bin" if weights is None else "Particles per bin"
    axes.set_ylabel(ylabel)

@openfigureifnoaxes(projection=None)
def hist2d(x, y, xbins=None, ybins=None, weights=None,
           log=[False, False, False], xlabel=None, ylabel=None, axesequal=False,
           axes=None, cax=None):
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
    log : [bool, bool, bool], optional
        Make [x, y, color] axes logarithmic.
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
    # Set log scales
    if log[0]: x = [np.log10(np.abs(x)) for x in x]
    if log[0] and ( xbins is not None and not isinstance(xbins, int) ):
        xbins = np.log10(xbins)

    if log[1]: y = [np.log10(np.abs(y)) for y in y]
    if log[1] and ( ybins is not None and not isinstance(ybins, int) ):
        ybins = np.log10(ybins)

    # Plot and colorbar
    h,_,_,m = axes.hist2d(x, y, bins=[xbins, ybins], weights=weights)

    norm = None
    if log[2]:
        norm = mpl.colors.LogNorm(np.amin(h), np.amax(h))
    cbar = plt.colorbar(m, norm=norm, ax=axes, cax=cax)

    # Set labels
    if log[0] and xlabel is not None: xlabel = "log10( "+xlabel+" )"
    axes.set_xlabel(xlabel)
    if log[1] and ylabel is not None: ylabel = "log10( "+ylabel+" )"
    axes.set_ylabel(ylabel)

    clabel = "Markers per bin" if weights is None else "Particles per bin"
    cbar.set_label(clabel)

    if axesequal:
        axes.set_aspect("equal", adjustable="box")

@openfigureifnoaxes(projection=None)
def mesh2d(x, y, z, xlabel=None, ylabel=None, clim=[None, None], axesequal=False,
         axes=None):
    """Make a mesh (surface) plot in 2D.

    Parameters
    ----------
    x : array_like
    y : array_like
    z : array_like
    """
    axes.image(z)

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    if axesequal:
        axes.set_aspect("equal", adjustable="box")
    """
    out = np.ma.masked_invalid(out)
    if clim[0] is None:
            clim[0] = np.nanmin(out)
        if clim[1] is None:
            clim[1] = np.nanmax(out)

        mesh = axes.pcolormesh(r, z, out, vmin=clim[0], vmax=clim[1])
        axes.patch.set(hatch='x', edgecolor=[0.9, 0.9, 0.9])

        plt.colorbar(mesh, cax=cax)
        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(r[0], r[-1])
        axes.set_ylim(z[0], z[-1])
    """

def contour2d():
    """Plot contour on 2D plane.
    """
    pass

def line2d():
    """Plot line segments on 2D plane.
    """
    pass

@openfigureifnoaxes(projection=None)
def poincare(x, y, ids, conlen=None, xlim=None, ylim=None, xlabel=None,
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
    conlen : array_like, optional
        Connection length at the position (x,y). Negative if the marker is
        confined. If given, the color scale shows connection length instead
        of marker ID. The confined markers are still shown with shades of red.
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
    axes : Axes, optional
        The Axes object to draw on.
    cax : Axes, optional
        The Axes object for the connection length (if conlen is given),
        otherwise taken from axes.
    """
    nc = 5 # How many colors are used

    if conlen is not None:
        # If we don't have lost markers then ignore the connection length.
        if np.argwhere(conlen > 0).size == 0:
            conlen = None

    if conlen is None:
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
        # negative conlen).
        idx = np.argsort(conlen)
        x   = x[idx]
        y   = y[idx]
        ids = ids[idx]
        conlen = conlen[idx]

        # Find where the line between confined and lost markers is.
        # Set conlen positive again for confined markers and rearrange again
        # by the connection length.
        idx    = np.argwhere(conlen > 0).ravel()[0]
        x      = np.append(x[idx:],      np.flip(x[:idx]))
        y      = np.append(y[idx:],      np.flip(y[:idx]))
        ids    = np.append(ids[idx:],    np.flip(ids[:idx]))
        conlen = np.append(conlen[idx:], np.flip(-conlen[:idx]))
        idx = (conlen.size - 1) - idx

        # The color has meaning only for lost markers so find the scale
        cmin = conlen[0]
        cmax = conlen[idx]

        logscale = False
        if cmin / cmax < 0.1:
            logscale = True
            cmin   = np.log10(cmin)
            cmax   = np.log10(cmax)
            conlen = np.log10(conlen)

        cmin = np.floor(cmin)
        cmax = np.ceil(cmax)
        nc_b = int(cmax - cmin)
        clim = np.linspace(cmin, cmax, nc_b+1)

        # Confined markers are plotted separately.
        conlen[idx+1:] = cmax + 1/nc
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
        idx  = np.searchsorted(conlen, clim) + 1

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
    axes : Axes, optional
        The Axes object to draw on.
    cax : Axes, optional
        The Axes object for the color data (if c contains data), otherwise
        taken from axes.
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
