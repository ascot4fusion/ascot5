import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .plothelpers import openfigureifnoaxes


@openfigureifnoaxes(projection=None)
def poincare(x, y, ids, conlen=None, xlim=None, ylim=None, xlabel=None,
             ylabel=None, clabel=None, axesequal=False, axes=None, cax=None):
    """
    Poincaré plot where color separates markers or shows the connection length.

    For better performance this function uses matplotlib's plot function instead
    of the scatter function. The only practical difference is that the marker
    the color can't be continuous.

    Args:
        x : array_like <br>
            Orbit x-coordinates.
        y : array_like <br>
            Orbit y-coordinates.
        ids : array_like <br>
            Array of marker IDs showing to which marker the points in x and y
            arrays correspond to.
        conlen : array_like, optional <br>
            Connection length at the position (x,y). Negative if the marker is
            confined. If given, the color scale shows connection length instead
            of marker ID. The confined markers are still shown with shades of
            red.
        xlim : tuple(float, float), optional <br>
            Min and max values for the x-axis.
        ylim : tuple(float, float), optional <br>
            Min and max values for the y-axis.
        xlabel : str, optional <br>
            Label for the x-axis.
        ylabel : str, optional <br>
            Label for the y-axis.
        clabel : str, optional <br>
            Label for the color axis.
        axes : Axes, optional <br>
            The Axes object to draw on.
        cax : Axes, optional <br>
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
