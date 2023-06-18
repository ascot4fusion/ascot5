import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .plothelpers import openfigureifnoaxes

@openfigureifnoaxes(projection=None)
def hist1d(x, xbins=None, weights=None, log=[False, False], xlabel=None,
           legend=None, axes=None):
    """
    Plot (stacked) histogram in 1D.

    Args:
        x : array_like <br>
            Array or list of arrays of data to be binned and plotted.
        xbins : int or array_like, optional <br>
            Number of bins or array storing bin edges for the x coordinate.
        weights : array_like, optional <br>
            Values the datapoints are weighted with. Same format as for x.
        log : [bool, bool], optional <br>
            Make [x-axis, y-axis] logarithmic.
        xlabel : str, optional <br>
            Label for the x-axis.
        legend : str, array_like <br>
            List of strings to label the data with. The length of the list must
            be same as the number of data arrays provided.
        axes : Axes, optional <br>
            The Axes object to draw on.
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
    """
    Plot histogram in 2D.

    Args:
        x : array_like <br>
            x-coordinates of the data to be binned and plotted.
        y : array_like <br>
            y-coordinates of the data to be binned and plotted.
        xbins : int or array_like, optional <br>
            Number of bins or array storing bin edges for the x coordinate.
        ybins : int or array_like, optional <br>
            Number of bins or array storing bin edges for the y coordinate.
        weights : array_like, optional <br>
            Values the datapoints are weighted with. Same size as x and y.
        log : [bool, bool, bool], optional <br>
            Make [x-axis, y-axis, color axis] logarithmic.
        xlabel : str, optional <br>
            Label for the x-axis.
        ylabel : str, optional <br>
            Label for the y-axis.
        axesequal : bool, optional <br>
            Flag to set aspect ratio of [x,y] axes equal.
        axes : Axes, optional <br>
            The Axes object to draw on.
        cax : Axes, optional <br>
            The Axes object for the color bar, otherwise taken from axes.
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
