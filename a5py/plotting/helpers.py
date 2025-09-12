import matplotlib as mpl
import matplotlib.pyplot as plt

from functools import wraps

"""General helper functions for plotting related matters."""

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
