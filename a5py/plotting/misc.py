"""Miscellaneous routines for plotting."""

from functools import wraps

from a5py.utils import OptionalDependency

mpl = OptionalDependency("matplotlib")
plt = OptionalDependency("matplotlib.pyplot")


def set_style(style: str, latex: bool = True):
    """Set default figure settings (label sizes etc.) to suit given style.

    Parameters
    ----------
    style : {"paper", "presentation"}
        Style to use.
    latex : bool, optional
        Use LaTex interpreter.
    """
    styles = {
        "paper": {
            "figure.autolayout": False,
            "font.family": "serif",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.labelsize": 14,
            "axes.titlesize": 16,
            "axes.titlepad": 6,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.labelpad": 6,
            "legend.fontsize": 12,
            "legend.numpoints": 1,
            "legend.scatterpoints": 1,
            "grid.linewidth": 0.8,
            "lines.linewidth": 1.4,
            "patch.linewidth": 0.24,
            "lines.markersize": 5.6,
            "lines.markeredgewidth": 0,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.pad": 5.6,
            "ytick.major.pad": 5.6,
            "savefig.dpi": 300,
            "axes.formatter.limits": [-2, 2],
        },
        "presentation": {
            "figure.autolayout": False,
            "font.family": "serif",
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "axes.titlepad": 12,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "axes.labelpad": 6,
            "legend.fontsize": 16,
            "legend.numpoints": 1,
            "legend.scatterpoints": 1,
            "grid.linewidth": 0.8,
            "lines.linewidth": 1.4,
            "patch.linewidth": 0.24,
            "lines.markersize": 5.6,
            "lines.markeredgewidth": 0,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.pad": 5.6,
            "ytick.major.pad": 5.6,
            "savefig.dpi": 300,
            "axes.formatter.limits": [-2, 2],
        },
    }
    mpl.style.use({styles[style]})
    if latex:
        mpl.style.use(
            {
                "font.serif": "ComputerModern",
                "text.usetex": True,
            }
        )


def create_figure(*args, aspectratio=3./2, doublecolumn=False, **kwargs):
    """Return figure that has a size suitable for printing in A4 assuming the
    paper has a single-column format.

    Parameters
    ----------
    *args : tuple
        Arguments for ``plt.figure``.
    aspectratio : float
        Width / height ratio of the returned figure.
    doublecolumn : bool
        If True, return figure whose width spans both columns.
    **kwargs : dict
        Keyword arguments for ``plt.figure``.

    Returns
    -------
    """
    if doublecolumn:
        return plt.figure(*args, figsize=(7.205, 7.205 / aspectratio), **kwargs)
    return plt.figure(*args, figsize=(3.504, 3.504 / aspectratio), **kwargs)


def open_figure_if_no_axes(projection: str="rectilinear"):
    """Display new figure if ``axes`` is None and pass *args and **kwargs for
    the plotter.

    Use this as a decorator for plotting routines.

    Parameters
    ----------
    projection : str, optional
        The axes projection type, see Matplotlib documentation for details.

        Available values are:
        "rectilinear", "polar", "3d", "aitoff", "hammer", "lambert", and
        "mollweide".

    Returns
    -------
    decorator : function
        Decorator that takes the plotting routine.
    """

    def actualdecorator(plotfun):
        """Decorator that takes the plotting routine."""

        @wraps(plotfun)
        def create_figure_if_no_axes(*args, axes=None, **kwargs):
            """Create a new figure if ``axes`` is None and pass *args and
            **kwargs for the plotter.
            """
            if axes is None:
                fig = plt.figure()
                axes = fig.add_subplot(111, projection=projection)
                plotfun(*args, axes=axes, **kwargs)
                plt.show()
            else:
                if axes.name != projection:
                    raise ValueError(
                        f"Invalid projection '{axes.name}' on axes: "
                        f"expected '{projection}'"
                    )
                plotfun(*args, axes=axes, **kwargs)

        return create_figure_if_no_axes

    return actualdecorator


def get_ax10b_formatter(format: str) -> mpl.ticker.Formatter:
    """Returns a label tick formatter that shows numbers in format "a x 10^b".

    Credit: https://stackoverflow.com/a/49330649

    Parameters
    ----------
    format : str
        Format string for the numbers.

    Returns
    -------
    formatter : mpl.ticker.Formatter
        A tick formatter that shows numbers in format "a x 10^b".

    Examples
    --------
    >>> plt.gca().yaxis.set_major_formatter(getmathtextsciformatter("%1.2e"))
    """

    class MathTextSciFormatter(mpl.ticker.Formatter):

        def __init__(self, format="%1.2e"):
            self.fmt = format

        def __call__(self, x):
            s = self.fmt % x
            decimal_point = "."
            positive_sign = "+"
            tup = s.split("e")
            significand = tup[0].rstrip(decimal_point)
            sign = tup[1][0].replace(positive_sign, "")
            exponent = tup[1][1:].lstrip("0")
            if exponent:
                exponent = "10^{%s%s}" % (sign, exponent)
            if significand and exponent:
                s = r"%s{\times}%s" % (significand, exponent)
            else:
                s = r"%s%s" % (significand, exponent)
            return "${}$".format(s)

    return MathTextSciFormatter(format)
