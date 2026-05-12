"""Helper module for plotting related matters."""

from a5py.utils import OptionalDependency

_plt = OptionalDependency("matplotlib.pyplot")

from .misc import (
    set_style,
    create_figure,
    open_figure_if_no_axes,
    get_ax10b_formatter,
    )
from .plotting import (
    scatter2d,
    scatter3d,
    hist1d,
    hist2d,
    mesh1d,
    mesh2d,
    contour2d,
    arrows2d,
    line2d,
    line3d,
    poincare,
    still,
    interactive,
    triangularpatch,
    )

openfigureifnoaxes = open_figure_if_no_axes

__all__ = [
    "line2d",
    "set_style",
    "create_figure",
    "open_figure_if_no_axes",
    "get_ax10b_formatter",
    ]

def __getattr__(name):
    return getattr(_plt, name)

def __dir__():
    return sorted(__all__ + dir(_plt))
