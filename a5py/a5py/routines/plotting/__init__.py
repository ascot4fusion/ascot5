"""Package for generating plots with ASCOT5.
"""
from .plothelpers import openfigureifnoaxes
from .scatter     import scatter2d, scatter3d
from .histogram   import hist1d, hist2d
from .poincare    import poincare
from .wall3d      import still, interactive, loadvsarea, defaultcamera
