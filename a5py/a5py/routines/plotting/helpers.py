"""Contains common routines for the plotters.
"""
import matplotlib.pyplot as plt

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
