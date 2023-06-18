import matplotlib.pyplot as plt


def openfigureifnoaxes(projection=None):
    """
    Decorator for creating and displaying a new figure if axes argument is None.

    Args:
        projection : {None, 'aitoff', 'hammer', 'lambert', 'mollweide', 'polar',
        'rectilinear', str}, optional <br>
            The axes projection type, see Matplotlib documentation for details.
    """
    def actualdecorator(plotfun):
        def wrapper(*args, axes=None, **kwargs):
            if not axes:
                fig  = plt.figure()
                axes = fig.add_subplot(111, projection=projection)
                plotfun(*args, axes=axes, **kwargs)
                plt.show()
            else:
                plotfun(*args, axes=axes, **kwargs)

        return wrapper

    return actualdecorator
