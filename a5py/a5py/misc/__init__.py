def openfigureifnoaxes(plotfun):
    import matplotlib.pyplot as plt
    def wrapper(*args, axes=None, **kwargs):
        if not axes:
            fig  = plt.figure()
            axes = fig.add_subplot(111)
            plotfun(*args, axes=axes, **kwargs)
            plt.show()
        else:
            plotfun(*args, axes=axes, **kwargs)

    return wrapper
