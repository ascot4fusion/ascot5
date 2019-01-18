"""
Routines for plotting markers.

File: marker/plot.py
"""
import numpy as np
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

def plot_orbit(x, y, ids=None, mask=None, axes=None, equal=False):
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    if ids is not None:
        uids = np.unique(ids)
        for i in uids:
            axes.plot(x[i==ids], y[i==ids])
    else:
        axes.plot(x, y)

    if equal:
        axes.axis("scaled")

    if newfig:
        plt.show(block=False)

def plot_histogram(arg, bins=10):
    pass

def plot_scatter(arg, bins=10):
    pass
