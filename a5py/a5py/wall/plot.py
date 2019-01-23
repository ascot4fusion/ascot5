"""
Routines for plotting walls.

File: wall/plot.py
"""
import numpy as np

import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
        import matplotlib.pyplot as plt


def plot_segments(x, y, axes=None):
    """
    """
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    x = np.append(x, x[0])
    y = np.append(y, y[0])
    axes.plot(x, y, color="black", linewidth=2)
    axes.axis("scaled")

    if newfig:
        plt.show(block=False)

def plot_projection(x1x2x3, y1y2y3, z1z2z3, axes=None):
    """
    """
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    r1r2r3 = np.sqrt(x1x2x3 * x1x2x3 + y1y2y3 * y1y2y3)

    for i in range(r1r2r3.size):
        axes.plot(r1r2r3[i], z1z2z3[i])
    axes.axis("scaled")

    if newfig:
        plt.show(block=False)
