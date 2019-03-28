"""
Routines for plotting walls.

File: wall/plot.py
"""
import numpy as np

import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
        import matplotlib.pyplot as plt

import a5py.postprocessing.mathlib as mathlib

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

    r1r2r3 = r1r2r3.ravel()
    z1z2z3 = z1z2z3.ravel()
    axes.plot(r1r2r3, z1z2z3, marker=".", markeredgecolor="black",
              linestyle="None")
    axes.axis("scaled")

    if newfig:
        plt.show(block=False)

def plot_intersection(x1x2x3, y1y2y3, z1z2z3, phi, axes=None):
    """
    """
    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    # Polar coordinates, p1p2p3 in range [-pi, pi]
    r1r2r3, p1p2p3 = mathlib.cart2pol(x1x2x3, y1y2y3)
    p1p2p3 = np.degrees(p1p2p3)
    # Intersection happens if
    # phi_i < phi and phi_j >= phi for some i, j = [1 2 3]
    # Phi to range [-180, 180)
    phi = np.degrees(phi)
    phi = np.mod(phi + 180, 360) - 180
    # Set p1p2p3 == phi to zero
    p1p2p3 = p1p2p3 - phi
    ind_crossing = np.amax(p1p2p3, 1) - np.amin(p1p2p3, 1) < 180
    ind_int = np.logical_and.reduce((np.sum(p1p2p3 < 0, 1) < 3,
                                     np.sum(p1p2p3 >= 0, 1) < 3,
                                     ind_crossing))

    plot_projection(x1x2x3[ind_int, :], y1y2y3[ind_int, :], z1z2z3[ind_int, :],
                    axes)
