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
    ind_int = np.logical_and(np.abs(np.sum(np.copysign(1, p1p2p3), 1)) < 3,
                             ind_crossing)
    # Remove extra data
    x1x2x3 = x1x2x3[ind_int, :]
    y1y2y3 = y1y2y3[ind_int, :]
    z1z2z3 = z1z2z3[ind_int, :]
    r1r2r3 = r1r2r3[ind_int, :]
    p1p2p3 = p1p2p3[ind_int, :]
    n = x1x2x3.shape[0]
    # Find intersection line for all crossing triangles
    r = np.zeros((n, 2))
    z = np.zeros((n, 2))
    for i in range(n):
        # For the crossing triangles, one of the points is always on the other
        # side of phi-plane as the other two
        polarity = np.sum(np.copysign(1, p1p2p3[i, :]))
        p_diff = np.flatnonzero(np.copysign(1, p1p2p3[i, :]) != polarity)
        # The triangle edge is linear, so we can easily calculate the
        # cross-section points
        p_same = np.setdiff1d([0, 1, 2], p_diff)
        for j in range(2):
            phi_diff = (
                    - p1p2p3[i, p_same[j]]  # Distance to phi = 0
                    / (p1p2p3[i, p_diff] - p1p2p3[i, p_same[j]]) # Width of phi
                    )
            r[i, j] = (r1r2r3[i, p_same[j]]
                       + (r1r2r3[i, p_diff] - r1r2r3[i, p_same[j]]) * phi_diff)
            z[i, j] = (z1z2z3[i, p_same[j]]
                       + (z1z2z3[i, p_diff] - z1z2z3[i, p_same[j]]) * phi_diff)

    axes.plot(r.T, z.T, 'k')
    axes.axis("scaled")

    if newfig:
        plt.show(block=False)
