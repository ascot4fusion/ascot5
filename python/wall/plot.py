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

    if newfig:
        axes.axis("scaled")
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

    if newfig:
        axes.axis("scaled")
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
    # Wrap phi to range [-pi, pi) and set p1p2p3 == phi to zero
    phi = np.mod(phi + np.pi, 2*np.pi) - np.pi
    p1p2p3 = p1p2p3 - phi
    # Mod p1p2p3 again to range [-pi, pi)
    p1p2p3 = np.mod(p1p2p3 + np.pi, 2*np.pi) - np.pi
    # Filter out triangles crossing the -pi/pi boundary
    ind_valid = np.amax(p1p2p3, 1) - np.amin(p1p2p3, 1) < np.pi
    # Intersection happens if any(p1p2p3) < phi = 0 and any(p1p2p3) >= phi = 0
    p1p2p3_sign = np.copysign(1, p1p2p3)
    p1p2p3_polarity = np.sum(p1p2p3_sign, 1)
    ind_inter = np.logical_and(np.abs(p1p2p3_polarity) < 3, ind_valid)

    # Find intersection line for all crossing triangles
    p1p2p3_polarity = np.tile(p1p2p3_polarity, (3,1)).T
    ind_inter = np.tile(ind_inter, (3,1)).T
    ind_diff = np.logical_and(p1p2p3_sign != p1p2p3_polarity, ind_inter)
    ind_same = np.logical_and(p1p2p3_sign == p1p2p3_polarity, ind_inter)

    phi_diff = p1p2p3[ind_diff]
    phi_same = np.reshape(p1p2p3[ind_same], (-1,2))
    phi_diffs = (- phi_same / (phi_diff[:,None] - phi_same))

    r_diff = r1r2r3[ind_diff]
    r_same = np.reshape(r1r2r3[ind_same], (-1,2))
    r = r_same + (r_diff[:,None] - r_same ) * phi_diffs

    z_diff = z1z2z3[ind_diff]
    z_same = np.reshape(z1z2z3[ind_same], (-1,2))
    z = z_same + (z_diff[:,None] - z_same ) * phi_diffs

    # Add triangles for which one side lies entirely on phi = 0
    ind_2zero = np.sum(p1p2p3 == 0, axis=1) == 2
    ind_2zero = np.tile(ind_2zero, (3,1)).T
    ind_zero = p1p2p3 == 0
    ind_valid = np.tile(ind_valid, (3,1)).T
    ind_onplane = np.logical_and.reduce((ind_2zero, ind_zero, ind_valid))

    r2 = np.reshape(r1r2r3[ind_onplane], (-1,2))
    z2 = np.reshape(z1z2z3[ind_onplane], (-1,2))

    # Add triangles for which all 3 sides are entirely on phi = 0
    ind_3zero = np.all(p1p2p3 == 0, axis=1)
    r3 = np.concatenate((r1r2r3[ind_3zero][:,(0,1)],
                         r1r2r3[ind_3zero][:,(1,2)],
                         r1r2r3[ind_3zero][:,(2,0)]))
    z3 = np.concatenate((z1z2z3[ind_3zero][:,(0,1)],
                         z1z2z3[ind_3zero][:,(1,2)],
                         z1z2z3[ind_3zero][:,(2,0)]))

    axes.plot(r.T, z.T, 'k')
    axes.plot(r2.T, z2.T, 'k')
    axes.plot(r3.T, z3.T, 'k')

    if newfig:
        axes.axis("scaled")
        plt.show(block=False)
