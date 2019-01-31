"""
Routines for plotting magnetic field.

File: field/plot.py
"""
import numpy as np
from numpy.random import randint
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt

def plot(r, z, val, axes=None):

    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    mesh = axes.pcolormesh(r, z, val)
    axes.axis("scaled")

    if newfig:
        plt.show(block=False)

    return axes
