"""
Ordinary coordinate transformations between different geometrical systems.
"""

import numpy as np

def xcoord(r, phi, z=0):
    """
    Cartesian x coordinate from cylindrical coordinates.
    """
    return r * np.cos(phi)


def ycoord(r, phi, z=0):
    """
    Cartesian y coordinate from cylindrical coordinates.
    """
    return r * np.sin(phi)


def anglemod(angle):
    """
    Transfer arbitrary angle to between interval [0, 2pi].
    """
    return np.mod(angle + 2*np.pi, 2*np.pi)
