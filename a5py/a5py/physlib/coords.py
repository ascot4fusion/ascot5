"""
Ordinary coordinate transformations between different geometrical systems.
"""
import numpy as np
import unyt

def cart2pol(x, y, z=None):
    """Cartesian x coordinate from cylindrical coordinates.
    """
    return np.sqrt(x**2 + y**2), np.arctan2(y, x) * unyt.rad, z

def pol2cart(r, phi, z=None):
    """Cartesian x coordinate from cylindrical coordinates.
    """
    return r * np.cos(phi), r * np.sin(phi), z


def anglemod(angle):
    """Transfer arbitrary angle to between interval [0, 2pi].
    """
    return np.mod(angle + 2*np.pi, 2*np.pi)
