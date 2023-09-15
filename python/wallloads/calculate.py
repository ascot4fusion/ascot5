"""
Created on Mar 9, 2020

@author: sjjamsa
"""

import numpy as np
import a5py.marker.endcond
from scipy.constants import physical_constants as const

def wallLoad3DEndstate(run):
    """
    Calculates the wall load from 3D wall for the specified run.
    @param hdf5: a5py.ascot5io.ascot5.Ascot(filename)

    @return: A float array of power load per triangle (W/m^2).
    """

    cW = run.wall  # The wall class
    es = run.endstate  # The Endstates

    A = cW.area()  # Triangle area
    nWallTris = cW.getNumberOfElements()

    E = es.get("endcond")
    wall_hit_endcond = a5py.marker.endcond.getbin("wall")
    # Pick the wall hit particles
    WH = (E & wall_hit_endcond) > 0

    T = es.get("walltile")[WH] - 1
    weight = es.get("weight")[WH]
    ene = es.get("energy")[WH]
    # The energies are now in eV. Need to convert to Joules.
    ene = ene * const["elementary charge"][0]

    # Calculate the power represented by each marker
    P = np.multiply(weight, ene)

    # Calculate the wall load
    wallLoad = np.bincount(
        T, weights=P, minlength=nWallTris
    )  # This is a histogram, where triangle indexes are the bins
    nonZero = wallLoad > 0.0
    wallLoad[nonZero] = np.divide(wallLoad[nonZero], A[nonZero])

    # This isn't accessed before function return. Should probably be removed
    # but leaving in case it is a useful variable that might be wanted in the
    # future
    # Nmarkers = np.bincount(T, minlength=nWallTris)

    return wallLoad
