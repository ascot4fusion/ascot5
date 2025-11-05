"""Histogram diagnostics.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin

HIST_R = 0
HIST_PHI = 1
HIST_Z = 2
HIST_RHO = 3
HIST_THETA = 4
HIST_PPAR = 5
HIST_PPERP = 6
HIST_PR = 7
HIST_PPHI = 8
HIST_PZ = 9
HIST_EKIN = 10
HIST_XI = 11
HIST_MU = 12
HIST_PTOR = 13
HIST_TIME = 14
HIST_CHARGE = 15
HIST_NDIM = 16

class HistAxis(ctypes.Structure):
    _fields_ = [
        ("min", ctypes.c_double),
        ("max", ctypes.c_double),
        ("n", ctypes.c_size_t),
        ("name", ctypes.c_int),
    ]

# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    _fields_ = [
        ("bins", ctypes.POINTER(ctypes.c_double)),
        ("nbin", ctypes.c_size_t),
        ("strides", ctypes.c_size_t * (HIST_NDIM - 1)),
        ("axes", HistAxis * HIST_NDIM),
    ]


def init(coordinate, binmin, binmax, nbin):
    hist = Struct()
    hist.axes[0].name = HIST_R
    hist.axes[1].name = HIST_PHI
    hist.axes[2].name = HIST_Z
    hist.axes[3].name = HIST_RHO
    hist.axes[4].name = HIST_THETA
    hist.axes[5].name = HIST_PPAR
    hist.axes[6].name = HIST_PPERP
    hist.axes[7].name = HIST_PR
    hist.axes[8].name = HIST_PPHI
    hist.axes[9].name = HIST_PZ
    hist.axes[10].name = HIST_EKIN
    hist.axes[11].name = HIST_XI
    hist.axes[12].name = HIST_MU
    hist.axes[13].name = HIST_PTOR
    hist.axes[14].name = HIST_TIME
    hist.axes[15].name = HIST_CHARGE

    hist.nbin = 1
    for i in np.range(HIST_NDIM).flip():

        hist.axes[i].n = 0
        hist.axes[i].min = 0
        hist.axes[i].max = 1
        for coord, bmax, bmin, nb in zip(coordinate, binmin, binmax, nbin):
            if hist.axes[i].name == coord:
                hist.axes[i].min = bmin
                hist.axes[i].max = bmax
                hist.axes[i].n = nb
                hist.nbin *= nb
        if i < HIST_NDIM - 1:
            hist.strides[i] = max(hist.axes[i + 1].n, 1)
            if i < HIST_NDIM - 2:
                hist.strides[i] *= hist.strides[i + 1]
    hist.bins = np.zeros()

