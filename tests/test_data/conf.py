"""Generates test data for all input variants."""
import unyt
import pytest
import numpy as np
from unyt import (T, Wb, V, m, s, deg, eV, rad, e, particles,)

from a5py import Ascot


unyt.define_unit("yomama", 1.e6 * unyt.kg**3)

test_data = {
    "BfieldCartesian": {
        "bxyz": np.ones((3,)) * T,
        "jacobian": np.zeros((3,3)) * T/m,
        "axisrz": (1.0, 0.0) * m,
        "rhoval": 0.5,
        "_psival": 0.5 * Wb/rad,
        },
    "BfieldAnalytical": {
        "coefficients": [2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                         6.200e-03, -1.205e-03, -3.701e-05,  0.000,
                         0.000,      0.000,      0.000,      0.000,
                        -0.155],
        "psiscaling": 100 * Wb/rad,
        "rmajor": 6.2 * m,
        "bphi": 5.3 * T,
        "_axisrz": [6.2, 0.0] * m,
        "_psilimits": [0.0, 1.0] * Wb/rad,
        "_nripple": 18,
        "_rminor": 2.0 * m,
        "_rippledamping": 0.05,
        "_ripplescaling": 0.5,
        },
    "BfieldSpline2D": {
        "rgrid": np.linspace(4.0, 8.0, 10) * m,
        "zgrid": np.linspace(-2.0, 2.0, 20) * m,
        "psi": np.full((10, 20), 1.0) * Wb/rad,
        "bphi": np.full((10, 20), 1.0) * T,
        "_br": np.full((10, 20), 1.0) * T,
        "_bz": np.full((10, 20), 1.0) * T,
        "axisrz": (6.2, 0.0) * m,
        "psilimits": (0.0, 1.0) * Wb/rad,
        },
    "BfieldSpline3D": {
        "rgrid": np.linspace(4.0, 8.0, 10) * m,
        "phigrid": np.linspace(0, 360, 16)[:-1] * deg,
        "zgrid": np.linspace(-2.0, 2.0, 20) * m,
        "psi": np.full((10, 20), 1.0) * Wb/rad,
        "bphi": np.full((10, 15, 20), 1.0) * T,
        "br": np.full((10, 15, 20), 1.0) * T,
        "bz": np.full((10, 15, 20), 1.0) * T,
        "_rgridpsi": np.linspace(4.0, 8.0, 10) * m,
        "_zgridpsi": np.linspace(-2.0, 2.0, 20) * m,
        "axisrz": (6.2, 0.0) * m,
        "psilimits": (0.0, 1.0) * Wb/rad,
        },
    "BfieldStellarator": {
        "rgrid": np.linspace(4.0, 8.0, 10) * m,
        "phigrid": np.linspace(0, 360, 16)[:-1] * deg,
        "zgrid": np.linspace(-2.0, 2.0, 20) * m,
        "psi": np.full((10, 15, 20), 1.0) * Wb/rad,
        "bphi": np.full((10, 15, 20), 1.0) * T,
        "br": np.full((10, 15, 20), 1.0) * T,
        "bz": np.full((10, 15, 20), 1.0) * T,
        "_rgridpsi": np.linspace(4.0, 8.0, 10) * m,
        "_phigridpsi": np.linspace(0, 360, 16)[:-1] * deg,
        "_zgridpsi": np.linspace(-2.0, 2.0, 20) * m,
        "axisgrid": np.linspace(0, 360, 4)[:-1] * deg,
        "axisrz": np.array([[6.2, 0.0], [6.2, 0.0], [6.2, 0.0]]) * m,
        "psilimits": (0.0, 1.0) * Wb/rad,
        },
    "EfieldCartesian": {
        "exyz": np.array([0, 1, 2]) * V/m,
        },
    "PlasmaLinear1D": {
        "species": ["H1", "He4"],
        "rhogrid": np.linspace(0, 1, 3),
        "ni": np.full((3,2), 1e20) * particles/m**3,
        "Ti": np.full(3, 1e3) * eV,
        "ne": np.full(3, 1e20) * particles/m**3,
        "Te": np.full(3, 1e3) * eV,
        "charge": np.array([1, 2]) * e,
        },
    "WallContour2D": {
        "r": np.array([1, 2, 1]) * m,
        "z": np.array([-1, 0, 1]) * m,
        "_flag": np.array([1, 1, 1]),
        "_labels": {"wood": 1},
        },
    "WallTriangular3D": {
        "vertices": np.array([1, 0, -1, 2, 0, 0, 1, 0, 1]) * m,
        "_flag": np.array([1,]),
        "_labels": {"wood": 1},
        },
    "FieldlineMarker": {
        "r": np.full(3, 6.2) * m,
        "z": np.full(3, 0.0) * m,
        "_phi": np.full(3, 0.0) * deg,
        "_direction": np.full(3, 1.0),
        "_time": np.full(3, 0.0) * s,
        "_ids": np.arange(1, 4),
        },
    "ParticleMarker": {
        "species": "H1",
        "charge": np.full(3, 1) * e,
        "r": np.full(3, 6.2) * m,
        "z": np.full(3, 0.0) * m,
        "vr": np.full(3, 1.0e6) * m/s,
        "vz": np.full(3, 0.0) * m/s,
        "vphi": np.full(3, 0.0) * m/s,
        "_weight": np.full(3, 1.0) * particles/s,
        "_phi": np.full(3, 0.0) * deg,
        "_time": np.full(3, 0.0) * s,
        "_ids": np.arange(1, 4),
        },
    "GuidingcenterMarker": {
        "species": "H1",
        "charge": np.full(3, 1) * e,
        "r": np.full(3, 6.2) * m,
        "z": np.full(3, 0.0) * m,
        "ekin": np.full(3, 1.0e3) * eV,
        "pitch": np.full(3, 1.0),
        "_weight": np.full(3, 1.0) * particles/s,
        "_phi": np.full(3, 0.0) * deg,
        "_gyroangle": np.full(3, 0.0) * rad,
        "_time": np.full(3, 0.0) * s,
        "_ids": np.arange(1, 4),
        },
    }
"""Test data for all input variants."""


def create_input(variant, **replace):
    a5 = Ascot()
    parameters = {
        param[1:] if param.startswith("_") else param: value
        for param, value in test_data[variant].items()
        }
    create = getattr(a5.data, f"create_{variant.lower()}")
    for name, val in replace.items():
        parameters[name] = val
    return create(**parameters)
