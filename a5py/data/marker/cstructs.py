import ctypes

from ...libascot import LIBASCOT

from ..cstructs import ENUM


INPUT_PARTICLE_TYPE = {
    "p":0,
    "gc":1,
    "ml":2,
    "s":3,
}
"""Enum values for the input particle types."""


# pylint: disable=too-few-public-methods
class particle(ctypes.Structure):
    """Python wrapper for the particle input struct in particle.h."""

    _pack_ = 1
    _fields_ = [
        ('r', ctypes.c_double),
        ('phi', ctypes.c_double),
        ('z', ctypes.c_double),
        ('p_r', ctypes.c_double),
        ('p_phi', ctypes.c_double),
        ('p_z', ctypes.c_double),
        ('mass', ctypes.c_double),
        ('charge', ctypes.c_double),
        ('anum', ctypes.c_int32),
        ('znum', ctypes.c_int32),
        ('weight', ctypes.c_double),
        ('time', ctypes.c_double),
        ('id', ctypes.c_int64),
        ]


# pylint: disable=too-few-public-methods
class particle_gc(ctypes.Structure):
    """Python wrapper for the guiding center input struct in particle.h."""

    _pack_ = 1
    _fields_ = [
        ('r', ctypes.c_double),
        ('phi', ctypes.c_double),
        ('z', ctypes.c_double),
        ('energy', ctypes.c_double),
        ('pitch', ctypes.c_double),
        ('zeta', ctypes.c_double),
        ('mass', ctypes.c_double),
        ('charge', ctypes.c_double),
        ('anum', ctypes.c_int32),
        ('znum', ctypes.c_int32),
        ('weight', ctypes.c_double),
        ('time', ctypes.c_double),
        ('id', ctypes.c_int64),
        ]


# pylint: disable=too-few-public-methods
class particle_ml(ctypes.Structure):
    """Python wrapper for the magnetic field line input struct in particle.h."""

    _pack_ = 1
    _fields_ = [
        ('r', ctypes.c_double),
        ('phi', ctypes.c_double),
        ('z', ctypes.c_double),
        ('pitch', ctypes.c_double),
        ('weight', ctypes.c_double),
        ('time', ctypes.c_double),
        ('id', ctypes.c_int64),
    ]


class particle_state(ctypes.Structure):
    """Python wrapper for the particle state struct in particle.h."""

    _pack_ = 1
    _fields_ = [
        ('r', ctypes.c_double),
        ('phi', ctypes.c_double),
        ('z', ctypes.c_double),
        ('ppar', ctypes.c_double),
        ('mu', ctypes.c_double),
        ('zeta', ctypes.c_double),
        ('rprt', ctypes.c_double),
        ('phiprt', ctypes.c_double),
        ('zprt', ctypes.c_double),
        ('p_r', ctypes.c_double),
        ('p_phi', ctypes.c_double),
        ('p_z', ctypes.c_double),
        ('mass', ctypes.c_double),
        ('charge', ctypes.c_double),
        ('anum', ctypes.c_int32),
        ('znum', ctypes.c_int32),
        ('weight', ctypes.c_double),
        ('time', ctypes.c_double),
        ('mileage', ctypes.c_double),
        ('cputime', ctypes.c_double),
        ('rho', ctypes.c_double),
        ('theta', ctypes.c_double),
        ('id', ctypes.c_int64),
        ('endcond', ctypes.c_int64),
        ('walltile', ctypes.c_int64),
        ('B_r', ctypes.c_double),
        ('B_phi', ctypes.c_double),
        ('B_z', ctypes.c_double),
        ('B_r_dr', ctypes.c_double),
        ('B_phi_dr', ctypes.c_double),
        ('B_z_dr', ctypes.c_double),
        ('B_r_dphi', ctypes.c_double),
        ('B_phi_dphi', ctypes.c_double),
        ('B_z_dphi', ctypes.c_double),
        ('B_r_dz', ctypes.c_double),
        ('B_phi_dz', ctypes.c_double),
        ('B_z_dz', ctypes.c_double),
        ('err', ctypes.c_uint64),
        ]


class union_input_particle(ctypes.Union):
    """Union of input particle types and particle state."""

    _pack_ = 1
    _fields_ = [
        ('p', particle),
        ('p_gc', particle_gc),
        ('p_ml', particle_ml),
        ('p_s', particle_state),
        ]


class input_particle(ctypes.Structure):
    """Python wrapper for the input particle struct in particle.h."""

    _pack_ = 1
    _anonymous_ = ('_0',)
    _fields_ = [
        ('type', ENUM),
        ('PADDING_0', ctypes.c_ubyte * 4),
        ('_0', union_input_particle),
        ]


def allocate_input_particles(n:int):
    """Allocate memory for an input particle array.

    Parameters
    ----------
    n : int
        Number of particles to allocate.

    Returns
    -------
    array : ctypes.POINTER(input_particle)
        Pointer to the allocated input particle array.
    """
    fun = LIBASCOT.libascot_allocate_input_particles
    fun.restype = ctypes.POINTER(input_particle)
    fun.argtypes = [ctypes.c_int32]
    return fun(n)


def allocate_particle_state(n:int):
    """Allocate memory for an input particle array.

    Parameters
    ----------
    n : int
        Number of particles to allocate.

    Returns
    -------
    array : ctypes.POINTER(input_particle)
        Pointer to the allocated input particle array.
    """
    fun = LIBASCOT.libascot_allocate_particle_states
    fun.restype = ctypes.POINTER(particle_state)
    fun.argtypes = [ctypes.c_int32]
    return fun(n)
