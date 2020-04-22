"""
Evaluate particle-related physical quantities.
"""

import numpy as np

def momentum_muppar(m, mu, ppar, b, zeta=None):
    """
    Evaluate momentum from mu and ppar.

    If gyroangle zeta is given, return the momentum vector. Otherwise return
    momentum norm.
    """
    if b.size == 3:
        b = np.sqrt( np.sum(b**2) )
    return ( 2 * mu * b + ppar / m ) / c


def mu_momentum(m, p, b):
    """
    Evaluate magnetic moment from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2) )
    pnorm = np.sqrt( np.sum(p**2) )
    pitch = np.sum( p * b ) / ( pnorm * bnorm )
    return ( 1 - pitch**2 ) * pnorm**2 / ( 2 * bnorm )


def mu_velocity(m, v, b):
    """
    Evaluate magnetic moment from velocity vector.
    """
    gamma = gamma_velocity(v=v)
    return (gamma * m * vperp)**2


def pitch_momentum(p, b):
    """
    Evaluate pitch from momentum vector.
    """
    return np.sum( p * b ) / np.sqrt( np.sum(p**2) * np.sum(b**2) )


def pitch_muppar(m, mu, ppar, b):
    """
    Evaluate pitch from mu and ppar.
    """
    if b.size == 3:
        b = np.sqrt( np.sum(b**2) )
    return ppar * m * c / ( 2 * m * mu * b + ppar )


def cptor_momentum():
    """
    Evaluate canonical angular toroidal momentum from momentum.
    """
    return


def cptor_muppar():
    """
    Evaluate canonical toroidal angular momentum from mu and ppar.
    """


def gyrolength(m, energy, pitch, b):
    """
    Evaluate gyrolength from energy and pitch.
    """


def gyrofrequency(m, energy, b):
    """
    Evaluate gyrofrequency from energy.
    """
