"""
Evaluate particle-related physical quantities (that are not done in other files)
"""

import numpy as np

from .gamma import gamma_velocity, gamma_momentum, velocity_momentum

def momentum_muppar(m, mu, ppar, b, zeta=None):
    """
    Evaluate momentum from mu and ppar.

    If gyroangle zeta is given, return the momentum vector. Otherwise return
    momentum norm.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )

    return np.sqrt( 2 * m * mu * b + ppar**2 )


def velocity_muppar(m, mu, ppar, b, zeta=None):
    """
    Evaluate velocity from mu and ppar.

    If gyroangle zeta is given, return the velocity vector. Otherwise return
    velocity norm.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )

    p = np.sqrt( 2 * m * mu * b + ppar**2 )
    return velocity_momentum(m=m, p=p)


def mu_momentum(m, p, b):
    """
    Evaluate magnetic moment from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    pnorm = np.sqrt( np.sum(p**2, axis=0) )
    pitch = np.sum( p * b, axis=0 ) / ( pnorm * bnorm )
    return ( 1 - pitch**2 ) * pnorm**2 / ( 2 * bnorm * m )


def ppar_momentum(p, b):
    """
    Evaluate parallel momentum from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    return np.sum( p * b, axis=0 ) / bnorm


def vpar_momentum(m, p, b):
    """
    Evaluate parallel velocity from momentum vector.
    """
    gamma = gamma_momentum(m=m, p=p)
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    return np.sum( p * b, axis=0 ) / ( gamma * m * bnorm)


def mu_velocity(m, v, b):
    """
    TODO: Evaluate magnetic moment from velocity vector.
    """
    gamma = gamma_velocity(v=v)
    return 0


def pitch_momentum(p, b):
    """
    Evaluate pitch from momentum vector.
    """
    return np.sum( p * b, axis=0 ) / np.sqrt(
        np.sum(p**2, axis=0) * np.sum(b**2, axis=0) )


def pitch_muppar(m, mu, ppar, b):
    """
    Evaluate pitch from mu and ppar.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )
    return ppar / np.sqrt( 2 * m * mu * b + ppar**2 )


def torcanangmom_momentum(q, r, p, psi):
    """
    Evaluate toroidal canonical angular momentum from momentum.
    """
    return p[1,:] * r + q * psi


def torcanangmom_ppar(q, r, ppar, b, psi):
    """
    Evaluate toroidal canonical angular momentum from mu and ppar.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    return ppar * r * b[1,:] / bnorm + q * psi


def gyrolength(m, energy, pitch, b):
    """
    Evaluate gyrolength from energy and pitch.
    """
    pass


def gyrofrequency(m, energy, b):
    """
    Evaluate gyrofrequency from energy.
    """
    pass
