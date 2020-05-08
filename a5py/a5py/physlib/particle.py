"""
Evaluate particle-related physical quantities (that are not done in other files)
"""

import numpy as np

from .gamma import gamma_velocity, gamma_momentum

def momentum_muppar(m, mu, ppar, b, zeta=None):
    """
    Evaluate momentum from mu and ppar.

    If gyroangle zeta is given, return the momentum vector. Otherwise return
    momentum norm.
    """
    if b.shape[1] == 3:
        b = np.sqrt( np.sum(b**2, axis=1) )
    return ( 2 * mu * b + ppar / m ) / c


def mu_momentum(m, p, b):
    """
    Evaluate magnetic moment from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=1) )
    pnorm = np.sqrt( np.sum(p**2, axis=1) )
    pitch = np.sum( p * b, axis=1 ) / ( pnorm * bnorm )
    return ( 1 - pitch**2 ) * pnorm**2 / ( 2 * bnorm * m )


def ppar_momentum(p, b):
    """
    Evaluate parallel momentum from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=1) )
    return np.sum( p * b, axis=1 ) / bnorm


def vpar_momentum(m, p, b):
    """
    Evaluate parallel velocity from momentum vector.
    """
    gamma = gamma_momentum(m=m, p=p)
    bnorm = np.sqrt( np.sum(b**2, axis=1) )
    return np.sum( p * b, axis=1 ) / ( gamma * m * bnorm)


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
    return np.sum( p * b, axis=1 ) / np.sqrt(
        np.sum(p**2, axis=1) * np.sum(b**2, axis=1) )


def pitch_muppar(m, mu, ppar, b):
    """
    Evaluate pitch from mu and ppar.
    """
    if b.shape[1] == 3:
        b = np.sqrt( np.sum(b**2, axis=1) )
    return ppar * m * c / ( 2 * m * mu * b + ppar )


def torcanangmom_momentum(q, r, p, psi):
    """
    Evaluate toroidal canonical angular momentum from momentum.
    """
    return p[:,1] * r + q * psi


def torcanangmom_ppar(q, r, ppar, b, psi):
    """
    Evaluate toroidal canonical angular momentum from mu and ppar.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=1) )
    return ppar * r * b[:,1] / bnorm + q * psi


def gyrolength(m, energy, pitch, b):
    """
    Evaluate gyrolength from energy and pitch.
    """


def gyrofrequency(m, energy, b):
    """
    Evaluate gyrofrequency from energy.
    """
