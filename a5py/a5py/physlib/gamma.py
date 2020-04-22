"""
Evaluate Lorentz factor and energy from various parameters.
"""

import numpy as np
from unyt import c

def velocity_momentum(m, p):
    """
    Evaluate velocity from momentum.

    If p is a vector (scalar) then returned value is a vector (scalar).
    """
    return p / ( m * gamma_momentum(m=m, p=p) )


def momentum_velocity(m, v):
    """
    Evaluate momentum from velocity.

    If v is a vector (scalar) then returned value is a vector (scalar).
    """
    return gamma_velocity(v=v) * m * v


def gamma_momentum(m, p):
    """
    Evaluate gamma from momentum.
    """
    if p.shape[1] == 3:
        p = np.sum(p**2, axis=1)
    else:
        p = p**2

    return np.sqrt( 1 + p / ( m * c )**2 )


def gamma_velocity(v):
    """
    Evaluate gamma from velocity.
    """
    if v.shape[1] == 3:
        v = np.sum(v**2, axis=1)
    else:
        v = v**2

    return np.sqrt( 1.0 / ( ( 1 - v / c ) * ( 1 + v / c ) ) )


def gamma_muppar(m, mu, ppar, b):
    """
    Evaluate gamma from magnetic moment and parallel momentum.
    """
    if b.shape[1] == 3:
        b = np.sqrt( np.sum(b**2, axis=1) )

    return np.sqrt( 1 + 2 * mu * b / ( m * c**2 ) + ( ppar / ( m * c ) )**2 )


def gamma_muvpar(m, mu, vpar, b):
    """
    Evaluate gamma from magnetic moment and parallel velocity.
    """
    if b.shape[1] == 3:
        b = np.sqrt( np.sum(b**2, axis=1) )

    a1 = 1 + 2 * mu * b / ( m * c**2 )
    a2 = ( 1 - vpar / c ) * ( 1 + vpar * c )
    return np.sqrt( a1 / a2 )


def gamma_energy(m, energy):
    """
    Evaluate gamma from kinetic energy.
    """
    return 1 + energy / (m * c**2)


def energy_momentum(m, p):
    """
    Evaluate kinetic energy from velocity.
    """
    gamma = gamma_momentum(m=m, p=p)
    return energy_gamma(m=m, gamma=gamma)


def energy_velocity(m, v):
    """
    Evaluate kinetic energy from velocity.
    """
    gamma = gamma_velocity(v=v)
    return gamma#energy_gamma(m=m, gamma=gamma)


def energy_muppar(m, mu, ppar, b):
    """
    Evaluate kinetic energy from magnetic moment and parallel momentum.
    """
    gamma = gamma_muppar(m=m, mu=mu, ppar=ppar, b=b)
    return energy_gamma(m=m, gamma=gamma)


def energy_muvpar(m, mu, vpar, b):
    """
    Evaluate kinetic energy from magnetic moment and parallel velocity.
    """
    gamma = gamma_muvpar(m=m, mu=mu, vpar=vpar, b=b)
    return energy_gamma(m=m, gamma=gamma)


def energy_gamma(m, gamma):
    """
    Evaluate kinetic energy from Lorentz factor.
    """
    return (gamma - 1.0) * m * c**2
