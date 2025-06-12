"""Tools for assessing physical quantities and units.
"""
import unyt
import numpy as np
import scipy.constants as constants
import warnings
from functools import wraps

from . import species
from a5py.exceptions import AscotUnitWarning

e     = constants.elementary_charge * unyt.C
m_e   = constants.physical_constants["electron mass"][0] * unyt.kg
m_p   = constants.physical_constants["proton mass"][0] * unyt.kg
m_a   = constants.physical_constants["alpha particle mass"][0] * unyt.kg
c     = unyt.c
eps_0 = unyt.eps_0

def cart2pol(x, y, z=None):
    """Convert a point in cartesian coordinates to polar coordinates.
    """
    return np.sqrt(x**2 + y**2), np.arctan2(y, x) * unyt.rad, z

def pol2cart(r, phi, z=None):
    """Convert a point in polar coordinates to cartesian coordinates.
    """
    if z is None:
        return r * np.cos(phi), r * np.sin(phi)
    return r * np.cos(phi), r * np.sin(phi), z

def sph2cart(r, phi, theta):
    """Convert a point in spherical coordinates to cartesian coordinates.
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def pol2cart_vec(vr, r, vphi, phi, vz=None, z=None):
    vx = vr * np.cos(phi) - vphi * np.sin(phi)
    vy = vr * np.sin(phi) + vphi * np.cos(phi)
    if vz is not None:
        return vx ,vy, vz
    return vx, vy

def cart2pol_vec(vx, x, vy, y, vz=None, z=None):
    phi = np.arctan2(y, x)
    vr  =  vx * np.cos(phi) + vy * np.sin(phi)
    vphi= -vx * np.sin(phi) + vy * np.cos(phi)
    if vz is not None:
        return vr, vphi, vz
    return vr, vphi

def anglemod(angle):
    """Transfer arbitrary angle to between interval [0, 2pi].
    """
    return np.mod(angle + 2*np.pi, 2*np.pi)

def velocity_momentum(m, p):
    """Evaluate velocity from momentum.

    If p is a vector (scalar) then returned value is a vector (scalar).
    """
    a = m * gamma_momentum(m=m, p=p)
    if p.shape[0] == 3:
        a = np.tile( a, (3,1))

    return p / a

def momentum_velocity(m, v):
    """Evaluate momentum from velocity.

    If v is a vector (scalar) then returned value is a vector (scalar).
    """
    return gamma_velocity(v=v) * m * v

def vpar_muppar(m, mu, ppar, b):
    """Evaluate vpar from mu and ppar.
    """
    gamma = gamma_muppar(m, mu, ppar, b)
    return ppar / (gamma * m)

def gamma_momentum(m, p):
    """Evaluate gamma from momentum.
    """
    if len(p.shape) > 0 and p.shape[0] == 3:
        p = np.sum(p**2, axis=0)
    else:
        p = p**2

    return np.sqrt( 1 + p / ( m * c )**2 )

def gamma_velocity(v):
    """Evaluate gamma from velocity.
    """
    if v.shape[0] == 3:
        v = np.sqrt(np.sum(v**2, axis=0))
    else:
        v = v
    return np.sqrt( 1.0 / ( ( 1 - v / c ) * ( 1 + v / c ) ) )

def gamma_muppar(m, mu, ppar, b):
    """Evaluate gamma from magnetic moment and parallel momentum.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )

    return np.sqrt( 1 + 2 * mu * b / ( m * c**2 ) + ( ppar / ( m * c ) )**2 )

def gamma_muvpar(m, mu, vpar, b):
    """
    Evaluate gamma from magnetic moment and parallel velocity.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )

    a1 = 1 + 2 * mu * b / ( m * c**2 )
    a2 = ( 1 - vpar / c ) * ( 1 + vpar * c )
    return np.sqrt( a1 / a2 )

def gamma_energy(m, energy):
    """Evaluate gamma from kinetic energy.
    """
    return 1 + energy / (m * c**2)

def energy_momentum(m, p):
    """Evaluate kinetic energy from momentum.
    """
    gamma = gamma_momentum(m=m, p=p)
    return energy_gamma(m=m, gamma=gamma)

def energy_velocity(m, v):
    """Evaluate kinetic energy from velocity.
    """
    gamma = gamma_velocity(v=v)
    return energy_gamma(m=m, gamma=gamma)

def energy_muppar(m, mu, ppar, b):
    """Evaluate kinetic energy from magnetic moment and parallel momentum.
    """
    gamma = gamma_muppar(m=m, mu=mu, ppar=ppar, b=b)
    return energy_gamma(m=m, gamma=gamma)

def energy_muvpar(m, mu, vpar, b):
    """Evaluate kinetic energy from magnetic moment and parallel velocity.
    """
    gamma = gamma_muvpar(m=m, mu=mu, vpar=vpar, b=b)
    return energy_gamma(m=m, gamma=gamma)

def energy_gamma(m, gamma):
    """Evaluate kinetic energy from Lorentz factor.
    """
    return (gamma - 1.0) * m * c**2

def vnorm_gamma(gamma):
    """Evaluate velocity norm from Lorentz factor.
    """
    return np.sqrt(1 - 1.0 / gamma**2) * c

def pnorm_gamma(m, gamma):
    """Evaluate momentum norm from Lorentz factor.
    """
    return np.sqrt(gamma**2 - 1) * m * c

def momentum_muppar(m, mu, ppar, b, zeta=None):
    """Evaluate momentum from mu and ppar.

    If gyroangle zeta is given, return the momentum vector. Otherwise return
    momentum norm.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )

    return np.sqrt( 2 * m * mu * b + ppar**2 )

def velocity_muppar(m, mu, ppar, b, zeta=None):
    """Evaluate velocity from mu and ppar.

    If gyroangle zeta is given, return the velocity vector. Otherwise return
    velocity norm.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )

    p = np.sqrt( 2 * m * mu * b + ppar**2 )
    return velocity_momentum(m=m, p=p)

def mu_momentum(m, p, b):
    """Evaluate magnetic moment from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    pnorm = np.sqrt( np.sum(p**2, axis=0) )
    pitch = np.sum( p * b, axis=0 ) / ( pnorm * bnorm )
    return ( 1 - pitch**2 ) * pnorm**2 / ( 2 * bnorm * m )

def ppar_momentum(p, b):
    """Evaluate parallel momentum from momentum vector.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    return np.sum( p * b, axis=0 ) / bnorm

def vpar_momentum(m, p, b):
    """Evaluate parallel velocity from momentum vector.
    """
    gamma = gamma_momentum(m=m, p=p)
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    return np.sum( p * b, axis=0 ) / ( gamma * m * bnorm)

def pitch_momentum(p, b):
    """Evaluate pitch from momentum vector.
    """
    return np.sum( p * b, axis=0 ) / np.sqrt(
        np.sum(p**2, axis=0) * np.sum(b**2, axis=0) )

def pitch_muppar(m, mu, ppar, b):
    """Evaluate pitch from mu and ppar.
    """
    if b.shape[0] == 3:
        b = np.sqrt( np.sum(b**2, axis=0) )
    return ppar / np.sqrt( 2 * m * mu * b + ppar**2 )

def torcanangmom_momentum(q, r, p, psi):
    """Evaluate toroidal canonical angular momentum from momentum.
    """
    return p[1,:] * r + q * psi

def torcanangmom_ppar(q, r, ppar, b, psi):
    """Evaluate toroidal canonical angular momentum from mu and ppar.
    """
    bnorm = np.sqrt( np.sum(b**2, axis=0) )
    return ppar * r * b[1,:] / bnorm + q * psi

def gyrolength(m, q, energy, pitch, bnorm):
    """Evaluate gyrolength from energy and pitch.
    """
    gamma = gamma_energy(m, energy)
    vnorm = vnorm_gamma(gamma)
    return ( gamma * np.sqrt(1 - pitch**2) * m * vnorm
             / (bnorm * np.abs(q)) ).to("m")

def gyrofrequency(m, q, energy, bnorm):
    """Evaluate gyrofrequency from energy.
    """
    gamma = gamma_energy(m, energy)
    return unyt.rad * np.abs(q) * bnorm / (gamma * m)

def bouncefrequency(m, ekin, minorradius, majorradius, safetyfactor):
    """Estimate bounce frequency from energy and aspect ratio.
    """
    eps = minorradius / majorradius
    gamma = gamma_energy(m, ekin)
    vnorm = vnorm_gamma(gamma)
    return np.sqrt(0.5 * eps) * vnorm / ( safetyfactor * majorradius )

def collfreq_ei(mi, qi, ne, Te, clog):
    """Evaluate electron-ion collision frequency.
    """
    return ( ( np.sqrt(2 / np.pi) / 3 )
             * ( unyt.e * qi / ( 4 * np.pi * unyt.eps_0 ) )**2 * ne * clog
             * 4 * np.pi / np.sqrt( unyt.me * Te**3 ) ).to("1/s")

def collfreq_ie(mi, qi, ne, Te, clog):
    """Evaluate ion-electron collision frequency.
    """
    return ( (unyt.me / mi) * ( np.sqrt(2 / np.pi) / 3 )
             * ( unyt.e * qi / ( 4 * np.pi * unyt.eps_0 ) )**2 * ne * clog
             * 4 * np.pi / np.sqrt( unyt.me * Te**3 ) ).to("1/s")

def parseunits(strip=False, **units):
    """Prepare arguments that are expected to have physical units.

    This decorator:

    - Makes sure every argument has expected physical dimensions.
    - Assigns units if units were not provided but they were expected.
    - Strips units if asked (after checking/assignment).

    Examples:

    .. code-block:: python

       @parseunits(x="m", strip=True)
       def fun(x):
           pass

    Parameters
    ----------
    strip : bool, optional
        Strip units so that dimensionless quantities (but in expected units)
        are passed to the function.
    **units
        Argument in wrapped function and the expected unit as a string.
    """
    for k in units.keys():
        # Convert unit strings to unyts
        units[k] = unyt.unyt_quantity.from_string(units[k]).units

    def actualdecorator(fun):
        """Parse fun arguments that are expected to have units.
        """
        argnames = list(fun.__code__.co_varnames)
        def checkandstrip(val, unit, name, assignedunits):
            """Check units of a given argument and strip if needed.
            """
            dim = unit.dimensions
            try:
                # Try to get argument units
                valdim = val.units.dimensions
            except AttributeError:
                # Argument doesn't have units, assign and add warning
                # except if the expected units were "dimensionless"
                val = val*unit
                valdim = dim
                if dim != 1:
                    assignedunits[name] = unit

            if valdim != dim:
                raise ValueError(
                    "\"%s\" has incorrect dimensions: expected %s but got %s" %
                    (name, dim, valdim))
            if strip:
                val = val.astype("f8", copy=False)
                return val.to(unit).v

            # Integers are not converted properly otherwise
            val = val.astype("f8", copy=False)
            return val.to(unit)

        @wraps(fun)
        def wrapper(*args, **kwargs):
            """Replace args and kwargs with parsed arguments when necessary.
            """
            parsedargs = [None]*len(args)
            assignedunits = {} # Dimensionless units with assigned units
            for i, name in enumerate(argnames):
                if i == len(args): break
                if name in units:
                    parsedargs[i] = checkandstrip(args[i], units[name],
                                                  name, assignedunits)
                else:
                    parsedargs[i] = args[i]
            for name, val in kwargs.items():
                if name in units:
                    kwargs[name] = checkandstrip(kwargs[name], units[name],
                                                 name, assignedunits)
            if len(assignedunits) > 0:
                msg1 = "Argument(s) "
                msg2 = " given without dimensions (assumed "
                for name, unit in assignedunits.items():
                    msg1 += name + ", "
                    msg2 += str(unit) + ", "
                msg = msg1[:-2] + msg2[:-2] + ")"
                warnings.warn(msg, AscotUnitWarning, stacklevel=2)

            return fun(*parsedargs, **kwargs)

        return wrapper


    return actualdecorator
