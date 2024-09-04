"""Tools for assessing physical quantities and units.
"""
import unyt
import numpy as np
import scipy.constants as constants
import warnings
from typing import Union, Tuple, List, Dict, Any
from functools import wraps

from ..utils import Numerical
from a5py.exceptions import AscotUnitWarning

e     = constants.elementary_charge * unyt.C
m_e   = constants.physical_constants["electron mass"][0] * unyt.kg
m_p   = constants.physical_constants["proton mass"][0] * unyt.kg
m_a   = constants.physical_constants["alpha particle mass"][0] * unyt.kg
c     = unyt.c
eps_0 = unyt.eps_0

def cart2pol(x, y, z=None):
    """Cartesian x coordinate from cylindrical coordinates.
    """
    return np.sqrt(x**2 + y**2), np.arctan2(y, x) * unyt.rad, z

def pol2cart(r, phi, z=None):
    """Cartesian x coordinate from cylindrical coordinates.
    """
    return r * np.cos(phi), r * np.sin(phi), z

def pol2cart_vec(vr,r,vphi,phi,vz=None,z=None):
    vx = vr * np.cos(phi) - vphi * np.sin(phi)
    vy = vr * np.sin(phi) + vphi * np.cos(phi)
    mps = unyt.m / unyt.s
    if vz is not None:
        return vx*mps,vy*mps,vz*mps
    return vx*mps,vy*mps,vz

def cart2pol_vec(vx,x,vy,y,vz=None,z=None):
    phi = np.arctan2( y, x )
    vr  =  vx * np.cos(phi) + vy * np.sin(phi)
    vphi= -vx * np.sin(phi) + vy * np.cos(phi)
    mps = unyt.m / unyt.s
    if vz is not None:
        return vr*mps,vphi*mps,vz*mps
    return vr*mps,vphi*mps,vz

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


def match_units(
        value: Numerical,
        unit: Union[str, unyt.Unit],
        strip: bool = False,
        auto_assign: bool = True,
        ) -> np.ndarray:
    """Check that the units of the value are what is expected.

    In case the value is given without units, it is assumed to be in the
    expected units already.

    Parameters
    ----------
    value : ndarray
        Numerical data with or without units.
    unit : str or unyt.Unit
        Expected units.
    strip : bool, optional
        After checking and converting the value to expected units, strip units
        and return an unitless nd.array instead.
    auto_assign : bool, optional
        Automatically assign the expected units to the value if it does not
        have any units.

    Returns
    -------
    matched : ndarray
        The value in expected units.

    Raises
    ------
    ValueError
        If the units of the value does not match the expected.
    """
    if isinstance(unit, str):
        unit = unyt.unyt_quantity.from_string(unit).units
    if not hasattr(value, "units"):
        if auto_assign:
            value *= unit
        else:
            raise ValueError(
                f"Value has incorrect units: expected {unit} but got "
                f"dimensionless"
            )
    if value.units.dimensions != unit.dimensions:
        raise ValueError(
            f"Value has incorrect dimensions: expected {unit} but got "
            f"{value.units.dimensions}"
        )
    # Converting integer quantities might fail silently (bug?), so we use floats
    value = value.astype("f8", copy=False)
    if strip:
        return value.to(unit).v
    return value.to(unit)


def match_units_arguments(
        args: List[np.ndarray],
        kwargs: Dict[str, np.ndarray],
        argument_names: List[str],
        expected: Dict[str, str],
        strip: bool,
        ) -> Tuple[List[unyt.unyt_array], Dict[str, unyt.unyt_array],
                   Dict[str, str]]:
    """Check that the function arguments have the expected units.

    This function is mostly meant to be used internally. Consider using the
    `parseunits` decorator instead.

    Parameters
    ----------
    args : List[np.ndarray]
        Positional arguments passed to the function.
    kwargs : Dict[str, np.ndarray]
        Keyword arguments passed to the function.
    argument_names : List[str]
        Names of args or all arguments (as they appear in the function
        signature).
    expected : Dict[str, str]
        Name of the argument and the expected unit.
    strip : bool
        After checking and converting the value to expected units, strip units
        and return an unitless nd.array instead.

    Returns
    -------
    parsed_args : Union[List[unyt.unyt_array], List[ndarray]]
        `*args` with the expected units.
    parsed_kwargs : Union[Dict[str, unyt.unyt_array], Dict[str, ndarray]]
        `**kwargs` with the expected units.
    assumedunits : Dict[str, str]
        Names of the arguments, that had no units, and which units they were
        assumed to be in (i.e. the expected units).
    """
    assumedunits = {}
    parsed_args = []
    for arg, name in zip(args, argument_names):
        if name in expected:
            unit = expected[name]
            parsed_args.append(match_units(arg, unit, strip=strip))
            if not hasattr(arg, "units"):
                assumedunits[name] = unit
        else:
            parsed_args.append(arg)
    parsed_kwargs = {}
    for name, arg in kwargs.items():
        if name in expected:
            unit = expected[name]
            parsed_kwargs[name] = match_units(arg, unit, strip=strip)
            if not hasattr(arg, "units"):
                assumedunits[name] = unit
        else:
            parsed_kwargs[name] = arg
    return parsed_args, parsed_kwargs, assumedunits


def parseunits(strip: bool = False, **expected: str) -> Any:
    """Check the units of the arguments before they are passed to the wrapped
    function.

    Parameters
    ----------
    strip : bool, optional
        Strip units after they have been checked so that the wrapped function
        receives unitless arguments (but implicitly in correct units).
    **expected : dict[str, str]
        Name of the argument and the expected unit, e.g. `r="m"`.
    """
    def actualdecorator(func):
        """Check that the arguments that are expected to have units have the
        correct units.

        If an argument is dimensionless, but it is expected to have units, the
        units are assigned and a warning is displayed.

        Parameters
        ----------
        fun : callable
            The wrapped function.
        """
        argument_names = list(func.__code__.co_varnames)
        @wraps(func)
        def wrapper(*args, **kwargs):
            """Wrapper that parses the units of the arguments before they are
            passed.
            """
            args, kwargs, assumed_units = match_units_arguments(
                args, kwargs, argument_names, expected, strip
                )
            if assumed_units:
                names, units = "", ""
                for name, unit in assumed_units.items():
                    if unit == "1":
                        # No need to warn if the quantity is dimensionless
                        continue
                    names += f"'{name}', "
                    units += f"'{unit}', "
                if names:
                    names, units = names[:-2], units[:-2]
                    warnings.warn(
                        f"Argument(s) {names} given without dimensions "
                        f"(assumed {units})", AscotUnitWarning, stacklevel=2,
                        )
            return func(*args, **kwargs)
        return wrapper
    return actualdecorator
