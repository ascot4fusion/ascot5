"""Simulation options input.
"""
import h5py
import numpy as np
import ast
import textwrap
import json
import warnings
import xmlschema

import xml.etree.ElementTree as ET

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class Opt(DataGroup):
    """Simulation options.

    Option parameters and default values are defined with properties and
    getters. These should not be accessed directly or altered in any way. Use
    :meth:`read` to obtain options from the HDF5 file.
    """

    def __init__(self, *args, **kwargs):
        """Initialize options attributes with default values.

        This class uses a hack that allows us to do two things simultaneously:

        1. Access and print the option parameters' meta data (name, description,
           and default value) when running the code.
        2. Show that meta data in documentation generated with Sphinx/autodoc.

        The hack is to implement all options attributes as properties even
        though they are not ever used other than to initialize default
        parameters. When adding new options, make sure both the attribute and
        the property is private (to discourage user calling them), and add
        "OPT" prefix for the attribute which allows us to parse options
        attributes from the rest.

        Note that the order in which attributes are declared is preserved.
        """
        super().__init__(*args, **kwargs)

        self._OPT_SIM_MODE                   = 2
        self._OPT_ENABLE_ADAPTIVE            = 1
        self._OPT_RECORD_MODE                = 0
        self._OPT_FIXEDSTEP_USE_USERDEFINED  = 0
        self._OPT_FIXEDSTEP_USERDEFINED      = 1.0e-8
        self._OPT_FIXEDSTEP_GYRODEFINED      = 20
        self._OPT_ADAPTIVE_TOL_ORBIT         = 1.0e-8
        self._OPT_ADAPTIVE_TOL_CCOL          = 1.0e-1
        self._OPT_ADAPTIVE_MAX_DRHO          = 0.1
        self._OPT_ADAPTIVE_MAX_DPHI          = 2.0
        self._OPT_ENDCOND_SIMTIMELIM         = 0
        self._OPT_ENDCOND_CPUTIMELIM         = 0
        self._OPT_ENDCOND_RHOLIM             = 0
        self._OPT_ENDCOND_ENERGYLIM          = 0
        self._OPT_ENDCOND_WALLHIT            = 0
        self._OPT_ENDCOND_MAXORBS            = 0
        self._OPT_ENDCOND_NEUTRALIZED        = 0
        self._OPT_ENDCOND_IONIZED            = 0
        self._OPT_ENDCOND_LIM_SIMTIME        = 1.0
        self._OPT_ENDCOND_MAX_MILEAGE        = 1.0
        self._OPT_ENDCOND_MAX_CPUTIME        = 3600.0
        self._OPT_ENDCOND_MAX_RHO            = 2.0
        self._OPT_ENDCOND_MIN_RHO            = 0.0
        self._OPT_ENDCOND_MIN_ENERGY         = 1.0e3
        self._OPT_ENDCOND_MIN_THERMAL        = 2.0
        self._OPT_ENDCOND_MAX_TOROIDALORBS   = 100
        self._OPT_ENDCOND_MAX_POLOIDALORBS   = 100
        self._OPT_ENABLE_ORBIT_FOLLOWING     = 0
        self._OPT_ENABLE_COULOMB_COLLISIONS  = 0
        self._OPT_ENABLE_MHD                 = 0
        self._OPT_ENABLE_ATOMIC              = 0
        self._OPT_DISABLE_FIRSTORDER_GCTRANS = 0
        self._OPT_DISABLE_ENERGY_CCOLL       = 0
        self._OPT_DISABLE_PITCH_CCOLL        = 0
        self._OPT_DISABLE_GCDIFF_CCOLL       = 0
        self._OPT_REVERSE_TIME               = 0
        self._OPT_ENABLE_DIST_5D             = 0
        self._OPT_ENABLE_DIST_6D             = 0
        self._OPT_ENABLE_DIST_RHO5D          = 0
        self._OPT_ENABLE_DIST_RHO6D          = 0
        self._OPT_ENABLE_RF                  = 0
        self._OPT_DIST_MIN_R                 = 0.1
        self._OPT_DIST_MAX_R                 = 10.0
        self._OPT_DIST_NBIN_R                = 1
        self._OPT_DIST_MIN_PHI               = 0.0
        self._OPT_DIST_MAX_PHI               = 360.0
        self._OPT_DIST_NBIN_PHI              = 1
        self._OPT_DIST_MIN_Z                 = -6.0
        self._OPT_DIST_MAX_Z                 = 6.0
        self._OPT_DIST_NBIN_Z                = 1
        self._OPT_DIST_MIN_RHO               = 0.0
        self._OPT_DIST_MAX_RHO               = 1.0
        self._OPT_DIST_NBIN_RHO              = 1
        self._OPT_DIST_MIN_THETA             = 0.0
        self._OPT_DIST_MAX_THETA             = 360.0
        self._OPT_DIST_NBIN_THETA            = 1
        self._OPT_DIST_MIN_PPA               = -10.0e-20
        self._OPT_DIST_MAX_PPA               = 10.0e-20
        self._OPT_DIST_NBIN_PPA              = 1
        self._OPT_DIST_MIN_PPE               = 0.0
        self._OPT_DIST_MAX_PPE               = 10.0e-20
        self._OPT_DIST_NBIN_PPE              = 1
        self._OPT_DIST_MIN_PR                = 0.0
        self._OPT_DIST_MAX_PR                = 10.0e-20
        self._OPT_DIST_NBIN_PR               = 1
        self._OPT_DIST_MIN_PPHI              = 0.0
        self._OPT_DIST_MAX_PPHI              = 10.0e-20
        self._OPT_DIST_NBIN_PPHI             = 1
        self._OPT_DIST_MIN_PZ                = 0.0
        self._OPT_DIST_MAX_PZ                = 10.0e-20
        self._OPT_DIST_NBIN_PZ               = 1
        self._OPT_DIST_MIN_TIME              = 0.0
        self._OPT_DIST_MAX_TIME              = 10.0
        self._OPT_DIST_NBIN_TIME             = 1
        self._OPT_DIST_MIN_CHARGE            = -2
        self._OPT_DIST_MAX_CHARGE            = 3
        self._OPT_DIST_NBIN_CHARGE           = 1
        self._OPT_ENABLE_DIST_COM            = 0
        self._OPT_DIST_MIN_MU                = 0.0
        self._OPT_DIST_MAX_MU                = 2.5e-13
        self._OPT_DIST_NBIN_MU               = 100
        self._OPT_DIST_MIN_EKIN              = 0.0
        self._OPT_DIST_MAX_EKIN              = 1.0e-12
        self._OPT_DIST_NBIN_EKIN             = 50
        self._OPT_DIST_MIN_PTOR              = -1.0e-18
        self._OPT_DIST_MAX_PTOR              = 1.0e-18
        self._OPT_DIST_NBIN_PTOR             = 200
        self._OPT_ENABLE_ORBITWRITE          = 0
        self._OPT_ORBITWRITE_MODE            = 1
        self._OPT_ORBITWRITE_NPOINT          = 100
        self._OPT_ORBITWRITE_POLOIDALANGLES  = [0.0]
        self._OPT_ORBITWRITE_TOROIDALANGLES  = [0.0]
        self._OPT_ORBITWRITE_RADIALDISTANCES = [1.0]
        self._OPT_ORBITWRITE_INTERVAL        = 0.0
        self._OPT_ENABLE_TRANSCOEF           = 0
        self._OPT_TRANSCOEF_INTERVAL         = 0.0
        self._OPT_TRANSCOEF_NAVG             = 5
        self._OPT_TRANSCOEF_RECORDRHO        = 0

    @property
    def _SIM_MODE(self):
        """Simulation mode (1, 2, 3, 4)

        - 1 Gyro-orbit
        - 2 Guiding center
        - 3 Hybrid
        - 4 Magnetic field lines
        """
        return self._OPT_SIM_MODE

    @property
    def _ENABLE_ADAPTIVE(self):
        """Use adaptive time-step (0, 1)

        This option is used only if SIM_MODE = 2 or 3. Gyro-orbit
        simulations are always done with fixed time-step and magnetic field line
        simulations with adaptive time-step.

        - 0 Use fixed time-step
        - 1 Use adaptive time-step
        """
        return self._OPT_ENABLE_ADAPTIVE

    @property
    def _RECORD_MODE(self):
        """Record GOs as GCs in diagnostics (0, 1)

        - 0 Record GOs as GOs
        - 1 Record GOs as GCs
        """
        return self._OPT_RECORD_MODE

    @property
    def _FIXEDSTEP_USE_USERDEFINED(self):
        """Define fixed time-step value explicitly (0,1)

        Note: The adaptive scheme uses fixed time-step value as an initial step.

        - 0 Calculate time-step from FIXEDSTEP_NSTEPS_PER_GYROTIME
        - 1 Use opt.opt.FIXEDSTEP_USERDEFINED as a time-step
        """
        return self._OPT_FIXEDSTEP_USE_USERDEFINED

    @property
    def _FIXEDSTEP_USERDEFINED(self):
        """User-defined time-step [s]
        """
        return self._OPT_FIXEDSTEP_USERDEFINED

    @property
    def _FIXEDSTEP_GYRODEFINED(self):
        """Time-step is 2pi / ( gyrofrequency * N ) where N is this parameter
        """
        return self._OPT_FIXEDSTEP_GYRODEFINED

    @property
    def _ADAPTIVE_TOL_ORBIT(self):
        """Relative error tolerance for orbit following in adaptive scheme
        """
        return self._OPT_ADAPTIVE_TOL_ORBIT

    @property
    def _ADAPTIVE_TOL_CCOL(self):
        """Relative error tolerance for Coulomb collisions in adaptive scheme
        """
        return self._OPT_ADAPTIVE_TOL_CCOL

    @property
    def _ADAPTIVE_MAX_DRHO(self):
        """Maximum allowed change in rho during one time-step in adaptive scheme
        """
        return self._OPT_ADAPTIVE_MAX_DRHO

    @property
    def _ADAPTIVE_MAX_DPHI(self):
        """Maximum allowed change in phi during one time-step in adaptive scheme
        """
        return self._OPT_ADAPTIVE_MAX_DPHI

    @property
    def _ENDCOND_SIMTIMELIM(self):
        """Terminate when marker time passes ENDCOND_LIM_SIMTIME or when marker
        time has advanced ENDCOND_MAX_MILEAGE in a simulation

        Marker is terminated if t > ENDCOND_MAX_MILEAGE or
        t0 + t > ENDCOND_LIM_SIMTIME where t0 is marker's initial time and t
        the time it has been simulated. Note that if time is reversed, the
        simulation is instead terminated when t0 + t < ENDCOND_LIM_SIMTIME.
        See also ENDCOND_CPUTIMELIM.
        """
        return self._OPT_ENDCOND_SIMTIMELIM

    @property
    def _ENDCOND_CPUTIMELIM(self):
        """Terminate marker when the computer has spent ENDCOND_MAX_CPUTIME
        amount of real time to simulate it

        This limit should be rarely used  as its intended use is in debugging
        to stop markers stuck in a loop.
        """
        return self._OPT_ENDCOND_CPUTIMELIM

    @property
    def _ENDCOND_RHOLIM(self):
        """Terminate if marker goes outside given rho boundaries

        The boundaries are defined by ENDCOND_MIN_RHO and ENDCOND_MAX_RHO.
        """
        return self._OPT_ENDCOND_RHOLIM

    @property
    def _ENDCOND_ENERGYLIM(self):
        """Terminate when marker energy is below a user-specified value

        The user specified values are ENDCOND_MIN_ENERGY and
        ENDCOND_MIN_ENERGY_TIMES_THERMAL. Marker is terminated when either
        of these limits is reached.
        """
        return self._OPT_ENDCOND_ENERGYLIM

    @property
    def _ENDCOND_WALLHIT(self):
        """Terminate when marker intersects a wall element
        """
        return self._OPT_ENDCOND_WALLHIT

    @property
    def _ENDCOND_MAXORBS(self):
        """Terminate when marker has completed user-specified number of orbits

        Limit ENDCOND_MAX_TOROIDALORBS is used for a number of toroidal
        and ENDCOND_MAX_POLOIDALORBS for poloidal orbits.

        - 0 The end condition is not active.
        - 1 Marker is terminated when either of these limits is reached.
        - 2 Marker is terminated when both limits are reached.
        """
        return self._OPT_ENDCOND_MAXORBS

    @property
    def _ENDCOND_NEUTRALIZED(self):
        """Terminate when the marker becomes neutral
        """
        return self._OPT_ENDCOND_NEUTRALIZED

    @property
    def _ENDCOND_IONIZED(self):
        """Terminate when the marker becomes ionized
        """
        return self._OPT_ENDCOND_IONIZED

    @property
    def _ENDCOND_LIM_SIMTIME(self):
        """Time when the simulation stops [s]
        """
        return self._OPT_ENDCOND_LIM_SIMTIME

    @property
    def _ENDCOND_MAX_MILEAGE(self):
        """The maximum amount of time this marker is simulated [s]
        """
        return self._OPT_ENDCOND_MAX_MILEAGE

    @property
    def _ENDCOND_MAX_CPUTIME(self):
        """Maximum cpu time per marker [s]
        """
        return self._OPT_ENDCOND_MAX_CPUTIME

    @property
    def _ENDCOND_MAX_RHO(self):
        """Maximum rho value
        """
        return self._OPT_ENDCOND_MAX_RHO

    @property
    def _ENDCOND_MIN_RHO(self):
        """Minimum rho value
        """
        return self._OPT_ENDCOND_MIN_RHO

    @property
    def _ENDCOND_MIN_ENERGY(self):
        """Minimum energy [eV]
        """
        return self._OPT_ENDCOND_MIN_ENERGY

    @property
    def _ENDCOND_MIN_THERMAL(self):
        """Minimum energy limit is local ion thermal energy times this value
        """
        return self._OPT_ENDCOND_MIN_THERMAL

    @property
    def _ENDCOND_MAX_TOROIDALORBS(self):
        """Maximum number of toroidal orbits
        """
        return self._OPT_ENDCOND_MAX_TOROIDALORBS

    @property
    def _ENDCOND_MAX_POLOIDALORBS(self):
        """Maximum number of poloidal orbits
        """
        return self._OPT_ENDCOND_MAX_POLOIDALORBS

    @property
    def _ENABLE_ORBIT_FOLLOWING(self):
        """Trace markers in an electromagnetic field
        """
        return self._OPT_ENABLE_ORBIT_FOLLOWING
    
    @property
    def _ENABLE_RF(self):
        """Include RF effects in orbit-following
        """
        return self._OPT_ENABLE_RF

    @property
    def _ENABLE_COULOMB_COLLISIONS(self):
        """Markers experience Coulomb collisions with background plasma
        """
        return self._OPT_ENABLE_COULOMB_COLLISIONS

    @property
    def _ENABLE_MHD(self):
        """Include MHD perturbations to orbit-following
        """
        return self._OPT_ENABLE_MHD

    @property
    def _ENABLE_ATOMIC(self):
        """Markers can undergo atomic reactions with background plasma
        or neutrals

        - 0 Atomic reactions are turned off.
        - 1 Atomic reactions are on but marker is terminated when outside
            the reaction data domain.
        - 2 Atomic reactions are on but they are ignored when marker is
            outside the reaction data domain.
        """
        return self._OPT_ENABLE_ATOMIC

    @property
    def _DISABLE_FIRSTORDER_GCTRANS(self):
        """Disable first order guiding center transformation in velocity space
        """
        return self._OPT_DISABLE_FIRSTORDER_GCTRANS

    @property
    def _DISABLE_ENERGY_CCOLL(self):
        """Disable guiding center energy collisions
        """
        return self._OPT_DISABLE_ENERGY_CCOLL

    @property
    def _DISABLE_PITCH_CCOLL(self):
        """Disable guiding center pitch collisions
        """
        return self._OPT_DISABLE_PITCH_CCOLL

    @property
    def _DISABLE_GCDIFF_CCOLL(self):
        """Disable guiding center spatial diffusion
        """
        return self._OPT_DISABLE_GCDIFF_CCOLL

    @property
    def _REVERSE_TIME(self):
         """Trace markers backwards in time.

         Collision operator isn't reversible so disable collisions if this
         option is used. Also when tracing markers, the simulation stops when
         marker time is below ENDCOND_LIM_SIMTIME.
         """
         return self._OPT_REVERSE_TIME


    @property
    def _ENABLE_DIST_5D(self):
        """Collect distribution histogram in [R, phi, z, ppa, ppe, t, q]

        The coordinates are:

        - R   major radius
        - phi toroidal angle
        - z   z-coordinate
        - ppa momentum component parallel to magnetic field
        - ppe momentum component perpendicular to magnetic field
        - t   time
        - q   charge
        """
        return self._OPT_ENABLE_DIST_5D

    @property
    def _ENABLE_DIST_6D(self):
        """Collect distribution histogram in [R, phi, z, pR, pphi, pz, t, q]

        The coordinates are:

        - R    major radius
        - phi  toroidal angle
        - z    z-coordinate
        - pR   momentum R-component
        - pphi momentum phi-component
        - pz   momentum z-component
        - t    time
        - q    charge
        """
        return self._OPT_ENABLE_DIST_6D

    @property
    def _ENABLE_DIST_RHO5D(self):
        """Collect distribution histogram in [rho, pol, phi, ppa, ppe, t, q]

        The coordinates are:

        - rho  flux surface
        - pol  poloidal angle
        - phi  toroidal angle
        - z    z-coordinate
        - ppa  momentum component parallel to magnetic field
        - ppe  momentum component perpendicular to magnetic field
        - t    time
        - q    charge
        """
        return self._OPT_ENABLE_DIST_RHO5D

    @property
    def _ENABLE_DIST_RHO6D(self):
        """Collect distribution histogram in [rho, pol, phi, pR, pphi, pz, t, q]

        The coordinates are:

        - rho  flux surface
        - pol  poloidal angle
        - phi  toroidal angle
        - z    z-coordinate
        - pR   momentum R-component
        - pphi momentum phi-component
        - pz   momentum z-component
        - t    time
        - q    charge
        """
        return self._OPT_ENABLE_DIST_RHO6D

    @property
    def _DIST_MIN_R(self):
        """Minimum bin edge for R coordinate [m]
        """
        return self._OPT_DIST_MIN_R

    @property
    def _DIST_MAX_R(self):
        """Maximum bin edge for R coordinate [m]
        """
        return self._OPT_DIST_MAX_R

    @property
    def _DIST_NBIN_R(self):
        """Number of bins the interval [MIN_R, MAX_R] is divided to
        """
        return self._OPT_DIST_NBIN_R

    @property
    def _DIST_MIN_PHI(self):
        """Minimum bin edge for phi coordinate [deg]
        """
        return self._OPT_DIST_MIN_PHI

    @property
    def _DIST_MAX_PHI(self):
        """Maximum bin edge for phi coordinate [deg]
        """
        return self._OPT_DIST_MAX_PHI

    @property
    def _DIST_NBIN_PHI(self):
        """Number of bins the interval [MIN_PHI, MAX_PHI] is divided to
        """
        return self._OPT_DIST_NBIN_PHI

    @property
    def _DIST_MIN_Z(self):
        """Minimum bin edge for z coordinate [m]
        """
        return self._OPT_DIST_MIN_Z

    @property
    def _DIST_MAX_Z(self):
        """Maximum bin edge for z coordinate [m]
        """
        return self._OPT_DIST_MAX_Z

    @property
    def _DIST_NBIN_Z(self):
        """Number of bins the interval [MIN_Z, MAX_Z] is divided to
        """
        return self._OPT_DIST_NBIN_Z

    @property
    def _DIST_MIN_RHO(self):
        """Minimum bin edge for rho coordinate [1]
        """
        return self._OPT_DIST_MIN_RHO

    @property
    def _DIST_MAX_RHO(self):
        """Maximum bin edge for rho coordinate [1]
        """
        return self._OPT_DIST_MAX_RHO

    @property
    def _DIST_NBIN_RHO(self):
        """Number of bins the interval [MIN_RHO, MAX_RHO] is divided to
        """
        return self._OPT_DIST_NBIN_RHO

    @property
    def _DIST_MIN_THETA(self):
        """Minimum bin edge for theta coordinate [deg]
        """
        return self._OPT_DIST_MIN_THETA

    @property
    def _DIST_MAX_THETA(self):
        """Maximum bin edge for theta coordinate [deg]
        """
        return self._OPT_DIST_MAX_THETA

    @property
    def _DIST_NBIN_THETA(self):
        """Number of bins the interval [MIN_THETA, MAX_THETA] is divided to
        """
        return self._OPT_DIST_NBIN_THETA

    @property
    def _DIST_MIN_PPA(self):
        """Minimum bin edge for ppa coordinate [kg m/s]
        """
        return self._OPT_DIST_MIN_PPA

    @property
    def _DIST_MAX_PPA(self):
        """Maximum bin edge for ppa coordinate [kg m/s]
        """
        return self._OPT_DIST_MAX_PPA

    @property
    def _DIST_NBIN_PPA(self):
        """Number of bins the interval [MIN_PPA, MAX_PPA] is divided to
        """
        return self._OPT_DIST_NBIN_PPA

    @property
    def _DIST_MIN_PPE(self):
        """Minimum bin edge for ppe coordinate [kg m/s]
        """
        return self._OPT_DIST_MIN_PPE

    @property
    def _DIST_MAX_PPE(self):
        """Maximum bin edge for ppe coordinate [kg m/s]
        """
        return self._OPT_DIST_MAX_PPE

    @property
    def _DIST_NBIN_PPE(self):
        """Number of bins the interval [MIN_PPE, MAX_PPE] is divided to
        """
        return self._OPT_DIST_NBIN_PPE

    @property
    def _DIST_MIN_PR(self):
        """Minimum bin edge for pR coordinate [kg m/s]
        """
        return self._OPT_DIST_MIN_PR

    @property
    def _DIST_MAX_PR(self):
        """Maximum bin edge for pR coordinate [kg m/s]
        """
        return self._OPT_DIST_MAX_PR

    @property
    def _DIST_NBIN_PR(self):
        """Number of bins the interval [MIN_PR, MAX_PR] is divided to
        """
        return self._OPT_DIST_NBIN_PR

    @property
    def _DIST_MIN_PPHI(self):
        """Minimum bin edge for pphi coordinate [kg m/s]
        """
        return self._OPT_DIST_MIN_PPHI

    @property
    def _DIST_MAX_PPHI(self):
        """Maximum bin edge for pphi coordinate [kg m/s]
        """
        return self._OPT_DIST_MAX_PPHI

    @property
    def _DIST_NBIN_PPHI(self):
        """Number of bins the interval [MIN_PPHI, MAX_PPHI] is divided to
        """
        return self._OPT_DIST_NBIN_PPHI

    @property
    def _DIST_MIN_PZ(self):
        """Minimum bin edge for pz coordinate [kg m/s]
        """
        return self._OPT_DIST_MIN_PZ

    @property
    def _DIST_MAX_PZ(self):
        """Maximum bin edge for pz coordinate [kg m/s]
        """
        return self._OPT_DIST_MAX_PZ

    @property
    def _DIST_NBIN_PZ(self):
        """Number of bins the interval [MIN_PZ, MAX_PZ] is divided to
        """
        return self._OPT_DIST_NBIN_PZ

    @property
    def _DIST_MIN_TIME(self):
        """Minimum bin edge for time coordinate [s]
        """
        return self._OPT_DIST_MIN_TIME

    @property
    def _DIST_MAX_TIME(self):
        """Maximum bin edge for time coordinate [s]
        """
        return self._OPT_DIST_MAX_TIME

    @property
    def _DIST_NBIN_TIME(self):
        """Number of bins the interval [MIN_TIME, MAX_TIME] is divided to
        """
        return self._OPT_DIST_NBIN_TIME

    @property
    def _DIST_MIN_CHARGE(self):
        """Minimum bin edge for charge coordinate [e]
        """
        return self._OPT_DIST_MIN_CHARGE

    @property
    def _DIST_MAX_CHARGE(self):
        """Maximum bin edge for charge coordinate [e]
        """
        return self._OPT_DIST_MAX_CHARGE

    @property
    def _DIST_NBIN_CHARGE(self):
        """Number of bins the interval [MIN_CHARGE, MAX_CHARGE] is divided to
        """
        return self._OPT_DIST_NBIN_CHARGE

    @property
    def _ENABLE_DIST_COM(self):
        """Collect constant-of-motion distribution histogram [mu, Ekin, Ptor]

        The coordinates are:

        - mu   magnetic moment (first adiabatic invariant)
        - Ekin kinetic energy
        - Ptor canonical toroidal angular momentum
        """
        return self._OPT_ENABLE_DIST_COM

    @property
    def _DIST_MIN_MU(self):
        """Minimum bin edge for the mu coordinate in COM-dist [J/T]
        """
        return self._OPT_DIST_MIN_MU

    @property
    def _DIST_MAX_MU(self):
        """Maximum bin edge for the mu coordinate in COM-dist [J/T]
        """
        return self._OPT_DIST_MAX_MU

    @property
    def _DIST_NBIN_MU(self):
        """Number of bins the interval [MIN_MU, MAX_MU] is divided to
        """
        return self._OPT_DIST_NBIN_MU

    @property
    def _DIST_MIN_EKIN(self):
        """Minimum bin edge for the energy coordinate [J]
        """
        return self._OPT_DIST_MIN_EKIN

    @property
    def _DIST_MAX_EKIN(self):
        """Maximum bin edge for the energy coordinate [J]
        """
        return self._OPT_DIST_MAX_EKIN

    @property
    def _DIST_NBIN_EKIN(self):
        """Number of bins the interval [MIN_EKIN, MAX_EKIN] is divided to
        """
        return self._OPT_DIST_NBIN_EKIN

    @property
    def _DIST_MIN_PTOR(self):
        """Minimum bin edge for the Ptor coordinate [kg m^2/s]
        """
        return self._OPT_DIST_MIN_PTOR

    @property
    def _DIST_MAX_PTOR(self):
        """Maximum bin edge for the Ptor coordinate [kg m^2/s]
        """
        return self._OPT_DIST_MAX_PTOR

    @property
    def _DIST_NBIN_PTOR(self):
        """Number of bins the interval [MIN_PTOR, MAX_PTOR] is divided to
        """
        return self._OPT_DIST_NBIN_PTOR

    @property
    def _ENABLE_ORBITWRITE(self):
        """Enable diagnostics that store marker orbit

        - 0 Marker orbit diagnostics are not collected
        - 1 Marker orbit diagnostics are collected
        """
        return self._OPT_ENABLE_ORBITWRITE

    @property
    def _ORBITWRITE_MODE(self):
        """What kind of marker orbit diagnostics are collected

        These are only used if ENABLE_ORBITWRITE is active.

        - 0 When marker crosses a plane (Poincare-plot)
        - 1 Between given time intervals
        """
        return self._OPT_ORBITWRITE_MODE

    @property
    def _ORBITWRITE_NPOINT(self):
        """Maximum number of points (per marker) to be written

        If this number is exceeded when marker is being simulated, the oldest
        points will be replaced as long as the simulation continues. Thus,
        this parameter is effectively the number of marker's last positions
        that are stored.
        """
        return self._OPT_ORBITWRITE_NPOINT

    @property
    def _ORBITWRITE_POLOIDALANGLES(self):
        """Poloidal angles of toroidal planes where toroidal plots are collected

        Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 0.
        """
        return self._OPT_ORBITWRITE_POLOIDALANGLES

    @property
    def _ORBITWRITE_TOROIDALANGLES(self):
        """Toroidal angles of poloidal planes where poloidal plots are collected

        Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 0.
        """
        return self._OPT_ORBITWRITE_TOROIDALANGLES

    @property
    def _ORBITWRITE_RADIALDISTANCES(self):
        """Minor radius coordinate where radial plots are collected

        Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 0.
        """
        return self._OPT_ORBITWRITE_RADIALDISTANCES

    @property
    def _ORBITWRITE_INTERVAL(self):
        """Time interval for writing marker state [s]

        Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 1.
        """
        return self._OPT_ORBITWRITE_INTERVAL

    @property
    def _ENABLE_TRANSCOEF(self):
        """Enable evaluation of transport coefficients.

        - 0 Transport coefficients are not collected
        - 1 Transport coefficients are collected
        """
        return self._OPT_ENABLE_TRANSCOEF

    @property
    def _TRANSCOEF_INTERVAL(self):
        """Time interval for recording data points

        The data points are recorded (at the outer mid-plane crossing) if this
        interval has passed from the previous recording.
        """
        return self._OPT_TRANSCOEF_INTERVAL

    @property
    def _TRANSCOEF_NAVG(self):
        """Number of subsequent data points that are averaged before calculating
        coefficients to reduce noise.
        """
        return self._OPT_TRANSCOEF_NAVG

    @property
    def _TRANSCOEF_RECORDRHO(self):
        """Record coefficients in terms of normalized poloidal flux instead of
        meters.
        """
        return self._OPT_TRANSCOEF_RECORDRHO

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.

        Raises
        ------
        AscotIOException
            If stored options were unviable.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        defopt = Opt.get_default()
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                if key not in defopt.keys():
                    warnings.warn("Unknown option " + key + " ignored.")
                    continue
                val = f[path][key][:]

                # Take type from the default parameter
                if isinstance(defopt[key], list):
                    try:
                        val[0]
                        out[key] = type(defopt[key])(val)
                    except Exception:
                        out[key] = val
                else:
                    out[key] = type(defopt[key])(val)

        for o in defopt.keys():
            if o not in out:
                raise ValueError("Missing parameter: " + o)

        return out

    def tostring(self, descriptions=True, aslist=False):
        """Convert options to string representation.

        Parameters
        ----------
        descriptions : bool, optional
            If True, section headers and descriptions are added above
            the parameters to make options more readable.
        aslist : bool, optional
            Instead of line breaks, separate each line as a list item.

        Returns
        -------
        opt : str or list [str]
            String where option parameters are written on each line as
            "<PARAMETER> = <VALUE>".
        """
        opt = self.read()

        def trim(val):
            """Make sure values are not too accurate or have too many leading
            zeroes.
            """
            if np.abs(val) >= 1e4 or (np.abs(val) <=1e-4 and val != 0):
                return "{0:e}".format(val)
            else:
                return "{0:g}".format(val)

        # Convert all parameters to strings first
        for o, val in opt.items():
            if isinstance(val, list):
                for i, v in enumerate(val):
                    val[i] = float(trim(v))
                opt[o] = str(val)
            else:
                opt[o] = trim(val)

        # Create an array of strings
        out = []
        for o in self.__dict__.keys():
            if len(o) < 5 or o[:5] != "_OPT_":
                continue

            # Extract parameter info
            attr = o[4:]
            desc = getattr(Opt, attr).__doc__
            name = getattr(Opt, attr).fget.__name__[1:]
            dval = opt[name]

            # Add decorations
            if descriptions:

                def makebanner(title):
                    """Makes a banner with given title.
                    """
                    nchar = len(title)
                    ljust = int((76 - nchar) / 2)
                    banner = []
                    banner += ["#"  + "*"*78 +  "#"]
                    banner += ["#*" + " "*ljust + title
                               + " "*(76 - nchar - ljust) + "*#"]
                    banner += ["#*" + " "*76 + "*#"]
                    banner += ["#"  + "*"*78 +  "#"]
                    banner += [""]
                    return banner

                if name == "SIM_MODE":
                    out.extend(makebanner("SIMULATION MODE AND TIME-STEP"))
                elif name == "ENDCOND_SIMTIMELIM":
                    out.extend(makebanner("END CONDITIONS"))
                elif name == "ENABLE_ORBIT_FOLLOWING":
                    out.extend(makebanner("PHYSICS"))
                elif name == "ENABLE_DIST_5D":
                    out.extend(makebanner("DISTRIBUTIONS"))
                elif name == "ENABLE_ORBITWRITE":
                    out.extend(makebanner("ORBIT WRITE"))
                elif name == "ENABLE_TRANSCOEF":
                    out.extend(makebanner("TRANSPORT COEFFICIENT"))

                # Clean docstrings a little bit by removing extra whispace and
                # empty lines.
                for d in desc.splitlines():
                    d = d.lstrip(" ")
                    if len(d) > 1: out.append("# " + d)

            # Write parameter and value after possible decorations
            out.append(name + " = " + str(dval))
            out.append("") # Empty line to separate parameters

        if not aslist:
            string = ""
            for o in out:
                string += o + "\n"

            out = string

        return out

    def new(self, desc=None, **kwargs):
        """Write new options with updated parameters.

        This method reads the current options, updates the given parameters,
        and writes the updated options as a new input.

        Parameters
        ----------
        desc : str, optional
            Input description.
        **kwargs
            <name> : <value> pairs for each options parameter that are updated.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If arguments contain unknown parameters.
        """
        options = self.read()
        for o, val in kwargs.items():
            if not o in options:
                raise ValueError("Unknown parameter: " + o)
            options[o] = val

        return self._root.create_input("opt", desc=desc, **options)

    @staticmethod
    def write_hdf5(fn, desc=None, **kwargs):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to HDF5 file.
        desc : str, optional
            Input description.
        **kwargs
            <name> : <value> pairs for each options parameter

            Hint: :meth:`get_default` generates a default options dictionary
            from which one can pick parameters to change.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If arguments contain unknown parameters or if there are parameters
            missing.
        """
        opt = Opt.get_default()
        for o in opt.keys():
            if not o in kwargs:
                raise ValueError("Missing parameter: " + o)

        for o in kwargs.keys():
            if not o in opt:
                raise ValueError("Unknown parameter: " + o)

        parent = "options"
        group  = "opt"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            # Convert all options to numpy float arrays and write
            for param, data in kwargs.items():
                data = np.asarray(data)

                data = data.astype("f8")
                d = g.create_dataset(param, (data.size,), data=data)

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The default options are created.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return Opt.get_default()

    @staticmethod
    def get_default():
        """Get default option parameters.

        The default options have all "ENABLE" settings disabled. The default
        values are not suitable for any proper simulation.

        Returns
        -------
        options : dict
            Default options in a format that can be passed to
            :meth:`write_hdf5`.
        """
        out = {}
        defopt = Opt(None, None)
        for opt in defopt.__dict__.keys():
            if len(opt) < 5 or opt[:5] != "_OPT_":
                continue

            attr = opt[4:]
            name = getattr(Opt, attr).fget.__name__[1:]
            dval = getattr(defopt, attr)
            out[name] = dval

        return out

    @staticmethod
    def validate():
        """
        Validate options and return names of faulty parameters.
        """
        opt  = Opt.get_default()
        opt0 = {}
        for o in opt:
            if len(o) > 2:
                opt0[o[0]] = o[3]

        msg = []
        opt = self.read()
        for o in opt0.keys():
            if o in opt.keys() and not opt0[o].validate(opt[o]):
                msg += [o]

        return msg

    @staticmethod
    def convert_string(lines):
        """Converts a string to :class:`Opt` input.

        Parameters
        ----------
        lines : str or [str]
            String with line breaks or list of strings where each item is
            a single line.

            Options are interpreted from a string assuming that the relevant
            lines have format "<PARAMETER> = <VALUE>" and nothing else.
            The parameter name must be in all uppercase. Lines that don't start
            with a recognizable parameter are ignored.

        Returns
        -------
        out : dict
            String converted as an input for :meth:`write_hdf5`.
        """
        if not isinstance(lines, list):
            lines = lines.splitlines()

        opt = {}
        defopt = Opt.get_default()
        for line in lines:
            # Remove empty lines and lines with comments
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue

            # Separate key-value pairs
            paramval = line.split("=")
            param = paramval[0].strip()
            val   = paramval[1].strip()

            # Make sure that a parameter is a PARAMETER
            if param != param.upper():
                continue
            if param not in defopt.keys():
                warnings.warn("Possibly unknown option " + param + " ignored.")
                continue
            if isinstance(defopt[param], list):
                # Value is (supposed to be) a list
                try:
                    opt[param] = [float(val)]
                except ValueError:
                   opt[param] = ast.literal_eval(val)
            else:
                # Value is a plain value
                opt[param] = type(defopt[param])(val)

        return opt

    @staticmethod
    def convert_xml(xml):
        """Convert XML file to :class:`Opt` input.

        Parameters
        ----------
        xml : str
            Path to XML file or contents of the file as a string.

        Returns
        -------
        out : dict
            XML converted to :class:`Opt` input.
        """
        # Is input file a string
        tree = None
        try:
            tree = ET.ElementTree().parse(xml)
        except:
            pass
        try:
            if tree is None: tree = ET.fromstring(xml)
        except:
            pass
        if tree is None:
            raise ValueError("Could not parse XML file or string.")

        opt = Opt.get_default()
        for tag in opt.keys():
            val = tree.findall("*"+tag)
            if len(val) == 0:
                warnings.warn(
                    f"Using default value for {tag} as it was not found in xml")
                continue
            opt[tag] = val[0].text
        return opt

    @staticmethod
    def schema(opt=None, fnxsd=None, fnxml=None):
        """Produce a schema and an XML file containing ASCOT5 option parameters.

        Parameters
        ----------
        opt : dict, optional
            Options used to fill the XML with parameter values.

            If None, the default options are used.
        fnxsd : str, optional
            None or filename to write XSD data.
        fnxml : str, optional
            None or filename to write XML data.

        Returns
        -------
        schema :
            Schema that contains parameter descriptions and which can validate
            XML files containing ASCOT5 options.
        xml :
            XML file containing the options parameters and their values.
        """
        if xmlschema is None:
            raise ImportError(
                "This feature requires package xmlschema:\n"
                "pip install xmlschema")

        def doc(name, stype):
            """Helper function to add parameters with docstrings, types, and
            proper indentation.

            Paramaters
            ----------
            name : str
                Name of the parameter as it appears in option attributes
                (without preceding underscore).
            stype : str
                This parameters simple type.

            Returns
            -------
            fstr
                Formatted string that specifies XML element.
            """
            attr = Opt.__dict__["_" + name]
            docstring = attr.__doc__.replace('<', '&lt;').replace('>', '&gt;')
            text = f"""
                <xs:element name="{name}" type="{stype}">
                    <xs:annotation>
                    <xs:documentation>
        {textwrap.indent(docstring, "            ")}
                    </xs:documentation>
                    </xs:annotation>
                </xs:element>"""
            return text

        xsd = textwrap.dedent(f"""\
            <?xml version="1.0"?>
            <xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">

                <xs:simpleType name="IntegerPositive">
                    <xs:restriction base="xs:integer">
                    <xs:minInclusive value="1"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="IntegerBinary">
                    <xs:restriction base="xs:integer">
                    <xs:minInclusive value="0"/>
                    <xs:maxInclusive value="1"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="Integer012">
                    <xs:restriction base="xs:integer">
                    <xs:minInclusive value="0"/>
                    <xs:maxInclusive value="2"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="Integer1234">
                    <xs:restriction base="xs:integer">
                    <xs:minInclusive value="1"/>
                    <xs:maxInclusive value="4"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="FloatPositive">
                    <xs:restriction base="xs:float">
                    <xs:minExclusive value="0.0"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="FloatNonNegative">
                    <xs:restriction base="xs:float">
                    <xs:minInclusive value="0.0"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="Float0to360">
                    <xs:restriction base="xs:float">
                    <xs:minInclusive value="0.0"/>
                    <xs:maxInclusive value="360.0"/>
                    </xs:restriction>
                </xs:simpleType>

                <xs:simpleType name="FloatNNList">
                    <xs:list itemType="FloatNonNegative"/>
                </xs:simpleType>

                <xs:simpleType name="FloatOrFloatNNList">
                    <xs:union memberTypes="xs:float FloatNNList"/>
                </xs:simpleType>

                <xs:element name="parameters">
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="SIMULATION_MODE_AND_TIMESTEP"/>
                        <xs:element ref="END_CONDITIONS"/>
                        <xs:element ref="PHYSICS"/>
                        <xs:element ref="DISTRIBUTIONS"/>
                        <xs:element ref="ORBIT_WRITE"/>
                        <xs:element ref="TRANSPORT_COEFFICIENT"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>

                <xs:element name="SIMULATION_MODE_AND_TIMESTEP">
                    <xs:annotation>
                    <xs:documentation>
                    Simulation mode and time-step
                    </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="SIM_MODE"/>
                        <xs:element ref="ENABLE_ADAPTIVE"/>
                        <xs:element ref="FIXEDSTEP_USE_USERDEFINED"/>
                        <xs:element ref="FIXEDSTEP_USERDEFINED"/>
                        <xs:element ref="FIXEDSTEP_GYRODEFINED"/>
                        <xs:element ref="ADAPTIVE_TOL_ORBIT"/>
                        <xs:element ref="ADAPTIVE_TOL_CCOL"/>
                        <xs:element ref="ADAPTIVE_MAX_DRHO"/>
                        <xs:element ref="ADAPTIVE_MAX_DPHI"/>
                        <xs:element ref="RECORD_MODE"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>
            {doc('SIM_MODE',                  'Integer1234')}
            {doc('ENABLE_ADAPTIVE',           'IntegerBinary')}
            {doc('FIXEDSTEP_USE_USERDEFINED', 'IntegerBinary')}
            {doc('FIXEDSTEP_USERDEFINED',     'FloatPositive')}
            {doc('FIXEDSTEP_GYRODEFINED',     'IntegerPositive')}
            {doc('ADAPTIVE_TOL_ORBIT',        'FloatPositive')}
            {doc('ADAPTIVE_TOL_CCOL',         'FloatPositive')}
            {doc('ADAPTIVE_MAX_DRHO',         'FloatPositive')}
            {doc('ADAPTIVE_MAX_DPHI',         'FloatPositive')}
            {doc('RECORD_MODE',               'IntegerBinary')}

                <xs:element name="END_CONDITIONS">
                    <xs:annotation>
                    <xs:documentation>
                    End conditions
                    </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="ENDCOND_SIMTIMELIM"/>
                        <xs:element ref="ENDCOND_CPUTIMELIM"/>
                        <xs:element ref="ENDCOND_RHOLIM"/>
                        <xs:element ref="ENDCOND_ENERGYLIM"/>
                        <xs:element ref="ENDCOND_WALLHIT"/>
                        <xs:element ref="ENDCOND_MAXORBS"/>
                        <xs:element ref="ENDCOND_NEUTRALIZED"/>
                        <xs:element ref="ENDCOND_IONIZED"/>
                        <xs:element ref="ENDCOND_LIM_SIMTIME"/>
                        <xs:element ref="ENDCOND_MAX_MILEAGE"/>
                        <xs:element ref="ENDCOND_MAX_CPUTIME"/>
                        <xs:element ref="ENDCOND_MAX_RHO"/>
                        <xs:element ref="ENDCOND_MIN_RHO"/>
                        <xs:element ref="ENDCOND_MIN_ENERGY"/>
                        <xs:element ref="ENDCOND_MIN_THERMAL"/>
                        <xs:element ref="ENDCOND_MAX_TOROIDALORBS"/>
                        <xs:element ref="ENDCOND_MAX_POLOIDALORBS"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>
            {doc('ENDCOND_SIMTIMELIM',       'IntegerBinary')}
            {doc('ENDCOND_CPUTIMELIM',       'IntegerBinary')}
            {doc('ENDCOND_RHOLIM',           'IntegerBinary')}
            {doc('ENDCOND_ENERGYLIM',        'IntegerBinary')}
            {doc('ENDCOND_WALLHIT',          'IntegerBinary')}
            {doc('ENDCOND_MAXORBS',          'IntegerBinary')}
            {doc('ENDCOND_NEUTRALIZED',      'IntegerBinary')}
            {doc('ENDCOND_IONIZED',          'IntegerBinary')}
            {doc('ENDCOND_LIM_SIMTIME',      'FloatNonNegative')}
            {doc('ENDCOND_MAX_MILEAGE',      'FloatNonNegative')}
            {doc('ENDCOND_MAX_CPUTIME',      'FloatNonNegative')}
            {doc('ENDCOND_MAX_RHO',          'FloatNonNegative')}
            {doc('ENDCOND_MIN_RHO',          'FloatNonNegative')}
            {doc('ENDCOND_MIN_ENERGY',       'FloatPositive')}
            {doc('ENDCOND_MIN_THERMAL',      'FloatPositive')}
            {doc('ENDCOND_MAX_TOROIDALORBS', 'IntegerPositive')}
            {doc('ENDCOND_MAX_POLOIDALORBS', 'IntegerPositive')}

                <xs:element name="PHYSICS">
                    <xs:annotation>
                    <xs:documentation>
                    Physics
                    </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="ENABLE_ORBIT_FOLLOWING"/>
                        <xs:element ref="ENABLE_COULOMB_COLLISIONS"/>
                        <xs:element ref="ENABLE_MHD"/>
                        <xs:element ref="ENABLE_ATOMIC"/>
                        <xs:element ref="DISABLE_FIRSTORDER_GCTRANS"/>
                        <xs:element ref="DISABLE_ENERGY_CCOLL"/>
                        <xs:element ref="DISABLE_PITCH_CCOLL"/>
                        <xs:element ref="DISABLE_GCDIFF_CCOLL"/>
                        <xs:element ref="REVERSE_TIME"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>
            {doc('ENABLE_ORBIT_FOLLOWING',     'IntegerBinary')}
            {doc('ENABLE_COULOMB_COLLISIONS',  'IntegerBinary')}
            {doc('ENABLE_MHD',                 'IntegerBinary')}
            {doc('ENABLE_ATOMIC',              'Integer012')}
            {doc('DISABLE_FIRSTORDER_GCTRANS', 'IntegerBinary')}
            {doc('DISABLE_ENERGY_CCOLL',       'IntegerBinary')}
            {doc('DISABLE_PITCH_CCOLL',        'IntegerBinary')}
            {doc('DISABLE_GCDIFF_CCOLL',       'IntegerBinary')}
            {doc('REVERSE_TIME',               'IntegerBinary')}

                <xs:element name="DISTRIBUTIONS">
                    <xs:annotation>
                    <xs:documentation>
                    Distribution output
                    </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="ENABLE_DIST_5D"/>
                        <xs:element ref="ENABLE_DIST_RHO5D"/>
                        <xs:element ref="ENABLE_DIST_6D"/>
                        <xs:element ref="ENABLE_DIST_RHO6D"/>
                        <xs:element ref="DIST_MIN_R"/>
                        <xs:element ref="DIST_MAX_R"/>
                        <xs:element ref="DIST_NBIN_R"/>
                        <xs:element ref="DIST_MIN_PHI"/>
                        <xs:element ref="DIST_MAX_PHI"/>
                        <xs:element ref="DIST_NBIN_PHI"/>
                        <xs:element ref="DIST_MIN_Z"/>
                        <xs:element ref="DIST_MAX_Z"/>
                        <xs:element ref="DIST_NBIN_Z"/>
                        <xs:element ref="DIST_MIN_RHO"/>
                        <xs:element ref="DIST_MAX_RHO"/>
                        <xs:element ref="DIST_NBIN_RHO"/>
                        <xs:element ref="DIST_MIN_THETA"/>
                        <xs:element ref="DIST_MAX_THETA"/>
                        <xs:element ref="DIST_NBIN_THETA"/>
                        <xs:element ref="DIST_MIN_PPA"/>
                        <xs:element ref="DIST_MAX_PPA"/>
                        <xs:element ref="DIST_NBIN_PPA"/>
                        <xs:element ref="DIST_MIN_PPE"/>
                        <xs:element ref="DIST_MAX_PPE"/>
                        <xs:element ref="DIST_NBIN_PPE"/>
                        <xs:element ref="DIST_MIN_PR"/>
                        <xs:element ref="DIST_MAX_PR"/>
                        <xs:element ref="DIST_NBIN_PR"/>
                        <xs:element ref="DIST_MIN_PPHI"/>
                        <xs:element ref="DIST_MAX_PPHI"/>
                        <xs:element ref="DIST_NBIN_PPHI"/>
                        <xs:element ref="DIST_MIN_PZ"/>
                        <xs:element ref="DIST_MAX_PZ"/>
                        <xs:element ref="DIST_NBIN_PZ"/>
                        <xs:element ref="DIST_MIN_TIME"/>
                        <xs:element ref="DIST_MAX_TIME"/>
                        <xs:element ref="DIST_NBIN_TIME"/>
                        <xs:element ref="DIST_MIN_CHARGE"/>
                        <xs:element ref="DIST_MAX_CHARGE"/>
                        <xs:element ref="DIST_NBIN_CHARGE"/>
                        <xs:element ref="ENABLE_DIST_COM"/>
                        <xs:element ref="DIST_MIN_MU"/>
                        <xs:element ref="DIST_MAX_MU"/>
                        <xs:element ref="DIST_NBIN_MU"/>
                        <xs:element ref="DIST_MIN_EKIN"/>
                        <xs:element ref="DIST_MAX_EKIN"/>
                        <xs:element ref="DIST_NBIN_EKIN"/>
                        <xs:element ref="DIST_MIN_PTOR"/>
                        <xs:element ref="DIST_MAX_PTOR"/>
                        <xs:element ref="DIST_NBIN_PTOR"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>
            {doc('ENABLE_DIST_5D',    'IntegerBinary')}
            {doc('ENABLE_DIST_RHO5D', 'IntegerBinary')}
            {doc('ENABLE_DIST_6D',    'IntegerBinary')}
            {doc('ENABLE_DIST_RHO6D', 'IntegerBinary')}
            {doc('DIST_MIN_R',        'FloatPositive')}
            {doc('DIST_MAX_R',        'FloatPositive')}
            {doc('DIST_NBIN_R',       'IntegerPositive')}
            {doc('DIST_MIN_PHI',      'Float0to360')}
            {doc('DIST_MAX_PHI',      'Float0to360')}
            {doc('DIST_NBIN_PHI',     'IntegerPositive')}
            {doc('DIST_MIN_Z',        'xs:float')}
            {doc('DIST_MAX_Z',        'xs:float')}
            {doc('DIST_NBIN_Z',       'IntegerPositive')}
            {doc('DIST_MIN_RHO',      'FloatNonNegative')}
            {doc('DIST_MAX_RHO',      'FloatNonNegative')}
            {doc('DIST_NBIN_RHO',     'IntegerPositive')}
            {doc('DIST_MIN_THETA',    'Float0to360')}
            {doc('DIST_MAX_THETA',    'Float0to360')}
            {doc('DIST_NBIN_THETA',   'IntegerPositive')}
            {doc('DIST_MIN_PPA',      'xs:float')}
            {doc('DIST_MAX_PPA',      'xs:float')}
            {doc('DIST_NBIN_PPA',     'IntegerPositive')}
            {doc('DIST_MIN_PPE',      'FloatNonNegative')}
            {doc('DIST_MAX_PPE',      'FloatNonNegative')}
            {doc('DIST_NBIN_PPE',     'IntegerPositive')}
            {doc('DIST_MIN_PR',       'xs:float')}
            {doc('DIST_MAX_PR',       'xs:float')}
            {doc('DIST_NBIN_PR',      'IntegerPositive')}
            {doc('DIST_MIN_PPHI',     'xs:float')}
            {doc('DIST_MAX_PPHI',     'xs:float')}
            {doc('DIST_NBIN_PPHI',    'IntegerPositive')}
            {doc('DIST_MIN_PZ',       'xs:float')}
            {doc('DIST_MAX_PZ',       'xs:float')}
            {doc('DIST_NBIN_PZ',      'IntegerPositive')}
            {doc('DIST_MIN_TIME',     'xs:float')}
            {doc('DIST_MAX_TIME',     'xs:float')}
            {doc('DIST_NBIN_TIME',    'IntegerPositive')}
            {doc('DIST_MIN_CHARGE',   'xs:integer')}
            {doc('DIST_MAX_CHARGE',   'xs:integer')}
            {doc('DIST_NBIN_CHARGE',  'IntegerPositive')}
            {doc('ENABLE_DIST_COM',   'IntegerBinary')}
            {doc('DIST_MIN_MU',       'FloatNonNegative')}
            {doc('DIST_MAX_MU',       'FloatNonNegative')}
            {doc('DIST_NBIN_MU',      'IntegerPositive')}
            {doc('DIST_MIN_EKIN',     'FloatNonNegative')}
            {doc('DIST_MAX_EKIN',     'FloatNonNegative')}
            {doc('DIST_NBIN_EKIN',    'IntegerPositive')}
            {doc('DIST_MIN_PTOR',     'xs:float')}
            {doc('DIST_MAX_PTOR',     'xs:float')}
            {doc('DIST_NBIN_PTOR',    'IntegerPositive')}

                <xs:element name="ORBIT_WRITE">
                    <xs:annotation>
                    <xs:documentation>
                    Orbit output
                    </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="ENABLE_ORBITWRITE"/>
                        <xs:element ref="ORBITWRITE_MODE"/>
                        <xs:element ref="ORBITWRITE_NPOINT"/>
                        <xs:element ref="ORBITWRITE_POLOIDALANGLES"/>
                        <xs:element ref="ORBITWRITE_TOROIDALANGLES"/>
                        <xs:element ref="ORBITWRITE_RADIALDISTANCES"/>
                        <xs:element ref="ORBITWRITE_INTERVAL"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>
            {doc('ENABLE_ORBITWRITE',          'IntegerBinary')}
            {doc('ORBITWRITE_MODE',            'IntegerBinary')}
            {doc('ORBITWRITE_NPOINT',          'IntegerPositive')}
            {doc('ORBITWRITE_POLOIDALANGLES',  'FloatOrFloatNNList')}
            {doc('ORBITWRITE_TOROIDALANGLES',  'FloatOrFloatNNList')}
            {doc('ORBITWRITE_RADIALDISTANCES', 'FloatOrFloatNNList')}
            {doc('ORBITWRITE_INTERVAL',        'FloatNonNegative')}

            <xs:element name="TRANSPORT_COEFFICIENT">
                    <xs:annotation>
                    <xs:documentation>
                    Transport coefficient output
                    </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                    <xs:all>
                        <xs:element ref="ENABLE_TRANSCOEF"/>
                        <xs:element ref="TRANSCOEF_INTERVAL"/>
                        <xs:element ref="TRANSCOEF_NAVG"/>
                        <xs:element ref="TRANSCOEF_RECORDRHO"/>
                    </xs:all>
                    </xs:complexType>
                </xs:element>
            {doc('ENABLE_TRANSCOEF',    'IntegerBinary')}
            {doc('TRANSCOEF_INTERVAL',  'FloatNonNegative')}
            {doc('TRANSCOEF_NAVG',      'IntegerPositive')}
            {doc('TRANSCOEF_RECORDRHO', 'IntegerBinary')}
            </xs:schema>""")
        schema = xmlschema.XMLSchema(xsd)
        if opt is None: opt = Opt.get_default()
        opt_hierarchy = {
            "SIMULATION_MODE_AND_TIMESTEP":{}, "END_CONDITIONS":{},
            "PHYSICS":{}, "DISTRIBUTIONS":{}, "ORBIT_WRITE":{},
            "TRANSPORT_COEFFICIENT":{}}
        for k, v in opt.items():
            if k == "SIM_MODE":
                grp = "SIMULATION_MODE_AND_TIMESTEP"
            elif k == "ENDCOND_SIMTIMELIM":
                grp = "END_CONDITIONS"
            elif k == "ENABLE_ORBIT_FOLLOWING":
                grp = "PHYSICS"
            elif k == "ENABLE_DIST_5D":
                grp = "DISTRIBUTIONS"
            elif k == "ENABLE_ORBITWRITE":
                grp = "ORBIT_WRITE"
            elif k == "ENABLE_TRANSCOEF":
                grp = "TRANSPORT_COEFFICIENT"
            opt_hierarchy[grp][k] = v

        data = json.dumps({"parameters":opt_hierarchy})
        xml = xmlschema.from_json(data, schema=schema, preserve_root=True)

        if fnxsd is not None:
            with open(fnxsd, 'w') as f:
                f.writelines(xsd)
        if fnxml is not None:
            ET.ElementTree(xml).write(fnxml)

        return schema, xml
