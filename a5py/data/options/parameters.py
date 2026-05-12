"""Defines option parameters and their documentation.

The option parameters are organized into dataclasses.
Each dataclass is named in lowercase to match the attribute name under which it
is exposed in the ``Options`` object.

This convention improves readability in the generated documentation.
For example, the documentation will show ``Options.simulation`` followed by
``simulation.simulation_mode`` rather than ``Simulation.simulation_mode``.

The parameters themselves are defined as properties since this allows us to
generate the documentation more easily while also making the parameters
immutable.
"""
from dataclasses import dataclass

import numpy as np
from .constraints import (
    summarize_property, parse_constraint, enforce_constraint
)

@dataclass
class simulation():

    _simulation_mode: int = 1
    _record_mode: int = 0
    _timestep: float = 1.0e-8
    _enable_adaptive: bool = True
    _adaptive_tolerance_orbit: float = 1.0e-8
    _adaptive_tolerance_collisions: float = 1.0e-1

    @property
    def simulation_mode(self):
        """Simulation mode: {1, 2, 3, 4}, default=1.

        - 1: Gyro-orbit
        - 2: Guiding center
        - 3: Hybrid
        - 4: Magnetic field lines
        """
        return self._simulation_mode

    @simulation_mode.setter
    def simulation_mode(self, value):
        self._simulation_mode = value

    @property
    def record_mode(self):
        """Change the physical picture before collecting diagnostics: {0, 1},
        default=0.

        This option only affects the gyro-orbit simulations.

        - 0: Record gyro-orbits as they are.
        - 1: Instead, record the guiding center position of the gyro-orbit.
        """
        return self._record_mode

    @record_mode.setter
    def record_mode(self, value):
        self._record_mode = value

    @property
    def timestep(self):
        """User-defined time-step [s]: (> 0), default=1.0e-8.

        This time-step is used as the fixed time step and as an initial step for
        adaptive time-stepping.
        """
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        self._timestep = value

    @property
    def enable_adaptive(self):
        """Use adaptive time-step: {0, 1}, default=1.

        This option is used only if ``simulation_mode`` is 2 or 3. Gyro-orbit
        simulations are always done with fixed time-step and magnetic field line
        simulations with adaptive time-step. Note: The adaptive scheme uses
        fixed time-step value as an initial step.

        - 0: Use fixed time-step.
        - 1: Use adaptive time-step.
        """
        return self._enable_adaptive

    @enable_adaptive.setter
    def enable_adaptive(self, value):
        self._enable_adaptive = value

    @property
    def adaptive_tolerance_orbit(self):
        """Relative error tolerance for orbit following in adaptive scheme:
        (> 0), default=1.0e-8.
        """
        return self._adaptive_tolerance_orbit

    @adaptive_tolerance_orbit.setter
    def adaptive_tolerance_orbit(self, value):
        self._adaptive_tolerance_orbit = value

    @property
    def adaptive_tolerance_collisions(self):
        """Relative error tolerance for Coulomb collisions in adaptive scheme:
        (> 0), default=1.0e-1.
        """
        return self._adaptive_tolerance_collisions

    @adaptive_tolerance_collisions.setter
    def adaptive_tolerance_collisions(self, value):
        self._adaptive_tolerance_collisions = value


@dataclass
class physics():

    _enable_orbit_following: bool = 0
    _enable_coulomb_collisions: bool = 0
    _enable_mhd: bool = 0
    _enable_atomic: bool = 0
    _enable_icrh: bool = 0
    _enable_aldforce: bool = 0
    _disable_first_order_gctransformation: bool = 0
    _disable_ccoll_gcenergy: bool = 0
    _disable_ccoll_gcpitch: bool = 0
    _disable_ccoll_gcspatial: bool = 0
    _reverse_time: bool = 0

    @property
    def enable_orbit_following(self):
        """Trace markers in an electromagnetic field: {0, 1}, default=0."""
        return self._enable_orbit_following

    @enable_orbit_following.setter
    def enable_orbit_following(self, value):
        self._enable_orbit_following = value

    @property
    def enable_coulomb_collisions(self):
        """Markers experience Coulomb collisions with background plasma:
        {0, 1}, default=0.
        """
        return self._enable_coulomb_collisions

    @enable_coulomb_collisions.setter
    def enable_coulomb_collisions(self, value):
        self._enable_coulomb_collisions = value

    @property
    def enable_mhd(self):
        """Include MHD perturbations to orbit-following: {0, 1}, default=0."""
        return self._enable_mhd

    @enable_mhd.setter
    def enable_mhd(self, value):
        self._enable_mhd = value

    @property
    def enable_atomic(self):
        """Markers can undergo atomic reactions with background plasma
        or neutrals: {0, 1, 2}, default=0.

        - 0: Atomic reactions are turned off.
        - 1: Atomic reactions are on but marker is terminated when outside
             the reaction data domain.
        - 2: Atomic reactions are on but they are ignored when marker is
             outside the reaction data domain.
        """
        return self._enable_atomic

    @enable_atomic.setter
    def enable_atomic(self, value):
        self._enable_atomic = value

    @property
    def enable_icrh(self):
        """Enable ion cyclotron resonance heating operator: {0, 1}, default=0.

        ICRH operator transfers energy (via "kicks") to ions when they are on
        the resonance. The code must be compiled with RFOF=1 and the RFOF
        library must be present in order to use the ICRH operator.
        """
        return self._enable_icrh

    @enable_icrh.setter
    def enable_icrh(self, value):
        self._enable_icrh = value

    @property
    def enable_aldforce(self):
        """Enable radiation reaction force (synchrotron losses): {0, 1},
        default=0.

        The radiation reaction force (a.k.a. Abraham-Lorentz-Dirac or ALD force)
        causes charged particles to lose energy via radiation. The losses are
        proportional to the particle energy and inversely proportional to the
        particle mass, making this option mostly relevant for (runaway)
        electrons.
        """
        return self._enable_aldforce

    @enable_aldforce.setter
    def enable_aldforce(self, value):
        self._enable_aldforce = value

    @property
    def disable_first_order_gctransformation(self):
        """Disable first order guiding center transformation in velocity space:
        {0, 1}, default=0.
        """
        return self._disable_first_order_gctransformation

    @disable_first_order_gctransformation.setter
    def disable_first_order_gctransformation(self, value):
        self._disable_first_order_gctransformation = value

    @property
    def disable_ccoll_gcenergy(self):
        """Disable guiding center energy collisions: {0, 1}, default=0."""
        return self._disable_ccoll_gcenergy

    @disable_ccoll_gcenergy.setter
    def disable_ccoll_gcenergy(self, value):
        self._disable_ccoll_gcenergy = value

    @property
    def disable_ccoll_gcpitch(self):
        """Disable guiding center pitch collisions: {0, 1}, default=0."""
        return self._disable_ccoll_gcpitch

    @disable_ccoll_gcpitch.setter
    def disable_ccoll_gcpitch(self, value):
        self._disable_ccoll_gcpitch = value

    @property
    def disable_ccoll_gcspatial(self):
        """Disable guiding center spatial diffusion: {0, 1}, default=0."""
        return self._disable_ccoll_gcspatial

    @disable_ccoll_gcspatial.setter
    def disable_ccoll_gcspatial(self, value):
        self._disable_ccoll_gcspatial = value

    @property
    def reverse_time(self):
         """Trace markers backwards in time: {0, 1}, default=0.

         Collision operator isn't reversible so disable collisions if this
         option is used. Also when tracing markers, the simulation stops when
         marker time is *below* ``simulation_time_limit``.
         """
         return self._reverse_time

    @reverse_time.setter
    def reverse_time(self, value):
        self._reverse_time = value

@dataclass
class endconditions():

    _activate_simulation_time_limits: bool = 0
    _activate_real_time_limit: bool = 0
    _activate_rho_limit: bool = 0
    _activate_energy_limits: bool = 0
    _activate_wall_hits: bool = 0
    _activate_orbit_limit: bool = 0
    _activate_neutralization: bool = 0
    _activate_ionization: bool = 0
    _lab_time_limit: float = 1.0
    _max_mileage: float = 1.0
    _max_real_time: float = 3600.0
    _rho_coordinate_limits: tuple[float, float] = (0.0, 2.0)
    _min_energy: float = 1.0e3
    _min_local_thermal_energy: float = 2.0
    _max_number_of_toroidal_orbits: float = 100
    _max_number_of_poloidal_orbits: float = 100

    @property
    def activate_simulation_time_limits(self):
        """Terminate marker based on laboratory time or its lifetime: {0, 1},
        default=0.

        Terminate when either of the following is true:

        1. Absolute time limit: The marker's current time (in laboratory time)
           exceeds ``lab_time_limit``.

        2. Relative lifetime limit: The marker's elapsed time since its birth
           exceeds ``max_mileage``.
        """
        return self._activate_simulation_time_limits

    @activate_simulation_time_limits.setter
    def activate_simulation_time_limits(self, value):
        self._activate_simulation_time_limits = value

    @property
    def activate_real_time_limit(self):
        """Terminate marker when the computer has spent specified
        amount of real time to simulate it: {0, 1}, default=0.

        This is not a "proper" end condition in a sense that it does not
        correspond to any physical process. This should be used just to control
        simulation duration or debugging. The limit is set by ``max_real_time``.
        """
        return self._activate_real_time_limit

    @activate_real_time_limit.setter
    def activate_real_time_limit(self, value):
        self._activate_real_time_limit = value

    @property
    def activate_rho_limit(self):
        """Terminate if marker goes outside given rho boundaries: {0, 1},
        default=0.

        The boundaries are defined by ``rho_coordinate_limits``.
        """
        return self._activate_rho_limit

    @activate_rho_limit.setter
    def activate_rho_limit(self, value):
        self._activate_rho_limit = value

    @property
    def activate_energy_limits(self):
        """Terminate when marker energy is below a user-specified value: {0, 1},
        default=0.

        The user specified values are ``min_energy`` and
        ``min_local_thermal_energy``. Marker is terminated when either
        of these limits is reached.
        """
        return self._activate_energy_limits

    @activate_energy_limits.setter
    def activate_energy_limits(self, value):
        self._activate_energy_limits = value

    @property
    def activate_wall_hits(self):
        """Terminate when marker intersects a wall element: {0, 1}, default=0.
        """
        return self._activate_wall_hits

    @activate_wall_hits.setter
    def activate_wall_hits(self, value):
        self._activate_wall_hits = value

    @property
    def activate_orbit_limit(self):
        """Terminate when marker has completed user-specified number of orbits:
        {0, 1, 2}, default=0.

        The number of toroidal and poloidal orbits is limited by
        ``max_number_of_toroidal_orbits`` and ``max_number_of_poloidal_orbits``,
        respectively.

        - 0: The end condition is not active.
        - 1: Marker is terminated when either of these limits is reached.
        - 2: Marker is terminated when both limits are reached.
        """
        return self._activate_orbit_limit

    @activate_orbit_limit.setter
    def activate_orbit_limit(self, value):
        self._activate_orbit_limit = value

    @property
    def activate_neutralization(self):
        """Terminate when the marker becomes neutral: {0, 1}, default=0."""
        return self._activate_neutralization

    @activate_neutralization.setter
    def activate_neutralization(self, value):
        self._activate_neutralization = value

    @property
    def activate_ionization(self):
        """Terminate when the marker becomes ionized: {0, 1}, default=0."""
        return self._activate_ionization

    @activate_ionization.setter
    def activate_ionization(self, value):
        self._activate_ionization = value

    @property
    def lab_time_limit(self):
        """Laboratory time when the simulation stops [s]: (> 0), default=1.0."""
        return self._lab_time_limit

    @lab_time_limit.setter
    def lab_time_limit(self, value):
        self._lab_time_limit = value

    @property
    def max_mileage(self):
        """The maximum amount of time this marker is simulated [s] or [m]:
        (> 0), default=1.0.
        """
        return self._max_mileage

    @max_mileage.setter
    def max_mileage(self, value):
        self._max_mileage = value

    @property
    def max_real_time(self):
        """Maximum real time spent simulating a marker [s]: (> 0),
        default=3600.0
        """
        return self._max_real_time

    @max_real_time.setter
    def max_real_time(self, value):
        self._max_real_time = value

    @property
    def rho_coordinate_limits(self):
        """Minimum and maximum values for rho: [a > 0, b > 0],
        default=[0.0, 1.0].
        """
        return self._rho_coordinate_limits

    @rho_coordinate_limits.setter
    def rho_coordinate_limits(self, value):
        self._rho_coordinate_limits = value

    @property
    def min_energy(self):
        """Minimum energy [eV]: (> 0), default=1000.0."""
        return self._min_energy

    @min_energy.setter
    def min_energy(self, value):
        self._min_energy = value

    @property
    def min_local_thermal_energy(self):
        """Minimum energy limit is local ion thermal energy times this value:
        (> 0), default=2.0.
        """
        return self._min_local_thermal_energy

    @min_local_thermal_energy.setter
    def min_local_thermal_energy(self, value):
        self._min_local_thermal_energy = value

    @property
    def max_number_of_toroidal_orbits(self):
        """Maximum number of toroidal orbits: (> 0), default=1000."""
        return self._max_number_of_toroidal_orbits

    @max_number_of_toroidal_orbits.setter
    def max_number_of_toroidal_orbits(self, value):
        self._max_number_of_toroidal_orbits = value

    @property
    def max_number_of_poloidal_orbits(self):
        """Maximum number of poloidal orbits: (> 0), default=1000."""
        return self._max_number_of_poloidal_orbits

    @max_number_of_poloidal_orbits.setter
    def max_number_of_poloidal_orbits(self, value):
        self._max_number_of_poloidal_orbits = value

@dataclass
class histograms():

    _abscissae: tuple[str, ...] = ()
    _bins: tuple[int, ...] = ()
    _interval: tuple[tuple[float, float], ...] = ()
    _charge_interval: tuple[int, int] = ()

    @property
    def abscissae(self):
        r"""Set of distribution abscissae coordinates.

        The available coordinates are:

        - :math:`r`: major radius
        - :math:`phi`: toroidal angle
        - :math:`z`: axial height
        - :math:`p_\parallel`: momentum component parallel to magnetic field
        - :math:`p_\perp`: momentum component perpendicular to magnetic field
        - :math:`t`: time
        """
        return self._abscissae

    @abscissae.setter
    def abscissae(self, value):
        self._abscissae = value

    @property
    def interval(self):
        r"""Abscissa limits for abscissa coordinates: [a > 0, b > 0],
        default=[0.1, 10.0].
        """
        return self._interval

    @interval.setter
    def interval(self, value):
        self._interval = value

    @property
    def bins(self):
        """Number of bins the corresponding ``interval`` is divided to: (> 0),
        default=10.
        """
        return self._bins

    @bins.setter
    def bins(self, value):
        self._bins = value

    @property
    def charge_interval(self):
        """Abscissa limits (inclusive) for test particle charge [e]: [a, b],
        default=[2, 2].

        For each charge state on the interval, the distribution is collected
        separately.
        """
        return self._charge_interval

    @charge_interval.setter
    def charge_interval(self, value):
        self._charge_interval = value

@dataclass
class orbit():

    _poincare: bool = 0
    _interval: float = 0.0
    _collect_orbit: bool = 0
    _poloidal_angles: list[float] | tuple[float] | float = (0.0,)
    _toroidal_angles: list[float] | tuple[float] | float = (0.0,)
    _radial_distances: list[float] | tuple[float] | float = (1.0,)
    _number_of_points_per_marker: int = 100

    @property
    def collect_orbit(self):
        """Enable diagnostics that store marker orbit: {0, 1}, default=0.

        - 0 Marker orbit diagnostics are not collected
        - 1 Marker orbit diagnostics are collected
        """
        return self._collect_orbit

    @collect_orbit.setter
    def collect_orbit(self, value):
        self._collect_orbit = value

    @property
    def poincare(self):
        """Collect data only when a (Poincaré) plane is crossed: {0, 1},
        default=0.

        Enable this option in order to generate Poincaré plots.
        """
        return self._poincare

    @poincare.setter
    def poincare(self, value):
        self._poincare = value

    @property
    def number_of_points_per_marker(self):
        """Maximum number of points (per marker) to be written: (> 0),
        default=100.

        If this number is exceeded when marker is being simulated, the oldest
        points will be replaced as long as the simulation continues. Thus,
        this parameter is effectively the number of marker's last positions
        that are stored.
        """
        return self._number_of_points_per_marker

    @number_of_points_per_marker.setter
    def number_of_points_per_marker(self, value):
        self._number_of_points_per_marker = value

    @property
    def poloidal_angles(self):
        """Poloidal angles of toroidal planes where toroidal Poincaré plots are
        collected [deg]: [0 <= a <= 360, ...], default=[0.0, 180.0].

        Used when ``poincare`` is enabled.
        """
        return np.asarray(self._poloidal_angles)

    @poloidal_angles.setter
    def poloidal_angles(self, value):
        self._poloidal_angles = value

    @property
    def toroidal_angles(self):
        """Toroidal angles of poloidal planes where poloidal Poincaré plots are
        collected [deg]: [0 <= a <= 360, ...], default=[0.0, 180.0].

        Used when ``poincare`` is enabled.
        """
        return np.asarray(self._toroidal_angles)

    @toroidal_angles.setter
    def toroidal_angles(self, value):
        self._toroidal_angles = value

    @property
    def radial_distances(self):
        """Minor radius coordinate where radial Poincaré plots are collected:
        [0 < a <= 1, ...], default=[1.0].

        Used when ``poincare`` is enabled.
        """
        return np.asarray(self._radial_distances).copy()

    @radial_distances.setter
    def radial_distances(self, value):
        self._radial_distances = value

    @property
    def interval(self):
        """Time interval for writing marker state [s]: (> 0), default=0.

        Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 1.
        """
        return self._interval

    @interval.setter
    def interval(self, value):
        self._interval = value
