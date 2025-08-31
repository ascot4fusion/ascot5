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

@dataclass
class simulation():

    _simulation_mode: int = 1
    _record_mode: int = 0
    _use_explicit_fixedstep: bool = 0
    _explicit_fixedstep: float = 1.0e-8
    _gyrodefined_fixedstep: int = 20
    _enable_adaptive: bool = True
    _adaptive_tolerance_orbit: float = 1.0e-8
    _adaptive_tolerance_collisions: float = 1.0e-1
    _adaptive_max_drho: float = 0.1
    _adaptive_max_dphi: float = 2.0

    @property
    def simulation_mode(self):
        """Simulation mode: {1, 2, 3, 4}, default=1.

        - 1: Gyro-orbit
        - 2: Guiding center
        - 3: Hybrid
        - 4: Magnetic field lines
        """
        return self._simulation_mode

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

    @property
    def record_mode(self):
        """Change the physical picture before collecting diagnostics: {0, 1},
        default=0.

        - 0: Record gyro-orbits as they are.
        - 1: Instead, record the guiding center position of the gyro-orbit.
        """
        return self._record_mode

    @property
    def use_explicit_fixedstep(self):
        """Define the fixed time-step value explicitly: {0, 1}, default=0.

        - 0: The time-step is defined by ``gyrodefined_fixedstep``.
        - 1: The time-step is defined by ``explicit_fixedstep``.
        """
        return self._use_explicit_fixedstep

    @property
    def explicit_fixedstep(self):
        """User-defined time-step [s]: (> 0), default=1.0e-8.
        """
        return self._explicit_fixedstep

    @property
    def gyrodefined_fixedstep(self):
        r"""How many gyro-periods are in one time-step: (> 0), default=20.

        In other words, the time-step is defined by
        :math:`2\pi / N\omega_\mathrm{gyro}` where :math:`N` is this parameter.
        """
        return self._gyrodefined_fixedstep

    @property
    def adaptive_tolerance_orbit(self):
        """Relative error tolerance for orbit following in adaptive scheme:
        (> 0), default=1.0e-8.
        """
        return self._adaptive_tolerance_orbit

    @property
    def adaptive_tolerance_collisions(self):
        """Relative error tolerance for Coulomb collisions in adaptive scheme:
        (> 0), default=1.0e-1.
        """
        return self._adaptive_tolerance_collisions

    @property
    def adaptive_max_drho(self):
        """Maximum allowed change in rho during one time-step in adaptive
        scheme: (> 0), default=0.1.
        """
        return self._adaptive_max_drho

    @property
    def adaptive_max_dphi(self):
        """Maximum allowed change in phi during one time-step in adaptive
        scheme: (> 0), default=2.0.
        """
        return self._adaptive_max_dphi


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

    @property
    def enable_coulomb_collisions(self):
        """Markers experience Coulomb collisions with background plasma:
        {0, 1}, default=0.
        """
        return self._enable_coulomb_collisions

    @property
    def enable_mhd(self):
        """Include MHD perturbations to orbit-following: {0, 1}, default=0."""
        return self._enable_mhd

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

    @property
    def enable_icrh(self):
        """Enable ion cyclotron resonance heating operator: {0, 1}, default=0.

        ICRH operator transfers energy (via "kicks") to ions when they are on
        the resonance. The code must be compiled with RFOF=1 and the RFOF
        library must be present in order to use the ICRH operator.
        """
        return self._enable_icrh

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

    @property
    def disable_first_order_gctransformation(self):
        """Disable first order guiding center transformation in velocity space:
        {0, 1}, default=0.
        """
        return self._disable_first_order_gctransformation

    @property
    def disable_ccoll_gcenergy(self):
        """Disable guiding center energy collisions: {0, 1}, default=0."""
        return self._disable_ccoll_gcenergy

    @property
    def disable_ccoll_gcpitch(self):
        """Disable guiding center pitch collisions: {0, 1}, default=0."""
        return self._disable_ccoll_gcpitch

    @property
    def disable_ccoll_gcspatial(self):
        """Disable guiding center spatial diffusion: {0, 1}, default=0."""
        return self._disable_ccoll_gcspatial

    @property
    def reverse_time(self):
         """Trace markers backwards in time: {0, 1}, default=0.

         Collision operator isn't reversible so disable collisions if this
         option is used. Also when tracing markers, the simulation stops when
         marker time is *below* ``simulation_time_limit``.
         """
         return self._reverse_time


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

    @property
    def activate_real_time_limit(self):
        """Terminate marker when the computer has spent specified
        amount of real time to simulate it: {0, 1}, default=0.

        This is not a "proper" end condition in a sense that it does not
        correspond to any physical process. This should be used just to control
        simulation duration or debugging. The limit is set by ``max_real_time``.
        """
        return self._activate_real_time_limit

    @property
    def activate_rho_limit(self):
        """Terminate if marker goes outside given rho boundaries: {0, 1},
        default=0.

        The boundaries are defined by ``rho_coordinate_limits``.
        """
        return self._activate_rho_limit

    @property
    def activate_energy_limits(self):
        """Terminate when marker energy is below a user-specified value: {0, 1},
        default=0.

        The user specified values are ``min_energy`` and
        ``min_local_thermal_energy``. Marker is terminated when either
        of these limits is reached.
        """
        return self._activate_energy_limits

    @property
    def activate_wall_hits(self):
        """Terminate when marker intersects a wall element: {0, 1}, default=0.
        """
        return self._activate_wall_hits

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

    @property
    def activate_neutralization(self):
        """Terminate when the marker becomes neutral: {0, 1}, default=0."""
        return self._activate_neutralization

    @property
    def activate_ionization(self):
        """Terminate when the marker becomes ionized: {0, 1}, default=0."""
        return self._activate_ionization

    @property
    def lab_time_limit(self):
        """Laboratory time when the simulation stops [s]: (> 0), default=1.0."""
        return self._lab_time_limit

    @property
    def max_mileage(self):
        """The maximum amount of time this marker is simulated [s] or [m]:
        (> 0), default=1.0.
        """
        return self._max_mileage

    @property
    def max_real_time(self):
        """Maximum real time spent simulating a marker [s]: (> 0),
        default=3600.0
        """
        return self._max_real_time

    @property
    def rho_coordinate_limits(self):
        """Minimum and maximum values for rho: [a > 0, b > 0],
        default=[0.0, 1.0].
        """
        return self._rho_coordinate_limits

    @property
    def min_energy(self):
        """Minimum energy [eV]: (> 0), default=1000.0."""
        return self._min_energy

    @property
    def min_local_thermal_energy(self):
        """Minimum energy limit is local ion thermal energy times this value:
        (> 0), default=2.0.
        """
        return self._min_local_thermal_energy

    @property
    def max_number_of_toroidal_orbits(self):
        """Maximum number of toroidal orbits: (> 0), default=1000."""
        return self._max_number_of_toroidal_orbits

    @property
    def max_number_of_poloidal_orbits(self):
        """Maximum number of poloidal orbits: (> 0), default=1000."""
        return self._max_number_of_poloidal_orbits


@dataclass
class distributions():

    _collect_dist5d: bool = 0
    _collect_dist6d: bool = 0
    _collect_dist5drho: bool = 0
    _collect_dist6drho: bool = 0

    _r_bins: int = 10
    _z_bins: int = 20
    _phi_bins: int = 1
    _rho_bins: int = 20
    _theta_bins: int = 1
    _r_interval: tuple[float, float] = 0.1, 10.0
    _z_interval: tuple[float, float] = -6.0, 6.0
    _phi_interval: tuple[float, float] = 0.0, 360.0
    _rho_interval: tuple[float, float] = 0.0, 1.0
    _theta_interval: tuple[float, float] = 0.0, 360.0

    _pr_bins: int = 20
    _pz_bins: int = 20
    _pphi_bins: int = 20
    _ppara_bins: int = 40
    _pperp_bins: int = 20
    _pr_interval: tuple[float, float] = 0.0, 10.0e-20
    _pz_interval: tuple[float, float] = 0.0, 10.0e-20
    _pphi_interval: tuple[float, float] = 0.0, 10.0e-20
    _ppara_interval: tuple[float, float] = -10.0e-20, 10.0e-20
    _pperp_interval: tuple[float, float] = 0.0, 10.0e-20

    _time_bins: int = 1
    _time_interval: tuple[float, float] = 0.0, 10.0
    _charge_interval: tuple[int, int] = 2, 2

    @property
    def collect_dist5d(self):
        r"""Collect distribution histogram in [R, phi, z, ppa, ppe, t]:
        {0, 1}, default=0.

        The coordinates are:

        - :math:`R`: major radius
        - :math:`phi`: toroidal angle
        - :math:`z`: axial height
        - :math:`p_\parallel`: momentum component parallel to magnetic field
        - :math:`p_\perp`: momentum component perpendicular to magnetic field
        - :math:`t`: time
        """
        return self._collect_dist5d

    @property
    def collect_dist6d(self):
        r"""Collect distribution histogram in [R, phi, z, pR, pphi, pz, t]:
        {0, 1}, default=0.

        The coordinates are:

        - :math:`R`: major radius
        - :math:`phi`: toroidal angle
        - :math:`z`: axial height
        - :math:`p_R` momentum R-component
        - :math:`p_phi` momentum phi-component
        - :math:`p_z` momentum z-component
        - :math:`t`: time
        """
        return self._collect_dist6d

    @property
    def collect_dist5drho(self):
        r"""Collect distribution histogram in [rho, theta, phi, ppa, ppe, t]:
        {0, 1}, default=0.

        The coordinates are:

        - :math:`\rho` flux surface
        - :math:`\theta` poloidal angle
        - :math:`phi` toroidal angle
        - :math:`p_\parallel`: momentum component parallel to magnetic field
        - :math:`p_\perp`: momentum component perpendicular to magnetic field
        - :math:`t`: time
        """
        return self._collect_dist5drho

    @property
    def collect_dist6drho(self):
        r"""Collect distribution histogram in
        [rho, theta, phi, pR, pphi, pz, t]: {0, 1}, default=0.

        The coordinates are:

        - :math:`\rho` flux surface
        - :math:`\theta` poloidal angle
        - :math:`phi` toroidal angle
        - :math:`p_R` momentum R-component
        - :math:`p_phi` momentum phi-component
        - :math:`p_z` momentum z-component
        - :math:`t`: time
        """
        return self._collect_dist6drho

    @property
    def r_interval(self):
        r"""Abscissa limits for :math:`R` coordinate [m]: [a > 0, b > 0],
        default=[0.1, 10.0].
        """
        return self._r_interval

    @property
    def r_bins(self):
        """Number of bins the ``r_interval`` is divided to: (> 0), default=10.
        """
        return self._r_bins

    @property
    def phi_interval(self):
        r"""Abscissa limits for :math:`phi` coordinate [deg]:
        [0 <= a < 360, 0 < b <= 360], default=[0.0, 360.0].
        """
        return self._phi_interval

    @property
    def phi_bins(self):
        """Number of bins the ``phi_interval`` is divided to: (> 0), default=1.
        """
        return self._phi_bins

    @property
    def z_interval(self):
        r"""Abscissa limits for :math:`z` coordinate [m]: [a, b],
        default=[-10.0, 10.0].
        """
        return self._z_interval

    @property
    def z_bins(self):
        """Number of bins the ``z_interval`` is divided to: (> 0), default=20.
        """
        return self._z_bins

    @property
    def rho_interval(self):
        r"""Abscissa limits for :math:`rho` coordinate: [a >= 0, b > 0],
        default=[0.0, 1.0].
        """
        return self._rho_interval

    @property
    def rho_bins(self):
        """Number of bins the ``rho_interval`` is divided to: (> 0), default=20.
        """
        return self._rho_bins

    @property
    def theta_interval(self):
        r"""Abscissa limits for :math:`\theta` coordinate [deg]:
        [0 <= a < 360, 0 < b <= 360], default=[0.0, 360.0].
        """
        return self._theta_interval

    @property
    def theta_bins(self):
        """Number of bins the ``theta_interval`` is divided to: (> 0),
        default=1.
        """
        return self._theta_bins

    @property
    def ppara_interval(self):
        r"""Abscissa limits for :math:`p_\parallel` coordinate [eV]:
        [a, b], default=[-4.5e6, 4.5e6].
        """
        return self._ppara_interval

    @property
    def ppara_bins(self):
        """Number of bins the ``ppara_interval`` is divided to: (> 0),
        default=1.
        """
        return self._ppara_bins

    @property
    def pperp_interval(self):
        r"""Abscissa limits for :math:`p_\perp` coordinate [eV]:
        [a >= 0, b > 0], default=[0.0, 4.5e6].
        """
        return self._pperp_interval

    @property
    def pperp_bins(self):
        """Number of bins the ``pperp_interval`` is divided to: (> 0),
        default=1.
        """
        return self._pperp_bins

    @property
    def pr_interval(self):
        r"""Abscissa limits for :math:`p_R` coordinate [eV]: [a >= 0, b > 0],
        default=[0.0, 4.5e6].
        """
        return self._pr_interval

    @property
    def pr_bins(self):
        """Number of bins the ``pr_interval`` is divided to: (> 0), default=1.
        """
        return self._pr_bins

    @property
    def pphi_interval(self):
        r"""Abscissa limits for :math:`p_phi` coordinate [eV]: [a >= 0, b > 0],
        default=[0.0, 4.5e6].
        """
        return self._pphi_interval

    @property
    def pphi_bins(self):
        """Number of bins the ``pphi_interval`` is divided to
        """
        return self._pphi_bins

    @property
    def pz_interval(self):
        r"""Abscissa limits for :math:`p_z` coordinate [eV]: [a >= 0, b > 0],
        default=[0.0, 4.5e6].
        """
        return self._pz_interval

    @property
    def pz_bins(self):
        """Number of bins the ``pz_interval`` is divided to: (> 0), default=1.
        """
        return self._pz_bins

    @property
    def time_interval(self):
        r"""Abscissa limits for :math:`t` coordinate [s]: [a, b],
        default=[0.0, 1.0].
        """
        return self._time_interval

    @property
    def time_bins(self):
        """Number of bins the ``time_interval`` is divided to: (> 0),
        default=1.
        """
        return self._time_bins

    @property
    def charge_interval(self):
        """Abscissa limits (inclusive) for test particle charge [e]: [a, b],
        default=[2, 2].

        For each charge state on the interval, the distribution is collected
        separately.
        """
        return self._charge_interval


@dataclass
class comdistribution():

    _mu_bins: int = 100
    _ekin_bins: int = 50
    _ptor_bins: int = 200
    _mu_interval: tuple[float, float] = 0.0, 2.5e-13
    _ekin_interval: tuple[float, float] = 0.0, 1.0e-12
    _ptor_interval: tuple[float, float] = -1.0e-18, 1.0e-18
    _collect_distcom: bool = 0

    @property
    def collect_distcom(self):
        r"""Collect constant-of-motion distribution histogram [mu, Ekin, Ptor]:
        {0, 1}, default=0.

        The coordinates are:

        - :math:`\mu` magnetic moment (first adiabatic invariant)
        - :math:`E_{kin}` kinetic energy
        - :math:`P_{tor}` canonical toroidal angular momentum
        """
        return self._collect_distcom

    @property
    def mu_interval(self):
        r"""Abscissa limits for the :math:`\mu` coordinate in COM-dist [eV/T]:
        [a > 0, b], default=[0.0, 2.5e6].
        """
        return self._mu_interval

    @property
    def mu_bins(self):
        """Number of bins the ``mu_interval`` is divided to: (> 0), default=100.
        """
        return self._mu_bins

    @property
    def ekin_interval(self):
        r"""Abscissa limits for the :math:`E_{kin}` coordinate [eV]:
        [a > 0, b > 0], default=[0.0, 4.5e6].
        """
        return self._ekin_interval

    @property
    def ekin_bins(self):
        """Number of bins the ``ekin_interval`` is divided to: (> 0),
        default=50.
        """
        return self._ekin_bins

    @property
    def ptor_interval(self):
        r"""Abscissa limits for the :math:`P_{tor}` coordinate [kg m^2/s]:
        [a, b], default=[0.0, 1.0].
        """
        return self._ptor_interval

    @property
    def ptor_bins(self):
        """Number of bins the ``ptor_interval`` is divided to: (> 0),
        default=200.
        """
        return self._ptor_bins


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

    @property
    def poincare(self):
        """Collect data only when a (Poincaré) plane is crossed: {0, 1},
        default=0.

        Enable this option in order to generate Poincaré plots.
        """
        return self._poincare

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

    @property
    def poloidal_angles(self):
        """Poloidal angles of toroidal planes where toroidal Poincaré plots are
        collected [deg]: [0 <= a <= 360, ...], default=[0.0, 180.0].

        Used when ``poincare`` is enabled.
        """
        return np.asarray(self._poloidal_angles)

    @property
    def toroidal_angles(self):
        """Toroidal angles of poloidal planes where poloidal Poincaré plots are
        collected [deg]: [0 <= a <= 360, ...], default=[0.0, 180.0].

        Used when ``poincare`` is enabled.
        """
        return np.asarray(self._toroidal_angles)

    @property
    def radial_distances(self):
        """Minor radius coordinate where radial Poincaré plots are collected:
        [0 < a <= 1, ...], default=[1.0].

        Used when ``poincare`` is enabled.
        """
        return np.asarray(self._radial_distances).copy()

    @property
    def interval(self):
        """Time interval for writing marker state [s]: (> 0), default=0.

        Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 1.
        """
        return self._interval


@dataclass
class transport_coefficient():

    _margin: float = 0.0
    _record_rho_instead_of_r: bool = 0
    _number_of_points_to_average: int = 5
    _collect_transport_coefficient: bool = 0

    @property
    def collect_transport_coefficient(self):
        """Enable evaluation of transport coefficients: {0, 1}, default=0.

        - 0 Transport coefficients are not collected
        - 1 Transport coefficients are collected
        """
        return self._collect_transport_coefficient

    @property
    def margin(self):
        """Time interval for recording data points [s]: (> 0), default=0.0.

        The data points are recorded (at the outer mid-plane crossing) if this
        margin has passed from the previous recording.
        """
        return self._margin

    @property
    def number_of_points_to_average(self):
        """Number of subsequent data points that are averaged before calculating
        the coefficients to reduce noise: (> 0), default=5.
        """
        return self._number_of_points_to_average

    @property
    def record_rho_instead_of_r(self):
        """Record coefficients in terms of normalized poloidal flux instead of
        meters: {0, 1}, default=0.
        """
        return self._record_rho_instead_of_r
