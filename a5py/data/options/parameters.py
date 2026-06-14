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
from dataclasses import dataclass, field, fields
from types import MappingProxyType

import unyt
from a5py.data.options import constraints

from typing import Any
from a5py.physlib import match_units
import inspect
import re


class Lockable:
    """Object (dataclass) that can be frozen to prevent modification.

    Raises an exception if an attempt is made to modify a frozen object.
    
    Parameters
    ----------
    _frozen : bool
        True if this object is frozen.
    """
    _frozen: bool = False

    def _freeze(self):
        """Recursively freeze this object and all nested dataclasses."""
        for f in fields(self):
            value = getattr(self, f.name)

            if isinstance(value, Lockable):
                value._freeze()

        object.__setattr__(self, "_frozen", True)

    def __setattr__(self, name: str, value: Any):
        if getattr(self, "_frozen", False):
            raise AttributeError(
                f"{self.__class__.__name__} is frozen; cannot modify '{name}'"
            )
        super().__setattr__(name, value)

    def _attribute_doc(self, attr):
        """Return the docstring for a given attribute.
        
        Parameters
        ----------
        attr : str
            Attribute name.

        Returns
        -------
        str
            Attribute docstring.
        """
        doc = inspect.getdoc(self.__class__)
        if not doc:
            raise ValueError(f"No docstring found for {self.__class__}")

        pattern = re.compile(
            rf"^\s*{re.escape(attr)}\s*:\s*.*?\n"
            rf"((?:^(?:\s{{4,}}|\s*$).*\n?)*)",
            re.MULTILINE,
        )

        match = pattern.search(doc)
        if not match:
            raise KeyError(f"Attribute '{attr}' not found in docstring")

        body = inspect.cleandoc(match.group(1))
        return body


class ParameterClass(Lockable):
    """Base class for option parameters.
    
    Raises an exception if trying to modify a parameter with an invalid value or
    if the object is frozen.
    """

    def __setattr__(self, name, value):
        try:
            doc = self._attribute_doc(name)
        except KeyError:
            raise ValueError(f"Unknown parameter '{name}'.") from None

        valid = constraints.enforce(value, constraints.parse(doc)["constraint"])
        if not valid:
            class_name = self.__class__.__name__.lower()
            raise ValueError(
                f"Invalid value for {class_name}.{name}: {value}."
                )
        super().__setattr__(name, value)


@dataclass
class Simulation(ParameterClass):
    """Parameters related to simulation mode and time-step.

    Attributes
    ----------
    mode : str
        Simulation mode: {"gyro-orbit", "guiding-center", "hybrid",
        "field-line"}, default="gyro-orbit".

        - "gyro-orbit": Solve the whole gyro-motion.
        - "guiding-center": Trace the guiding center of the gyro-motion.
        - "hybrid": Hybrid mode where the simulation mode changes from
           guiding-center to gyro-orbit when the marker exits the plasma region.
        - "field-line": Trace magnetic field lines.

    record_mode : int
        Change the physical picture before collecting diagnostics: {0, 1},
        default=0.

        This option only affects the gyro-orbit simulations.

        - 0: Record gyro-orbits as they are.
        - 1: Instead, record the guiding center position of the gyro-orbit.

    timestep : float
        User-defined time-step [s]: (> 0), default=1.0e-8.

        This time-step is used as the fixed time step and as an initial step for
        adaptive time-stepping.

    enable_adaptive : bool
        Use adaptive time-step: {0, 1}, default=1.

        This option is used only if ``simulation_mode`` is 2 or 3. Gyro-orbit
        simulations are always done with fixed time-step and magnetic field line
        simulations with adaptive time-step. Note: The adaptive scheme uses
        fixed time-step value as an initial step.

        - 0: Use fixed time-step.
        - 1: Use adaptive time-step.

    adaptive_tolerance_orbit : float
        Relative error tolerance for orbit following in adaptive scheme:
        (> 0), default=1.0e-8.

    adaptive_tolerance_collisions : float
        Relative error tolerance for Coulomb collisions in adaptive scheme:
        (> 0), default=1.0e-1.
    """

    mode: str = "gyro-orbit"
    record_mode: int = 0
    timestep: float = 1.0e-9
    enable_adaptive: bool = True
    adaptive_tolerance_orbit: float = 1.0e-8
    adaptive_tolerance_collisions: float = 1.0e-1


@dataclass
class Physics(ParameterClass):
    """Physics included in the simulation.

    Attributes
    ----------
    enable_orbit_following : bool
        Trace markers in an electromagnetic field: {0, 1}, default=0.

    enable_coulomb_collisions : bool
        Markers experience Coulomb collisions with background plasma:
        {0, 1}, default=0.

    enable_mhd : bool
        Include MHD perturbations to orbit-following: {0, 1}, default=0.

    enable_atomic : bool
        Markers can undergo atomic reactions with background plasma
        or neutrals: {0, 1, 2}, default=0.

        - 0: Atomic reactions are turned off.
        - 1: Atomic reactions are on but marker is terminated when outside
             the reaction data domain.
        - 2: Atomic reactions are on but they are ignored when marker is
             outside the reaction data domain.

    enable_icrh : bool
        Enable ion cyclotron resonance heating operator: {0, 1}, default=0.

        ICRH operator transfers energy (via "kicks") to ions when they are on
        the resonance. The code must be compiled with RFOF=1 and the RFOF
        library must be present in order to use the ICRH operator.

    enable_aldforce : bool
        Enable radiation reaction force (synchrotron losses): {0, 1},
        default=0.

        The radiation reaction force (a.k.a. Abraham-Lorentz-Dirac or ALD force)
        causes charged particles to lose energy via radiation. The losses are
        proportional to the particle energy and inversely proportional to the
        particle mass, making this option mostly relevant for (runaway)
        electrons.

    disable_first_order_gctransformation : bool
        Disable first-order guiding center transformation: {0, 1}, default=0.

    disable_ccoll_gcenergy : bool
        Disable guiding center energy collisions: {0, 1}, default=0.

    disable_ccoll_gcpitch : bool
        Disable guiding center pitch angle collisions: {0, 1}, default=0.

    disable_ccoll_gcspatial : bool
        Disable guiding center spatial collisions: {0, 1}, default=0.

    reverse_time : bool
        Trace markers backwards in time: {0, 1}, default=0.

         Collision operator isn't reversible so disable collisions if this
         option is used. Also when tracing markers, the simulation stops when
         marker time is *below* ``simulation_time_limit``.
    """

    enable_orbit_following: bool = 0
    enable_coulomb_collisions: bool = 0
    enable_mhd: bool = 0
    enable_atomic: bool = 0
    enable_icrh: bool = 0
    enable_aldforce: bool = 0
    disable_first_order_gctransformation: bool = 0
    disable_ccoll_gcenergy: bool = 0
    disable_ccoll_gcpitch: bool = 0
    disable_ccoll_gcspatial: bool = 0
    reverse_time: bool = 0


@dataclass
class Endconditions(ParameterClass):
    """Conditions for terminating marker simulation.

    Attributes
    ----------
    activate_simulation_time_limits : bool
        Terminate marker based on laboratory time or its lifetime: {0, 1},
        default=0.

        Terminate when either of the following is true:
        - Laboratory time exceeds ``lab_time_limit``.
        - Marker lifetime exceeds ``max_mileage``.

    activate_real_time_limit : bool
        Terminate marker when the computer has spent specified
        amount of real time to simulate it: {0, 1}, default=0.

        This is not a "proper" end condition in a sense that it does not
        correspond to any physical process. This should be used just to control
        simulation duration or debugging. The limit is set by ``max_real_time``.

    activate_rho_limit : bool
        Terminate if marker goes outside given rho boundaries: {0, 1},
        default=0.

        The boundaries are defined by ``rho_coordinate_limits``.

    activate_energy_limits : bool
        Terminate when marker energy is below a user-specified value: {0, 1},
        default=0.

        The user specified values are ``min_energy`` and
        ``local_thermal_limit``. Marker is terminated when either
        of these limits is reached.

    activate_wall_hits : bool
        Terminate when marker intersects a wall element: {0, 1}, default=0.

    activate_orbit_limit : bool
        Terminate when marker has completed user-specified number of orbits:
        {"no", "either", "both"}, default="no".

        - "either": The number of toroidal and poloidal orbits is limited by
          ``max_number_of_toroidal_orbits`` or
          ``max_number_of_poloidal_orbits``, whichever is reached first.
        - "both": Both limits must be reached.

    activate_neutralization : bool
        Terminate when marker becomes neutral: {0, 1}, default=0.

    activate_ionization : bool
        Terminate when marker becomes ionized: {0, 1}, default=0.

    lab_time_limit : float
        Maximum laboratory time in seconds: (> 0), default=1.0.

    max_mileage : float
        The maximum amount of time this marker is simulated [s] or [m]:
        (> 0), default=1.0.

    max_real_time : float
        Maximum real time spent simulating a marker [s]: (> 0),
        default=3600.0

    rho_coordinate_limits : tuple[float, float]
        Minimum and maximum values for rho: [a >= 0, b > 0],
        default=[0., 1.0].

    min_energy : float
        Minimum energy [eV]: (> 0), default=1000.0.

    local_thermal_limit : float
        Minimum energy is local ion thermal energy times this value:
        (> 0), default=2.0.

    max_number_of_toroidal_orbits : float
        Maximum number of toroidal orbits: (> 0), default=100.

    max_number_of_poloidal_orbits : float
        Maximum number of poloidal orbits: (> 0), default=100.
    """

    activate_simulation_time_limits: bool = 0
    activate_real_time_limit: bool = 0
    activate_rho_limit: bool = 0
    activate_energy_limits: bool = 0
    activate_wall_hits: bool = 0
    activate_orbit_limit: str = "no"
    activate_neutralization: bool = 0
    activate_ionization: bool = 0
    lab_time_limit: float = 1.0
    max_mileage: float = 1.0
    max_real_time: float = 3600.0
    rho_coordinate_limits: tuple[float, float] = (0., 1.)
    min_energy: float = 1.0e3
    local_thermal_limit: float = 2.0
    max_number_of_toroidal_orbits: float = 100
    max_number_of_poloidal_orbits: float = 100


@dataclass
class OrbitParams(ParameterClass):
    """Diagnostic that records the exact marker trajectory.

    Attributes
    ----------
    collect : str
        Enable diagnostics that store marker orbit: {"no", "interval",
        "poincare"}, default="no".
    buffer_size : int
        Maximum number of orbit points stored per marker: (> 0),
        default=100.

        When the limit is exceeded, the oldest points are overwritten, so only
        the most recent positions are retained.
    interval : float
        Time interval for writing marker state [s]: (>= 0), default=0.

        Used when ``collect='interval'``.
    poloidal_angles : list[float] | tuple[float] | float
        Poloidal angles of toroidal planes where toroidal Poincaré plots are
        collected [deg]: [0 <= a <= 360, ...], default=[0.0, 180.0].

        Used when ``collect='poincare'``.
    toroidal_angles : list[float] | tuple[float] | float
        Toroidal angles of poloidal planes where poloidal Poincaré plots are
        collected [deg]: [0 <= a <= 360, ...], default=[0.0, 180.0].

        Used when ``collect='poincare'``.
    radial_distances : list[float] | tuple[float] | float
        Minor radius coordinate where radial Poincaré plots are collected:
        [0 < a <= 1, ...], default=[1.0].

        Used when ``collect='poincare'``.
    """

    collect: str = "no"
    buffer_size: int = 100
    interval: float = 0.0
    poloidal_angles: list[float] | tuple[float] | float = (0.0,)
    toroidal_angles: list[float] | tuple[float] | float = (0.0,)
    radial_distances: list[float] | tuple[float] | float = (1.0,)


DIMENSIONS: MappingProxyType = MappingProxyType({
    "r": "m",
    "phi": "rad",
    "z": "m",
    "rho": "1",
    "theta": "rad",
    "ppar": "kg*m/s",
    "pperp": "kg*m/s",
    "pr": "kg*m/s",
    "pphi": "kg*m/s",
    "pz": "kg*m/s",
    "ekin": "J",
    "pitch": "1",
    "mu": "J/T",
    "ptor": "kg*m*m/s",
    "time": "s",
    "charge": 1,
})
"""All possible coordinate dimensions for the histogram diagnostics."""

@dataclass(frozen=True)
class Dimension():
    r"""Histogram axis.

    Attributes
    ----------
    name : str
        Name of the dimension.
    min : float
        Minimum value of the axis interval: [x].
    max : float
        Maximum value of the axis interval: [x].
    bins : int
        Number of bins the interval [``min``, ``max``] is divided to: (> 0).
    """
    name: str
    min: unyt.unyt_array
    max: unyt.unyt_array
    bins: int

    def __post_init__(self):
        if self.name not in DIMENSIONS or self.name == "charge":
            raise ValueError(
                f"Invalid histogram dimension {self.name}."
                )

        if self.bins <= 0:
            raise ValueError(
                f"Negative number of bins in dimension {self.name}."
                )
        if self.bins != int(self.bins):
            raise ValueError(
                f"Non-integer number of bins in dimension {self.name}."
                )

        object.__setattr__(self, "min", match_units(self.min, DIMENSIONS[self.name]))
        object.__setattr__(self, "max", match_units(self.max, DIMENSIONS[self.name]))

        if self.min >= self.max:
            raise ValueError(
                "Minimum value of the interval is larger than the maximum value"
                f" in dimension {self.name}."
                )


@dataclass
class HistParams(Lockable):
    r"""A histogram.

    Attributes
    ----------
    dimensions : list
        Set of distribution abscissae coordinates.

        The available coordinates are:

        - 'r' (:math:`r`): Major radius [m].
        - 'phi' (:math:`\phi`): Toroidal angle [def].
        - 'z' (:math:`z`): Axial height [m].
        - 'rho' (:math:`\rho`): Square root of normalized poloidal flux [1].
        - 'theta' (:math:`\theta`): Geometrical poloidal angle [def].
        - 'ppara' (:math:`p_\parallel`): Momentum component parallel to magnetic
          field [kg m/s].
        - 'pperp' (:math:`p_\perp`): Momentum component perpendicular to magnetic
          field [kg m/s].
        - 'pr' (:math:`p_r`): Radial component of momentum [kg m/s].
        - 'pphi' (:math:`p_\phi`): Toroidal component of momentum [kg m/s].
        - 'pz' (:math:`p_z`): Axial component of momentum [kg m/s].
        - 'ekin' (:math:`E_{kin}`): Kinetic energy [eV].
        - 'pitch' (:math:`\lambda`): Pitch angle [1].
        - 'mu' (:math:`\mu`): Magnetic moment [eV/T].
        - 'ptor' (:math:`P_{tor}`): Canonical toroidal angular momentum
          [kg m**2/s].
        - 'time' (:math:`t`): Time [s].

    charge_interval : set
        Histogram limits (inclusive) for test particle charge [e]: [a, b],
        default=[2, 2].

        For each charge state on the interval, the distribution is collected
        separately.
    """
    dimensions: tuple[Dimension] = field(default_factory=tuple[Dimension])
    charge_interval: tuple[int, int] = (2, 2)

    def __setattr__(self, name, value):
        if name == "dimensions":
            arr = []
            for dim in value:
                if isinstance(dim, dict):
                    arr.append(Dimension(**dim))
                elif isinstance(dim, Dimension):
                    arr.append(dim)
                else:
                    arr.append(Dimension(*dim))

            if len(value) != len(set(d.name for d in arr)):
                raise ValueError(
                    "Multiple dimensions with the same name."
                    )
            object.__setattr__(self, name, arr)
        elif name == "charge_interval":
            try:
                doc = self._attribute_doc(name)
            except KeyError:
                raise ValueError(f"Unknown parameter '{name}'.") from None

            valid = constraints.enforce(value, constraints.parse(doc)["constraint"])
            if not valid:
                class_name = self.__class__.__name__.lower()
                raise ValueError(
                    f"Invalid value for {class_name}.{name}: {value}."
                    )
            object.__setattr__(self, name, value)
        else:
            raise ValueError(
                f"Uknown parameter {name}."
                )

    @classmethod
    def from_dict(cls, params):
        """Initialize parameters from a dictionary.
        
        Parameters
        ----------
        params : dict
            Dictionary of parameters.

            ``charge_interval`` is optional but there must be at least one
            dimension.
        """
        hist = cls()
        hist.dimensions = params["dimensions"]
        if len(hist.dimensions) == 0:
            raise ValueError(
                "Histogram must have at least one dimension."
                )
        if "charge_interval" in params:
            hist.charge_interval = params["charge_interval"]
        return hist
    
    @classmethod
    def from_hdf5(cls, file, index):
        hist = {"dimensions": []}
        try:
            value = file.read(f"hist_{index}__charge_interval")
            hist["charge_interval"] = (int(value[0]), int(value[1]))
        except KeyError:
            return None

        for dim in DIMENSIONS.keys():
            try:
                hist["dimensions"].append((
                    dim,
                    file.read(f"hist_{index}_{dim}_min"),
                    file.read(f"hist_{index}_{dim}_max"),
                    int(file.read(f"hist_{index}_{dim}_bins")),
                ))
            except KeyError:
                continue
        return hist
