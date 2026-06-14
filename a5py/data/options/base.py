"""Defines simulation parameters."""

from __future__ import annotations

import copy
import ctypes
from pathlib import Path
import textwrap
import toml
import tomli_w
import numpy as np
from typing import Optional
from dataclasses import fields, field, dataclass, asdict

from a5py.engine.functions import get_endcond

from a5py.libascot import DataStruct

from .parameters import (
    Simulation,
    Physics,
    Endconditions,
    HistParams,
    OrbitParams,
    Lockable,
)


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in options.h."""

    _fields_ = [
        ("simulation_mode", ctypes.c_int),
        ("enable_adaptive", ctypes.c_int),
        ("record_mode", ctypes.c_int),
        ("timestep", ctypes.c_double),
        ("adaptive_tolerance_orbit", ctypes.c_double),
        ("adaptive_tolerance_collisions", ctypes.c_double),
        ("enable_orbit_following", ctypes.c_int),
        ("enable_coulomb_collisions", ctypes.c_int),
        ("enable_mhd", ctypes.c_int),
        ("enable_atomic", ctypes.c_int),
        ("enable_icrh", ctypes.c_int),
        ("enable_aldforce", ctypes.c_int),
        ("disable_first_order_gctransformation", ctypes.c_int),
        ("disable_ccoll_gcenergy", ctypes.c_int),
        ("disable_ccoll_gcpitch", ctypes.c_int),
        ("disable_ccoll_gcspatial", ctypes.c_int),
        ("reverse_time", ctypes.c_int),
        ("endcond_active", ctypes.c_int),
        ("require_both_tor_and_pol", ctypes.c_int),
        ("lab_time_limit", ctypes.c_double),
        ("max_mileage", ctypes.c_double),
        ("max_real_time", ctypes.c_double),
        ("rho_coordinate_limits", ctypes.c_double * 2),
        ("min_energy", ctypes.c_double),
        ("local_thermal_limit", ctypes.c_double),
        ("max_number_of_toroidal_orbits", ctypes.c_double),
        ("max_number_of_poloidal_orbits", ctypes.c_double),
    ]

    @classmethod
    def from_params(cls, params: SimulationOptions):
        map_endcond_to_param = {
            "reached_time_limit": "activate_simulation_time_limits",
            "below_min_energy": "activate_energy_limits",
            "thermalized": "activate_energy_limits",
            "hit_wall": "activate_wall_hits",
            "below_rho_limit": "activate_rho_limit",
            "above_rho_limit": "activate_rho_limit",
            "completed_poloidal_orbits": "activate_orbit_limit",
            "completed_toroidal_orbits": "activate_orbit_limit",
            "simulation_not_finished": "activate_real_time_limit",
            "neutralized": "activate_neutralization",
            "ionized": "activate_ionization",
        }

        endcond_active = 0
        for ec, param in map_endcond_to_param.items():
            endcond = getattr(params.endconditions, param)
            if param == "activate_orbit_limit":
                endcond = endcond != "no"
            endcond_active += get_endcond(ec) * endcond
        require_both_tor_and_pol = int(
            params.endconditions.activate_orbit_limit == "both"
        )

        cdata = Struct()
        cdata.simulation_mode = [
            "gyro-orbit",
            "guiding-center",
            "hybrid",
            "field-line",
        ].index(params.simulation.mode) + 1
        cdata.enable_adaptive = params.simulation.enable_adaptive
        cdata.record_mode = params.simulation.record_mode
        cdata.timestep = params.simulation.timestep
        cdata.adaptive_tolerance_orbit = params.simulation.adaptive_tolerance_orbit
        cdata.adaptive_tolerance_collisions = (
            params.simulation.adaptive_tolerance_collisions
        )
        cdata.enable_orbit_following = params.physics.enable_orbit_following
        cdata.enable_coulomb_collisions = params.physics.enable_coulomb_collisions
        cdata.enable_mhd = params.physics.enable_mhd
        cdata.enable_atomic = params.physics.enable_atomic
        cdata.enable_icrh = params.physics.enable_icrh
        cdata.enable_aldforce = params.physics.enable_aldforce
        cdata.disable_first_order_gctransformation = (
            params.physics.disable_first_order_gctransformation
        )
        cdata.disable_ccoll_gcenergy = params.physics.disable_ccoll_gcenergy
        cdata.disable_ccoll_gcpitch = params.physics.disable_ccoll_gcpitch
        cdata.disable_ccoll_gcspatial = params.physics.disable_ccoll_gcspatial
        cdata.reverse_time = params.physics.reverse_time
        cdata.endcond_active = endcond_active
        cdata.require_both_tor_and_pol = require_both_tor_and_pol
        cdata.lab_time_limit = params.endconditions.lab_time_limit
        cdata.max_mileage = params.endconditions.max_mileage
        cdata.max_real_time = params.endconditions.max_real_time
        cdata.rho_coordinate_limits[0] = params.endconditions.rho_coordinate_limits[0]
        cdata.rho_coordinate_limits[1] = params.endconditions.rho_coordinate_limits[1]
        cdata.min_energy = params.endconditions.min_energy
        cdata.local_thermal_limit = params.endconditions.local_thermal_limit
        cdata.max_number_of_toroidal_orbits = (
            params.endconditions.max_number_of_toroidal_orbits
        )
        cdata.max_number_of_poloidal_orbits = (
            params.endconditions.max_number_of_poloidal_orbits
        )
        return cdata


@dataclass
class SimulationOptions(Lockable):
    """Simulation options.

    Attributes
    ----------
    simulation : Simulation
        Parameters related to simulation mode and time-step.
    physics : Physics
        Physics included in the simulation.
    endconditions : Endconditions
        End conditions when the marker simulation is ceased.
    orbit : Orbit
        Diagnostic that records the exact marker trajectory.
    histograms : tuple[HistParams]
        Diagnostics that collect data for reproducing the particle distribution function.
    """

    simulation: Simulation = field(default_factory=Simulation)
    physics: Physics = field(default_factory=Physics)
    endconditions: Endconditions = field(default_factory=Endconditions)
    orbit: OrbitParams = field(default_factory=OrbitParams)
    histograms: tuple[HistParams] = field(default_factory=tuple[HistParams])
    

    def __setattr__(self, name, value):
        if name == "histograms":
            arr = []
            for v in value:
                if isinstance(v, HistParams):
                    arr.append(v)
                elif isinstance(v, dict):
                    arr.append(HistParams.from_dict(v))
                else:
                    raise TypeError(
                        f"Invalid type for histograms: {type(v)}. "
                        "Must be either a HistParam or a dictionary."
                        )
            return super().__setattr__("histograms", tuple(arr))

        try:
            self._attribute_doc(name)
        except KeyError:
            raise ValueError(f"Unknown parameter '{name}'.") from None

        if (
            (name == "simulation" and not isinstance(value, Simulation)) or
            (name == "physics" and not isinstance(value, Physics)) or
            (name == "endconditions" and not isinstance(value, Endconditions)) or
            (name == "orbit" and not isinstance(value, OrbitParams))
            ):
            raise TypeError(f"Invalid type for {name}: {type(value)}.")

        super().__setattr__(name, value)


    def _write_hdf5(self, file):
        """Write parameters to file."""
        params = asdict(self)

        for group, settings in params.items():
            if group == "histograms":
                for i, hist in enumerate(params["histograms"]):
                    file.write(
                        f"hist_{i}__charge_interval",
                        np.asanyarray(hist["charge_interval"]),
                        )
                    for dim in hist["dimensions"]:
                        name = dim["name"]
                        file.write(
                            f"hist_{i}_{name}_min", np.asanyarray(dim["min"]),
                            )
                        file.write(
                            f"hist_{i}_{name}_max", np.asanyarray(dim["max"]),
                            )
                        file.write(
                            f"hist_{i}_{name}_bins", np.asanyarray(dim["bins"]),
                            )

                continue

            for setting, value in settings.items():
                if isinstance(value, str):
                    value = np.bytes_(value)
                value = np.asanyarray(value)
                file.write(f"{group}__{setting}", value)

    def copy(self) -> "SimulationOptions":
        """Make a deep copy of the parameters."""
        return copy.deepcopy(self)

    def write_toml(self, path: str | Path, include_comments="minimal") -> None:
        """Write parameters to TOML file.

        Parameters
        ----------
        path : str
            Path to the TOML file.
        include_comments : {"no", "minimal", "detailed"}, optional
            Whether to include descriptions of the parameters in the output.

            - "no": No comments.
            - "minimal": Include only the single-line description for each
              parameter.
            - "detailed": Include the full description for each parameter.
        """
        path = Path(path)

        def format_docstring(doc):
            """Convert docstring into TOML comments."""

            doc = textwrap.dedent(doc).strip()

            lines = []
            if include_comments == "minimal":
                parts = doc.split("\n\n", 1)
                doc = parts[0].replace("\n", " ")
            else:
                lines.append("")

            for line in doc.splitlines():

                if line.strip() == "":
                    lines.append("#")
                else:
                    lines.append(f"# {line}")

            return lines

        params = asdict(self)
        for hist in params["histograms"]:
            for dim in hist["dimensions"]:
                dim["min"] = float(dim["min"])
                dim["max"] = float(dim["max"])
 
        if len(params["histograms"]) == 0:
            del params["histograms"]
        toml_text = tomli_w.dumps(params)
        current_section = None

        out, wrote_histogram_docs = [], False
        for line in toml_text.splitlines():

            stripped = line.strip()

            if stripped.startswith("[") and stripped.endswith("]"):
                section = stripped[1:-1]
                current_section = section

                if include_comments == "detailed" and section != "[histograms]":
                    section_doc = self._attribute_doc(section)
                    if section_doc:
                        out.append(f"# {section_doc}")

                if (include_comments == "detailed" and section == "[histograms]"
                    and not wrote_histogram_docs):
                    section_doc = self._attribute_doc("histograms")
                    out.append(f"# {section_doc[:-1]}")
                    out.append("#")
                    import inspect
                    doc = inspect.cleandoc(inspect.getdoc(HistParams)).split("----------")[-1]
                    out.extend(format_docstring(doc)[1:])
                    wrote_histogram_docs = True

                out.append(line)
                continue

            if "=" in line and current_section:
                key = line.split("=", 1)[0].strip()

                if include_comments != "no" and current_section != "[histograms]":
                    group = getattr(self, current_section)
                    key_doc = group._attribute_doc(key)
                    out.extend(format_docstring(key_doc))

            out.append(line)

        final_text = "\n".join(out)

        with path.open("wb") as f:
            f.write(final_text.encode("utf-8"))

    @classmethod
    def from_dict(
        cls: "SimulationOptions", **params: dict[str, dict]
    ) -> "SimulationOptions":
        """Initialize parameters from keyword arguments (dictionary).

        Parameters
        ----------
        **params : dict
            Name of the parameter group and dictionary of settings related to
            that group.

        Returns
        -------
        SimulationOptions
            Initialized parameters.
        """
        opt = cls()
        for group, settings in params.items():
            if group == "histograms":
                arr = [HistParams.from_dict(hist) for hist in settings]
                opt.histograms = tuple(arr)
                continue
            if group not in opt.__dict__:
                raise ValueError(f"Unknown parameter group: {group}")

            for setting, value in settings.items():
                if setting not in getattr(opt, group).__dict__:
                    raise ValueError(f"Unknown parameter in group {group}: {setting}")
                setattr(getattr(opt, group), setting, value)
        return opt

    @classmethod
    def from_toml(cls: "SimulationOptions", path: str | Path) -> "SimulationOptions":
        """Initialize parameters from TOML file.

        Parameters
        ----------
        filename : str
            Name of the TOML file.

        Returns
        -------
        SimulationOptions
            Initialized parameters.
        """
        # path = Path(path)
        with open(path, "r") as f:
            data = toml.load(f)

        return SimulationOptions.from_dict(**data)

    @classmethod
    def from_json(cls: "SimulationOptions", filename: str) -> "SimulationOptions":
        return cls()

    @classmethod
    def from_hdf5(cls: "SimulationOptions", file) -> "SimulationOptions":
        """Read parameters from HDF5 file."""
        params = {
            "simulation": {}, "physics": {}, "endconditions": {}, "orbit": {},
            "histograms": [],
            }
        def read_parameter(group, setting, convert=None):
            try:
                value = file.read(f"{group}__{setting}")
            except KeyError:
                # Todo: implement any backwards compatibility here
                raise ValueError(f"Unknown parameter: {group}__{setting}") from None
            if convert == "tuple":
                value = tuple(value)
            if convert == "str":
                value = np.astype(value, "bytes_").decode("utf-8")
            params[group][setting] = value


        group = "simulation"
        read_parameter(group, "mode", convert="str")
        read_parameter(group, "record_mode")
        read_parameter(group, "timestep")
        read_parameter(group, "enable_adaptive")
        read_parameter(group, "adaptive_tolerance_orbit")
        read_parameter(group, "adaptive_tolerance_collisions")

        group = "physics"
        read_parameter(group, "enable_orbit_following")
        read_parameter(group, "enable_coulomb_collisions")
        read_parameter(group, "enable_mhd")
        read_parameter(group, "enable_atomic")
        read_parameter(group, "enable_icrh")
        read_parameter(group, "enable_aldforce")
        read_parameter(group, "disable_first_order_gctransformation")
        read_parameter(group, "disable_ccoll_gcenergy")
        read_parameter(group, "disable_ccoll_gcpitch")
        read_parameter(group, "disable_ccoll_gcspatial")
        read_parameter(group, "reverse_time")

        group = "endconditions"
        read_parameter(group, "activate_simulation_time_limits")
        read_parameter(group, "activate_real_time_limit")
        read_parameter(group, "activate_rho_limit")
        read_parameter(group, "activate_energy_limits")
        read_parameter(group, "activate_wall_hits")
        read_parameter(group, "activate_orbit_limit", convert="str")
        read_parameter(group, "activate_neutralization")
        read_parameter(group, "activate_ionization")
        read_parameter(group, "lab_time_limit")
        read_parameter(group, "max_mileage")
        read_parameter(group, "max_real_time")
        read_parameter(group, "rho_coordinate_limits", convert="tuple")
        read_parameter(group, "min_energy")
        read_parameter(group, "local_thermal_limit")
        read_parameter(group, "max_number_of_toroidal_orbits")
        read_parameter(group, "max_number_of_poloidal_orbits")

        group = "orbit"
        read_parameter(group, "collect", convert="str")
        read_parameter(group, "buffer_size")
        read_parameter(group, "interval")
        read_parameter(group, "poloidal_angles", convert="tuple")
        read_parameter(group, "toroidal_angles", convert="tuple")
        read_parameter(group, "radial_distances", convert="tuple")

        i = 0
        while True:
            hist = HistParams.from_hdf5(file, i)
            if hist is None:
                break

            params["histograms"].append(hist)
            i += 1

        return cls.from_dict(**params)
