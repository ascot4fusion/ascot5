"""Defines simulation options input class and the corresponding factory method.
"""
from __future__ import annotations

import ctypes
from typing import Tuple, Optional, TypeVar
from dataclasses import fields

import inspect
import re

import unyt
import numpy as np

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils
from ...exceptions import AscotIOException
from a5py.engine.functions import END_CONDITIONS

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from a5py import Ascot


from .parameters import (
    simulation,
    physics,
    endconditions,
    distributions,
    comdistribution,
    orbit,
    transport_coefficient,
)

parameter_groups = [
    "simulation", "physics", "endconditions", "orbit", "distributions",
    "comdistribution", "transport_coefficient",
    ]
"""The groups in which options parameters are categorized."""

require_both_tor_and_pol = 2
"""Options setting that requires both toroidal and poloidal endconditions to
be met.
"""


class Options(InputVariant):
    """Simulation options."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in options.h."""
        MAXPOINCARE = 32
        _pack_ = 1
        _fields_ = [
            ("simulation_mode", ctypes.c_int),
            ("enable_adaptive", ctypes.c_int),
            ("record_mode", ctypes.c_int),
            ("use_explicit_fixedstep", ctypes.c_int),
            ("gyrodefined_fixedstep", ctypes.c_int),
            ('PADDING_0', ctypes.c_ubyte * 4),
            ("explicit_fixedstep", ctypes.c_double),
            ("adaptive_tolerance_orbit", ctypes.c_double),
            ("adaptive_tolerance_collisions", ctypes.c_double),
            ("adaptive_max_drho", ctypes.c_double),
            ("adaptive_max_dphi", ctypes.c_double),
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
            ('PADDING_1', ctypes.c_ubyte * 4),
            ("lab_time_limit", ctypes.c_double),
            ("max_mileage", ctypes.c_double),
            ("max_real_time", ctypes.c_double),
            ("rho_coordinate_limits", ctypes.c_double * 2),
            ("min_energy", ctypes.c_double),
            ("min_local_thermal_energy", ctypes.c_double),
            ("max_number_of_toroidal_orbits", ctypes.c_double),
            ("max_number_of_poloidal_orbits", ctypes.c_double),
            ("collect_orbit", ctypes.c_int),
            ("collect_dist5d", ctypes.c_int),
            ("collect_dist6d", ctypes.c_int),
            ("collect_dist5drho", ctypes.c_int),
            ("collect_dist6drho", ctypes.c_int),
            ("collect_distcom", ctypes.c_int),
            ("collect_transport_coefficient", ctypes.c_int),
            ("poincare", ctypes.c_int),
            ("number_of_points_per_marker", ctypes.c_int),
            ("ntoroidalplots", ctypes.c_int),
            ("npoloidalplots", ctypes.c_int),
            ("nradialplots", ctypes.c_int),
            ("interval", ctypes.c_double),
            ("toroidal_angles", ctypes.c_double * MAXPOINCARE),
            ("poloidal_angles", ctypes.c_double * MAXPOINCARE),
            ("radial_distances", ctypes.c_double * MAXPOINCARE),
            ("number_of_points_to_average", ctypes.c_int),
            ("record_rho_instead_of_r", ctypes.c_int),
            ("margin", ctypes.c_double),
            ("r_bins", ctypes.c_int),
            ("phi_bins", ctypes.c_int),
            ("z_bins", ctypes.c_int),
            ("ppara_bins", ctypes.c_int),
            ("pperp_bins", ctypes.c_int),
            ("time_bins", ctypes.c_int),
            ("charge_bins", ctypes.c_int),
            ('PADDING_2', ctypes.c_ubyte * 4),
            ("r_interval", ctypes.c_double * 2),
            ("phi_interval", ctypes.c_double * 2),
            ("z_interval", ctypes.c_double * 2),
            ("ppara_interval", ctypes.c_double * 2),
            ("pperp_interval", ctypes.c_double * 2),
            ("time_interval", ctypes.c_double * 2),
            ("charge_interval", ctypes.c_double * 2),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Options",
            struct=Options.Struct(),
            )
        self._orbit: orbit
        self._physics: physics
        self._simulation: simulation
        self._endconditions: endconditions
        self._distributions: distributions
        self._comdistribution: comdistribution
        self._transport_coefficient: transport_coefficient

    @property
    def simulation(self):
        """Parameters related to simulation mode and time-step."""
        return self._simulation

    @property
    def physics(self):
        """Physics included in the simulation."""
        return self._physics

    @property
    def endconditions(self):
        """End conditions when the marker simulation is ceased."""
        return self._endconditions

    @property
    def distributions(self):
        """Diagnostics that collect data for reproducing the particle
        distribution function."""
        return self._distributions

    @property
    def comdistribution(self):
        """Diagnostic that reproduces the particle distribution in
        constants-of-motion phase-space."""
        return self._comdistribution

    @property
    def orbit(self):
        """Diagnostic that records the exact marker trajectory."""
        return self._orbit

    @property
    def transport_coefficient(self):
        """Diagnostic that evaluates advection and diffusion coefficients for
        the radial transport of the simulated particle population."""
        return self._transport_coefficient


    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        """Return a dictionary with sufficient data to duplicate this instance.

        Returns
        -------
        data : dict[str, np.ndarray or unyt.unyt_array]
            Data that can be passed to create_btc to duplicate this instance.
        """
        def dataclass_to_dict(obj):
            collection = {}
            for f in fields(obj):
                no_underscore = f.name[1:]
                collection[no_underscore] = getattr(obj, no_underscore)
            return collection

        data = {}
        for param_group in parameter_groups:
            dataclass = getattr(self, param_group)
            data.update(dataclass_to_dict(dataclass))
        return data

    def stage(self):
        parameters = self.export()
        map_endcond_to_param = {
            "reached_time_limit":"activate_simulation_time_limits",
            "below_min_energy":"activate_energy_limits",
            "thermalized":"activate_energy_limits",
            "hit_wall":"activate_wall_hits",
            "below_rho_limit":"activate_rho_limit",
            "above_rho_limit":"activate_rho_limit",
            "completed_poloidal_orbits":"activate_orbit_limit",
            "completed_toroidal_orbits":"activate_orbit_limit",
            "simulation_not_finished":"activate_real_time_limit",
            "neutralized":"activate_neutralization",
            "ionized":"activate_ionization",
        }
        parameters["endcond_active"] = 0
        for ec, param in map_endcond_to_param.items():
            parameters["endcond_active"] += (
                END_CONDITIONS[ec] * parameters[param]
                )
        parameters.update({
            "charge_bins":np.abs(np.diff(parameters["charge_interval"])[0]) + 1,
            "nradialplots":parameters["radial_distances"].size,
            "ntoroidalplots":parameters["poloidal_angles"].size,
            "npoloidalplots":parameters["toroidal_angles"].size,
            "require_both_tor_and_pol":parameters["activate_orbit_limit"] == require_both_tor_and_pol,
            })
        for field, _ in self._struct_._fields_:
            if field in parameters:
                val = parameters[field]
                if isinstance(getattr(self._struct_, field), ctypes.Array):
                    arr = getattr(self._struct_, field)
                    for i, v in enumerate(val):
                        arr[i] = v
                else:
                    setattr(self._struct_, field, val)

    def unstage(self):
        pass

    def configure(self, desc=None, **params):
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

    def export_as_string(self, descriptions=True, aslist=False):
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
        def make_banner(title: str, desc: str, width: int = 80) -> str:
            border = "#" + "*" * (width - 2) + "#"
            empty = "#*".ljust(width - 2) + "*#"

            def center_line(text: str) -> str:
                content = text.center(width - 4)
                return f"#*{content}*#"

            return "\n".join([
                border,
                center_line(title),
                empty,
                center_line(desc),
                border,
                "\n"
            ])

        text = ""
        def clean_docstring(docstring):
            docstring = inspect.cleandoc(docstring) + "\n"
            docstring = "\n".join("# " + line for line in docstring.splitlines())
            docstring = re.sub(r"``([^`]+)``", lambda m: m.group(1), docstring)
            docstring = re.sub(r":math:`([^`]+)`", lambda m: m.group(1), docstring)
            if docstring.splitlines()[-1].strip() != "#":
                docstring += "\n#"
            return docstring

        for param_group in parameter_groups:
            title = param_group
            text += make_banner(title, getattr(Options, param_group).__doc__)
            for name, value in globals()[param_group].__dict__.items():
                if isinstance(value, property):
                    attribute = getattr(self, param_group)
                    text += clean_docstring(value.__doc__) + "\n"
                    text += name + " = " + str(getattr(attribute, name)) + "\n"
                    text += "\n"
            text += "\n"
        return text

    @classmethod
    def from_text(
        cls,
        ascot: Ascot,
        text: str | list[str],
        note: Optional[str] = None,
        activate: bool = False,
        dryrun: bool = False,
        store_hdf5: Optional[bool] = None,
        ) -> Options:
        """Read options parameter from a text.

        Parameters
        ----------
        text : str or [str]
            Parameters either as a list of strings or in a single string with
            one parameter per line.

            The string may contain empty lines or comment lines (starting with
            '#') as those are ignored. The parameters are expected to have the
            format <parameter name> = <value>. The string doesn't have to
            contain all possible parameters.
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        options : ~a5py.data.options.Options
            Freshly minted options input.
        """
        if not isinstance(text, list):
            text = text.splitlines()

        params_found = {}
        for line in text:
            is_parameter_line = (
                len(line) > 0 and
                not line.startswith("#") and
                "=" in line
            )
            if is_parameter_line:
                parameter, value = line.strip().split("=")
                params_found[parameter] = value

        return ascot.data.create_options(
            note=note, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            **params_found
            )


# pylint: disable=too-few-public-methods
class CreateOptionsMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create Options input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_options(
            self,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            **parameters,
            ) -> Options:
        r"""Create simulation options.

        This method creates a simulation options input.

        Parameters
        ----------
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.
        **parameters
            Options parameters and the values to be set.

            If a parameter is not set, the default value is used.

        Returns
        -------
        inputdata : ~a5py.data.options.Options
            Freshly minted input data object.
        """
        meta = variants.new_metadata("Options", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        parameters_found = []
        for paramgroup in parameter_groups:
            params_for_this_group = {}
            for k, v in parameters.items():
                if k in globals()[paramgroup].__dict__:
                    parameters_found.append(k)
                    params_for_this_group["_" + k] = v
            dataclass = globals()[paramgroup](**params_for_this_group)
            setattr(obj, "_" + paramgroup, dataclass)
        for k in parameters.keys():
            if k not in parameters_found:
                raise ValueError(f"Unknown options parameter '{k}'.")

        if store_hdf5:
            obj._export_hdf5()
        return obj
