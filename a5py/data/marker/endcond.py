from typing import Literal

import numpy as np

from a5py.engine.functions import get_enum_value


EndConds = Literal[
    "reached_time_limit",
    "below_min_energy",
    "thermalized",
    "hit_wall",
    "below_rho_limit",
    "above_rho_limit",
    "completed_poloidal_orbits",
    "completed_toroidal_orbits",
    "simulation_not_finished",
    "finished_gc_in_hybrid_mode",
    "neutralized",
    "ionized",
    ]


endconds = {
    k:get_enum_value(v) for k, v in {
        "reached_time_limit": "endcond_tlim",
        "below_min_energy": "endcond_emin",
        "thermalized": "endcond_therm",
        "hit_wall": "endcond_wall",
        "below_rho_limit": "endcond_rhomin",
        "above_rho_limit": "endcond_rhomax",
        "completed_poloidal_orbits": "endcond_polmax",
        "completed_toroidal_orbits": "endcond_tormax",
        "simulation_not_finished": "endcond_cpumax",
        "finished_gc_in_hybrid_mode": "endcond_hybrid",
        "neutralized": "endcond_neutr",
        "ionized": "endcond_ioniz",
        }.items()}
"""Mapping from end condition name to end condition code."""


def endcond_check(bitarr, string):
    """Check if the binary end condition matches the human-readable.

    Parameters
    ----------
    bitarr : :obj:`np.uint32`
        Value of the end condition as it is stored in the HDF5 file.
    string : str
        Human-readable end condition (case-insensitive).

        Markers may have multiple end conditions active simultaneously.
        If just the name of the end condition e.g. "MAXPOL" is passed,
        then all markers with the ``MAXPOL`` end condition are returned.

        If the end cond is preceded by "NOT", e.g. "NOT MAXPOL", then
        markers that don't have that end condition are returned.

        Passing multiple end conditions returns markers that have all listed
        end conditions active, e.g. "MAXPOL MAXTOR" returns markers that
        have both ``MAXPOL`` and ``MAXTOR`` active simultaneously.

    Returns
    -------
    match : bool
        True if the two representations of end conditions match.
    """
    endconds = string.upper().split()
    ec_yes = np.array(0, dtype=np.uint32)
    ec_non = np.array(0, dtype=np.uint32)

    i = 0
    while i < len(endconds):
        ec = endconds[i]

        NOT = False
        if ec == "NOT":
            NOT = True
            i += 1
            ec = endconds[i]

        try:
            ec = getattr(State, "_" + ec)
        except AttributeError:
            raise ValueError("Unknown end condition: " + ec)

        i += 1
        if NOT:
            ec_non = ec_non | ec
        else:
            ec_yes = ec_yes | ec

    return (bitarr & ec_yes) == ec_yes and (bitarr & ec_non) == 0


def endcond_tostring(bitarr):
    """Convert end condition bitarray to human readable.

    Parameters
    ----------
    bitarr : :obj:`np.uint32`
        Value of the end condition as it is stored in the HDF5 file.

    Returns
    -------
    string : str
        End condition in a human readable format.
    """
    string = ""
    for ec in endcond:
        if bitarr & getattr(State, "_" + ec):
            string += ec if string == "" else " and " + ec
    return string
