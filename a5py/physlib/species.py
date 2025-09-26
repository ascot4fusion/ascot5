"""Library of commonly encountered plasma species and their properties."""
from __future__ import annotations

from typing import NamedTuple

import unyt

plasma_species = {
    "e": (0, 0, 0.0005486),
    "n": (0, 1, 1.009),
    "H": (1, 1, 1.007),
    "D": (1, 2, 2.014),
    "T": (1, 3, 3.016),
    "He3": (2, 3, 3.016),
    "He4": (2, 4, 4.003),
    "Li6": (3, 6, 6.015),
    "Li7": (3, 7, 7.016),
    "Be9": (4, 9, 9.012),
    "B10": (5, 10, 10.012),
    "B11": (5, 11, 11.009),
    "C12": (6, 12, 12.011),
    "C13": (6, 13, 13.003),
    "N14": (7, 14, 14.003),
    "N15": (7, 15, 15.000),
    "O16": (8, 16, 15.994),
    "O17": (8, 17, 16.999),
    "O18": (8, 18, 17.999),
    "Ne20": (10, 20, 19.992),
    "Ne22": (10, 22, 21.991),
    "Ar36": (18, 36, 35.967),
    "Ar40": (18, 40, 39.962),
    "Ni59": (28, 59, 58.934),
    "Sn118": (50, 118, 117.901),
    "Xe132": (54, 132, 131.904),
    "W183": (74, 183, 182.950),
    "W184": (74, 184, 183.950),
    }
"""Names of the recognized species and their properties (Z, A, m)."""


species_aliases = {
    "e": ("electron", "Electron", "e-"),
    "n": ("neutron", "Neutron", "n0",),
    "H": ("hydrogen", "Hydrogen", "H1", "proton", "p",),
    "D": ("deuterium", "Deuterium", "H2",),
    "T": ("tritium", "Tritium", "H3",),
    "He3": ("helium-3", "Helium-3",),
    "He4": ("helium-4", "Helium-4", "alpha",),
}
"""Aliases for the plasma species."""


class Species(NamedTuple):
    """Named tuple containing species properties."""

    name : str
    """Name of the species."""

    znum: int
    """Charge number."""

    anum: int
    """Atomic mass number."""

    mass: unyt.unyt_quantity
    """Species mass."""

    @staticmethod
    def from_string(species: str) -> Species:
        """Get Species of given name.

        Available names are listed by `speciesdict`.

        Parameters
        ----------
        species : str
            Name of the species.

        Returns
        -------
        obj : Species
            Queried species.

        Raises
        ------
        ValueError
            If the species is not recognized.
        """
        for species_name, aliases in species_aliases.items():
            if species in aliases:
                species = species_name

        for name, properties in plasma_species.items():
            if name == species:
                return Species(
                    name, properties[0], properties[1], properties[2]*unyt.amu,
                    )
        raise ValueError(
            f"Unknown species: {species}. List of known species are: "
            f"{list(plasma_species.keys())}"
            )

    @staticmethod
    def from_znumanum(znum: int, anum: int) -> Species:
        """Get Species of given charge number and atomic mass number.

        Parameters
        ----------
        znum : int
            Charge number.
        anum : int
            Atomic mass number.

        Returns
        -------
        obj : Species
            Queried species.

        Raises
        ------
        ValueError
            If the species is not recognized.
        """
        for name, properties in plasma_species.items():
            if properties[0] == znum and properties[1] == anum:
                return Species(name, znum, anum, properties[2]*unyt.amu)
        raise ValueError(
            f"Unknown species: Z={znum} A={anum}. List of known species are: "
            f"{list(plasma_species.keys())}"
            )
