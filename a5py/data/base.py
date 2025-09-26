"""Contains the main tree-like data structure class.

The tree maintains the tree structure and manages the metadata within the Tree
object itself, input categories, and input and output variants. Only the latter
two access the actual data directly, which is why they are also referred as
"data".
"""
from typing import Optional

from .access import Tree
from . import (
    options, marker, bfield, efield, plasma, wall,
    #efield, plasma, neutral, wall, mhd, boozer, nbi, asigma, options,
    )

input_categories = {
    "options", "bfield", "efield", "marker", "plasma", "neutral", "wall",
    "boozer", "mhd", "asigma", "nbi",
}

class AscotData(
    Tree, bfield.CreateBfieldMixin,
    efield.CreateEfieldMixin,
    plasma.CreatePlasmaMixin,
    # neutral.CreateNeutralMixin,
    wall.CreateWallMixin,
    marker.CreateMarkerMixin,
    #mhd.CreateMhdMixin, boozer.CreateBoozerMixin,
    #nbi.CreateNbiMixin, asigma.CreateAsigmaLocMixin,
    options.CreateOptionsMixin,
    ):
    """Stores and manages simulation inputs and outputs.

    The tree-like structure in which data is accessed looks like this::

        Data
        ├── bfield
        │   ├── BfieldCartesian_1
        |   |   ├── bxyz
        │   │   └── Other *datasets*
        │   └── Other input *variants* in <variant>_<id> format
        ├── Other input *categories*
        ├── run_1
        │   ├── inistate
        │   │   ├── weight
        │   │   └── Other *datasets*
        │   └── Other *diagnostics*
        └── Other run *variants* in <variant>_<id> format

    Data is referred to as input or output *variants*. Inputs are divided to
    different *categories* which act akin to interfaces: simulations and some
    methods assume that there's some data present in requested categories but
    the actual implementation is defined by the variant.

    This class manages the overall data structure, reading the data from file,
    and creation of new inputs.
    """

    def __init__(self, hdf5file: Optional[tuple[str, bool]]=None):
        super().__init__(input_categories=input_categories, hdf5file=hdf5file)
