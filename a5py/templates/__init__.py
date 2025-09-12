"""Package for generating inputs from templates and imported data.
"""

from .analyticalinputs import (
    PremadeMagneticField,
    FlatPlasma,
    )
from .optionsxml import make_simple_type, make_element_block, doc, make_schema

def create_input(ascot, template, note=None, activate=False, dryrun=False, store_hdf5=False):
    """Create an input object from the template output."""
    factory = getattr(ascot, "create_" + template[0])
    return factory(**template[1], note=note, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5)

__all__ = [
    "PremadeMagneticField",
    "FlatPlasma",
    ]
