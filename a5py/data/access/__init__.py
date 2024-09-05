"""An interface between the user, the data objects, and the actual data.

The user interface is provided by `Tree` instance and `InputVariant` instances.
The latter two provide also the interface between the actual data (whose storage
is specified by `Format` and the data objects.
"""
from .tree import Tree
from .variants import InputVariant

from .dataholder import Format

__all__  = ["Tree", "InputVariant", "Format"]
