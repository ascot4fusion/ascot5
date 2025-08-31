"""An interface between the user, the data objects, and the actual data.

The user interface is provided by `Tree` instance and `InputVariant` instances.
The latter two provide also the interface between the actual data (whose storage
is specified by `Format` and the data objects.
"""
from typing import Protocol

from .tree import Tree
from .treeparts import TreeManager
from .variants import InputVariant, RunVariant, Diagnostic
from .dataholder import Format


#pylint: disable=too-few-public-methods
class TreeCreateClassMixin(Protocol):
    """Base class for mixin classes that introduce create* methods.

    This class suppresses linter warnings by providing access to the
    `_treemanager` attribute of the child (main) class.
    """
    _treemanager: TreeManager

__all__  = [
    "Tree",
    "Format",
    "Diagnostic",
    "RunVariant",
    "InputVariant",
    "TreeCreateClassMixin",
    ]
