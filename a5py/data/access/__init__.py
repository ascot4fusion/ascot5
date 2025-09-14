"""An interface between the user, the data objects, and the actual data.

The user interface is provided by `Tree` instance and `InputVariant` instances.
The latter two provide also the interface between the actual data (whose storage
is specified by `Format` and the data objects.
"""
from typing import Protocol

from .leaf import Leaf, InputLeaf
from .nodes import OutputLeaf
from .tree import Tree, TreeManager
from .dataholder import Format

from .variants import InputVariant, RunVariant


#pylint: disable=too-few-public-methods
class CreateLeafMixin(Protocol):
    """Base class for mixin classes that introduce ``create*`` methods.

    This class suppresses linter warnings by providing access to the
    ``_treemanager`` attribute of the child (main) class.
    """
    _treemanager: TreeManager


TreeCreateClassMixin = CreateLeafMixin

class Diagnostic():
    pass

__all__  = [
    "Tree",
    "Leaf",
    "InputLeaf",
    "OutputLeaf",
    "Format",
    "CreateLeafMixin",
    ]
