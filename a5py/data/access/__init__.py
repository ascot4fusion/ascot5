"""An interface between the user, the data objects, and the actual data.

The user interface is provided by `Tree` instance and `InputVariant` instances.
The latter two provide also the interface between the actual data (whose storage
is specified by `Format` and the data objects.
"""
from .leaves import Leaf, InputVariant, OutputVariant, Status, DataStruct
from .nodes import Tree


class TreeManager():
    pass

class Format():
    pass

#pylint: disable=too-few-public-methods
class CreateLeafMixin():
    """Base class for mixin classes that introduce ``create*`` methods.

    This class suppresses linter warnings by providing access to the
    ``_treemanager`` attribute of the child (main) class.
    """
    _treemanager: TreeManager


TreeCreateClassMixin = CreateLeafMixin
RunVariant = OutputVariant

class Diagnostic():
    pass

__all__  = [
    "Tree",
    "Leaf",
    "Status",
    "DataStruct",
    "InputVariant",
    "OutputVariant",
    "Format",
    "CreateLeafMixin",
    ]
