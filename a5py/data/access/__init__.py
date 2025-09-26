"""An interface between the user, the data objects, and the actual data.

The user interface is provided by :class:`.Tree` instance, which manages
the overall data, and :class:`.InputVariant` and :class:`.OutputVariant` which
manage the actual datasets whil themselves being managed by :class:`.Tree`.

.. rubric:: Public

.. autosummary::
    :nosignatures:

    ~leaves.Status
    ~nodes.Tree
    ~leaves.Leaf
    ~leaves.InputVariant
    ~leaves.OutputVariant

a5py.data.access.leaves
***********************

.. automodule:: a5py.data.access.leaves

.. autoclass:: a5py.data.access.leaves.Status
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.leaves.Leaf
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.leaves.InputVariant
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.leaves.OutputVariant
    :members:
    :undoc-members:
    :private-members:

a5py.data.access.nodes
**********************

.. autoclass:: a5py.data.access.nodes.ImmutableNode
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.nodes.InputCategory
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.nodes.Tree
    :members:
    :undoc-members:
    :private-members:

a5py.data.access.hdf5io
***********************

.. automodule:: a5py.data.access.hdf5io

.. autodata:: a5py.data.access.hdf5io.RESULTGROUP

.. autoclass:: a5py.data.access.hdf5io.TreeFileManager
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.hdf5io.DataAccess
    :members:
    :undoc-members:
    :private-members:

.. autoclass:: a5py.data.access.hdf5io.TreeFile
    :members:
    :undoc-members:
    :private-members:

a5py.data.access.tree
*********************

.. autoclass:: a5py.data.access.tree.TreeManager
    :members:
    :undoc-members:
    :private-members:

"""
from typing import Protocol

from .tree import TreeManager
from .nodes import Tree
from .leaves import Leaf, InputVariant, OutputVariant, Status

class TreeMixin(Protocol):
    """Supplies the tree manager for type checking.

    This is used by the mixin classes that make up the :class:`.AscotData`
    class that inherits from :class:`.Tree`.
    """

    _treemanager: TreeManager


__all__  = [
    "Tree",
    "Leaf",
    "Status",
    "InputVariant",
    "OutputVariant",
    ]
