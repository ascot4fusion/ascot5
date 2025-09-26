"""This module contains all marker input variants and combines their
factory methods to a single class.

.. autosummary::
    :nosignatures:

    ~ParticleMarker
    ~GuidingcenterMarker
    ~FieldlineMarker
    ~CreateMarkerMixin.create_particlemarker
    ~CreateMarkerMixin.create_guidingcentermarker
    ~CreateMarkerMixin.create_fieldlinemarker

.. rubric:: Classes

.. autoclass:: ParticleMarker
    :members:

.. autoclass:: GuidingcenterMarker
    :members:

.. autoclass:: FieldlineMarker
    :members:

.. autoclass:: CreateMarkerMixin
    :members:
    :inherited-members:
"""
from . import particle
from . import fieldline
from . import guidingcenter
from .particle import ParticleMarker
from .fieldline import FieldlineMarker
from .guidingcenter import GuidingcenterMarker

# pylint: disable=too-many-ancestors
class CreateMarkerMixin(
    fieldline.CreateMixin,
    particle.CreateMixin,
    guidingcenter.CreateMixin,
    ):
    """Mixin class used by :class:`.AscotData` to create marker input.

    This class just combines all the marker mixin classes.
    """

__all__  = [
    "CreateMarkerMixin",
    "ParticleMarker",
    "FieldlineMarker",
    "GuidingcenterMarker",
    ]
