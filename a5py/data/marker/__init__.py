"""This module contains all marker input variants and combines their
factory methods to a single class."""
from .particle import ParticleMarker, CreateParticleMixin
from .guidingcenter import GuidingCenterMarker, CreateGuidingcenterMixin
from .fieldline import FieldlineMarker, CreateFieldlineMixin

# pylint: disable=too-many-ancestors
class CreateMarkerMixin(
    CreateParticleMixin, CreateGuidingcenterMixin, CreateFieldlineMixin,
    ):
    """Mixin class used by `Data` to create marker input.

    This class just combines all the marker mixin classes.
    """

__all__  = [
    "CreateMarkerMixin", "ParticleMarker", "GuidingCenterMarker",
    "FieldlineMarker",
    ]
