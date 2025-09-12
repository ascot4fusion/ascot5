"""This module contains all magnetic field input variants."""
import ctypes

from .cartesian import BfieldCartesian, CreateBfieldCartesianMixin
from .analytical import BfieldAnalytical, CreateBfieldAnalyticalMixin
from .axisymmetric import Bfield2D, CreateBfield2DMixin
from .perturbed import Bfield3D, CreateBfield3DMixin
from .stellarator import BfieldStellarator, CreateBfieldStellaratorMixin


def with_use_method(cls):
    # derive mapping from _fields_
    mapping = {
        f[0]: i for i, f in enumerate(cls._fields_)
        if f[0] != "type"  # exclude the type field
    }

    def use(self, variant):
        variant.stage()
        for name in mapping:
            if name == variant.variant:
                setattr(self, name, ctypes.pointer(variant._struct_))
                self.type = ctypes.c_uint32(mapping[name])
                break
        else:
            raise ValueError(f"Unknown variant {variant.variant}")

    cls.use = use
    return cls


@with_use_method
class Bfield(ctypes.Structure):
    """Wrapper for the magnetic field data in libascot.so."""

    _fields_ = [
        ('BTC', ctypes.POINTER(BfieldCartesian.Struct)),
        ('BGS', ctypes.POINTER(BfieldAnalytical.Struct)),
        ('B2DS', ctypes.POINTER(Bfield2D.Struct)),
        ('B3DS', ctypes.POINTER(Bfield3D.Struct)),
        ('BSTS', ctypes.POINTER(BfieldStellarator.Struct)),
        ('type', ctypes.c_uint32),
    ]


# pylint: disable=too-many-ancestors
class CreateBfieldMixin(
    CreateBfieldCartesianMixin,
    CreateBfieldAnalyticalMixin,
    CreateBfield2DMixin,
    CreateBfield3DMixin,
    CreateBfieldStellaratorMixin,
    ):
    """Mixin class used by `Data` to create magnetic field input.

    This class just combines all the magnetic field mixin classes.
    """

__all__  = [
    "CreateBfieldMixin",
    "BfieldCartesian",
    "BfieldAnalytical",
    "BfieldStellarator",
    "Bfield2D",
    "Bfield3D",
    ]
