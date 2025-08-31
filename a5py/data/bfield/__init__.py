"""This module contains all magnetic field input variants."""
import ctypes

from .cartesian import BfieldCartesian, CreateBfieldCartesianMixin
from .analytical import BfieldAnalytical, CreateBfieldAnalyticalMixin
from .axisymmetric import Bfield2D, CreateBfield2DMixin
from .perturbed import Bfield3D, CreateBfield3DMixin
from .stellarator import BfieldStellarator, CreateBfieldStellaratorMixin


class Bfield(ctypes.Structure):
    """Wrapper for the magnetic field data in libascot.so."""

    _pack_ = 1
    _fields_ = [
        ('BTC', ctypes.POINTER(BfieldCartesian.Struct)),
        ('BGS', ctypes.POINTER(BfieldAnalytical.Struct)),
        ('B2DS', ctypes.POINTER(Bfield2D.Struct)),
        ('B3DS', ctypes.POINTER(Bfield3D.Struct)),
        ('BSTS', ctypes.POINTER(BfieldStellarator.Struct)),
        ('type', ctypes.c_uint32),
    ]

    def use(self, variant, mpirank, mpisize):
        """Initialize the pointer and set the type corresponding to the data."""
        variant_vs_name_and_enum = {
            "BfieldCartesian": ("BTC", 0),
            "BfieldAnalytical": ("BGS", 1),
            "Bfield2D": ("B2DS", 2),
            "Bfield3D": ("B3DS", 3),
            "BfieldStellarator": ("BSTS", 4),
            }
        if mpirank == 0:
            variant.stage()
            variant_name = variant.variant
            variant_struct = variant._struct_
        else:
            variant_name = None
            variant_struct = None
        name, enum = variant_vs_name_and_enum[variant_name]
        setattr(self, name, ctypes.pointer(variant_struct))
        self.type = ctypes.c_uint32(enum)


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
