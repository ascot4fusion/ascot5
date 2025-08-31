"""Defines EfieldCartesian electric field input class and the corresponding
factory method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class EfieldCartesian(InputVariant):
    """Electric field in Cartesian basis for testing purposes."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in E_TC.h."""
        _pack_ = 1
        _fields_ = [
            ('Exyz', ctypes.c_double * 3),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="EfieldCartesian",
            struct=EfieldCartesian.Struct(),
            )
        self._exyz: unyt.unyt_array

    @property
    def exyz(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._staged:
            return self._from_struct_("Exyz", shape=(3,), units="V/m")
        if self._format == Format.HDF5:
            return self._read_hdf5("exyz")
        return self._exyz.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "exyz":self.exyz,
        }
        return data

    def stage(self):
        init = LIBASCOT.E_TC_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.exyz,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._exyz
            self._staged = True

    def unstage(self):
        free = LIBASCOT.E_TC_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._exyz = self.exyz
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateEfieldCartesianMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create EfieldCartesian input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_efieldcartesian(
            self,
            exyz: utils.ArrayLike | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> EfieldCartesian:
        r"""Create a constant electric field input which is defined on
        a Cartesian basis.

        This method creates an electric field input which is simple and mainly
        used when the electric field is to be disabled in a simulation or when
        testing the orbit integration.

        To disable the electric field, simply call this method without setting
        any input data and the resulting field will be zero everywhere.

        Parameters
        ----------
        exyz : array_like (3,)
            Value of the electric field vector in the simulation domain.
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        inputdata : ~a5py.data.efield.EfieldCartesian
            Freshly minted input data object.

        Notes
        -----
        The electric field is simply a constant vector in Cartesian basis:

        .. math::
            \mathbf{E} = E_x \hat{x} + E_y \hat{y} + E_z \hat{z}.

        Note that this input cannot be used to set a constant electric field
        in cylindrical coordinates (unless the electric field is zero).
        """
        parameters = variants.parse_parameters(
            exyz,
        )
        variants.validate_required_parameters(
            parameters,
            names=["exyz",],
            units=["V/m",],
            shape=[(3,),],
            dtype="f8",
            default=[np.zeros((3,)),],
        )
        meta = variants.new_metadata("EfieldCartesian", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
