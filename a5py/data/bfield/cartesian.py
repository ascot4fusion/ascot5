"""Defines `BfieldCartesian` a Cartesian magnetic field input class and the
corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import _variants, InputVariant, Status, DataStruct
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class BfieldCartesian(InputVariant):
    """Magnetic field in Cartesian basis for testing purposes."""

    # pylint: disable=too-few-public-methods
    class Struct(DataStruct):
        """Python wrapper for the struct in B_TC.h."""
        _pack_ = 1
        _fields_ = [
            ('axisr', ctypes.c_double),
            ('axisz', ctypes.c_double),
            ('psival', ctypes.c_double),
            ('rhoval', ctypes.c_double),
            ('B', ctypes.c_double * 3),
            ('dB', ctypes.c_double * 9),
            ]

    @property
    def bxyz(self) -> unyt.unyt_array:
        """Magnetic field vector in Cartesian basis at origo."""
        if self.status is Status.SAVED:
            return self._file.read_hdf5("bxyz")
        elif self.status & Status.STAGED:
            return self._cdata.readonly_carray("B", (3,), "T")

    @property
    def jacobian(self) -> unyt.unyt_array:
        """Magnetic field Jacobian :math:`jac[i,j] = dB_i/dx_j`."""
        if self._staged:
            return self._from_struct_("dB", shape=(3,3), units="T/m")
        if self.status is Status.SAVED:
            return self._read_hdf5("jacobian")

    @property
    def rhoval(self) -> unyt.unyt_array:
        """Square root of the normalized poloidal flux which has a constant
        value."""
        if self._staged:
            return self._from_struct_("rhoval", shape=(1,), units="1")
        if self.status is Status.SAVED:
            return self._read_hdf5("rhoval")

    @property
    def psival(self) -> unyt.unyt_array:
        """Normalized poloidal flux which has a constant value."""
        if self._staged:
            return self._from_struct_("psival", shape=(1,), units="Wb/m")
        if self.status is Status.SAVED:
            return self._read_hdf5("psival")

    @property
    def axisrz(self) -> unyt.unyt_array:
        """Magnetic axis R and z coordinates."""
        if self._staged:
            return unyt.unyt_array((
                self._from_struct_("axisr", shape=(1,), units="m"),
                self._from_struct_("axisz", shape=(1,), units="m")
            )).T
        if self.status is Status.SAVED:
            return unyt.unyt_array((
                self._read_hdf5("axisr"), self._read_hdf5("axisz")
                ))

    def _stage(self, bxyz):
        cdata = BfieldCartesian.Struct()
        cdata.bxyz = bxyz
        self._cdata = cdata

    def _save_data(self):
        self._file.write("bxyz", self.bxyz)

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        axisrz = data["axisrz"]
        data["axisr"], data["axisz"] = axisrz[0], axisrz[1]
        del data["axisrz"]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "bxyz":self.bxyz,
            "jacobian":self.jacobian,
            "axisrz":self.axisrz,
            "rhoval":self.rhoval,
            "psival":self.psival,
        }
        return data

    def stage(self):
        init = LIBASCOT.B_TC_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.axisrz[0].v,
                self.axisrz[1].v,
                self.psival.v,
                self.rhoval.v,
                self.bxyz.v,
                self.jacobian.v,
                ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._bxyz
                del self._jacobian
            self._staged = True


    def unstage(self):
        free = LIBASCOT.B_TC_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._bxyz = self.bxyz
                self._jacobian = self.jacobian
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateBfieldCartesianMixin():
    """Mixin class used by `Data` to create :class:`BfieldCartesian` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_bfieldcartesian(
            self,
            bxyz: utils.ArrayLike | None = None,
            jacobian: utils.ArrayLike | None = None,
            axisrz: Tuple[float, float] | None = None,
            rhoval: float | None = None,
            psival: Optional[float] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> BfieldCartesian:
        r"""Create a magnetic field input which is defined on a Cartesian basis.

        This method creates a magnetic field input which is simple and mainly
        serves to validate the orbit-integrators.

        The input also defines the magnetic axis coordinates and the value of
        (normalized) poloidal flux, but these are meaningless in this case and
        are only provided because the code assumes the existence of some flux
        surfaces.

        Parameters
        ----------
        bxyz : array_like (3,)
            Magnetic field vector at origo in Cartesian basis.
        jacobian : array_like (3,3)
            Magnetic field Jacobian :math:`J[i,j] = \partial B_i/ \partial x_j`.
        axisrz : array_like (2,)
            Magnetic axis :math:`(R,z)` coordinates.

            In reality, there is no magnetic axis so these values can be
            almost anything.
        rhoval : float
            Square root of the normalized poloidal flux which has a constant
            value.
        psival : float, optional
            Normalized poloidal flux which has a constant value.

            Has the same value as `rhoval` by default.
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
        inputdata : ~a5py.data.bfield.BfieldCartesian
            Freshly minted input data object.

        Notes
        -----
        The magnetic field is given by

        .. math::

            \mathbf{B}(\mathbf{r}) = \mathbf{B}_0 + \mathbf{r}\cdot\mathbf{J},

        where both the vector,
        :math:`\mathbf{B}_0 = b_x \hat{x} + b_z \hat{y} + b_z \hat{z}`, and
        the Jacobian :math:`J=\partial B_i/\partial x_i` are given in Cartesian
        basis. The vector given as input is the value of the magnetic field on
        origo.
        """
        parameters = _variants.parse_parameters(
            bxyz, jacobian, axisrz, rhoval, psival,
        )
        _variants.validate_required_parameters(
            parameters,
            names=["bxyz", "jacobian", "axisrz", "rhoval"],
            units=["T", "T/m", "m", "1"],
            shape=[(3,), (3,3), (2,), ()],
            dtype="f8",
            default=[
                np.ones((3,)), np.zeros((3,3)), (1.0, 0.0), 0.5,
                ],
        )
        _variants.validate_optional_parameters(
            parameters,
            ["psival"], ["Wb/m"], [()], "f8", [parameters["rhoval"]],
        )
        leaf = BfieldCartesian(note=note)
        leaf._stage(bxyz=parameters["bxyz"])
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=store_hdf5, category="bfield",
            )
        return leaf
