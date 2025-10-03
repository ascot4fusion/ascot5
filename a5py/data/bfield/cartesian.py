"""Defines the Cartesian magnetic field input class and the corresponding
factory method.
"""
import ctypes
from typing import Optional

import unyt
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin

# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in B_TC.h."""

    _fields_ = [
        ("axisrz", ctypes.c_double * 2),
        ("psival", ctypes.c_double),
        ("rhoval", ctypes.c_double),
        ("bxyz", ctypes.c_double * 3),
        ("jacobian", ctypes.c_double * 9),
        ]

init_fun(
    "BfieldCartesian_init",
    ctypes.POINTER(Struct),
    *(2*[ctypes.c_double]),
    *(3*[ndpointer(ctypes.c_double)]),
    restype=ctypes.c_int32,
    )

init_fun("BfieldCartesian_free", ctypes.POINTER(Struct))


@Leaf.register
class BfieldCartesian(InputVariant):
    """Magnetic field in Cartesian basis for testing purposes."""

    _cdata: Optional[Struct]

    @property
    def bxyz(self) -> unyt.unyt_array:
        """Magnetic field vector in Cartesian basis at origo."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("bxyz", (3,), "T")
        assert self._file is not None
        return self._file.read("bxyz")

    @property
    def jacobian(self) -> unyt.unyt_array:
        r"""Magnetic field Jacobian :math:`jac[i,j] = dB_i/dx_j`."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("jacobian", (3,3), "T/m")
        assert self._file is not None
        return self._file.read("jacobian")

    @property
    def rhoval(self) -> unyt.unyt_array:
        """Square root of the normalized poloidal flux which has a constant
        value.
        """
        if self._cdata is not None:
            return self._cdata.readonly_carray("rhoval", (), "1")
        assert self._file is not None
        return self._file.read("rhoval")

    @property
    def psival(self) -> unyt.unyt_array:
        """Normalized poloidal flux which has a constant value."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("psival", (), "Wb/rad")
        assert self._file is not None
        return self._file.read("psival")

    @property
    def axisrz(self) -> unyt.unyt_array:
        r"""Magnetic axis :math:`R` and :math:`z` coordinates."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("axisrz", (2,), "m")
        assert self._file is not None
        return self._file.read("axisrz")

    # pylint: disable=too-many-arguments
    def _stage(self, axisrz: unyt.unyt_array,
               psival: unyt.unyt_array,
               rhoval: unyt.unyt_array,
               bxyz: unyt.unyt_array,
               jacobian: unyt.unyt_array,
               ) -> None:
        self._cdata = Struct()
        if LIBASCOT.BfieldCartesian_init(
            ctypes.byref(self._cdata),
            psival.v,
            rhoval.v,
            axisrz.v,
            bxyz.v,
            jacobian.v,
            ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for field in ["bxyz", "jacobian", "axisrz", "rhoval", "psival"]:
            self._file.write(field, getattr(self, field))

    def export(self) -> dict[str, unyt.unyt_array]:
        fields = ["bxyz", "jacobian", "axisrz", "rhoval", "psival"]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.BfieldCartesian_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_bfieldcartesian(
            self,
            bxyz: utils.ArrayLike,
            jacobian: utils.ArrayLike,
            axisrz: utils.ArrayLike,
            rhoval: utils.ArrayLike,
            psival: Optional[utils.ArrayLike]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
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
        psival : float, *optional*
            Normalized poloidal flux which has a constant value.

            Has the same value as `rhoval` by default.
        note : str, *optional*
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, *optional*
            Set this input as active on creation.
        preview : bool, *optional*
            If ``True``, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        data : :class:`.BfieldCartesian`
            Input variant created from the given parameters.

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
        with utils.validate_variables() as v:
            bxyz = v.validate("bxyz", bxyz, (3,), "T")
            jacobian = v.validate("jacobian", jacobian, (3,3), "T/m")
            axisrz = v.validate("axisrz", axisrz, (2,), "m")
            rhoval = v.validate("rhoval", rhoval, (), "1")

        with utils.validate_variables() as v:
            psival = v.validate("psival", psival, (), "Wb/rad", default=rhoval.v)

        leaf = BfieldCartesian(note=note)
        leaf._stage(axisrz, psival, rhoval, bxyz, jacobian)
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="bfield",
            )
        return leaf
