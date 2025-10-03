"""Defines the Cartesian electric field input class and the corresponding
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
    """Python wrapper for the struct in E_TC.h."""

    _fields_ = [
        ("exyz", ctypes.c_double * 3),
    ]


init_fun(
    "EfieldCartesian_init",
    ctypes.POINTER(Struct),
    ndpointer(ctypes.c_double),
)

init_fun("EfieldCartesian_free", ctypes.POINTER(Struct))


@Leaf.register
class EfieldCartesian(InputVariant):
    """Electric field in Cartesian basis for testing purposes."""

    @property
    def exyz(self) -> unyt.unyt_array:
        """Electric field vector (uniform in time and space)."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("exyz", (3,), "V/m")
        assert self._file is not None
        return self._file.read("exyz")

    # pylint: disable=too-many-arguments
    def _stage(self, exyz: unyt.unyt_array) -> None:
        self._cdata = Struct()
        if LIBASCOT.EfieldCartesian_init(ctypes.byref(self._cdata), exyz.v):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("exyz", self.exyz)

    def export(self) -> dict[str, unyt.unyt_array]:
        return {"exyz": self.exyz}

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.EfieldCartesian_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    # pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_efieldcartesian(
        self,
        exyz: utils.ArrayLike,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
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
        inputdata : :class:`.EfieldCartesian`
            Input variant created from the given parameters.

        Notes
        -----
        The electric field is simply a constant vector in Cartesian basis:

        .. math::
            \mathbf{E} = E_x \hat{x} + E_y \hat{y} + E_z \hat{z}.

        Note that this input cannot be used to set a constant electric field
        in cylindrical coordinates (unless the electric field is zero).
        """
        with utils.validate_variables() as v:
            exyz = v.validate("exyz", exyz, (3,), "V/m")

        leaf = EfieldCartesian(note=note)
        leaf._stage(exyz=exyz)
        if preview:
            return leaf

        self._treemanager.enter_leaf(
            leaf,
            activate=activate,
            save=save,
            category="efield",
        )
        return leaf
