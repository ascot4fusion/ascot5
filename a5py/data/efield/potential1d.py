"""Defines 1D potential electric field input class and the corresponding factory
method.
"""

import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, Linear1D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in E_1DS.h."""

    _fields_ = [
        ("dvdrho", Linear1D),
    ]


init_fun(
    "EfieldPotential1D_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
)

init_fun("EfieldPotential1D_free", ctypes.POINTER(Struct))


@Leaf.register
class EfieldPotential1D(InputVariant):
    """Radial electric field evaluated from the gradient of a 1D potential."""

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "1", "dvdrho")
        assert self._file is not None
        return self._file.read("rhogrid")

    @property
    def dvdrho(self) -> unyt.unyt_array:
        """Derivative of the electric potential with respect to minor radius."""
        if self._cdata is not None:
            return self._cdata.readonly_interp("dvdrho", "V")
        assert self._file is not None
        return self._file.read("dvdrho")


    def _stage(
        self, dvdrho: unyt.unyt_array, rhogrid: unyt.unyt_array
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.EfieldPotential1D_init(
            ctypes.byref(self._cdata),
            rhogrid.size,
            rhogrid[[0, -1]],
            dvdrho,
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("exyz", self.exyz)

    def export(self) -> dict[str, unyt.unyt_array]:
        data = {
            "rhogrid": self.rhogrid,
            "dvdrho": self.dvdrho,
        }
        return data

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.EfieldPotential1D_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    # pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_efieldpotential1d(
        self,
        rhogrid: utils.ArrayLike,
        dvdrho: utils.ArrayLike,
        note: Optional[str]=None,
        activate: bool=False,
        preview: bool=False,
        save: Optional[bool]=None,
    ) -> EfieldPotential1D:
        r"""Create radial electric field input that is evaluated from the
        gradient of a 1D potential.

        This input was designed to use NEOTRANSP output.

        Parameters
        ----------
        rhogrid : array_like (nrho,)
            Radial grid in rho in which the data is tabulated.
        dvdrho : array_like (nrho,)
            Derivative of electric potential with respect to minor radius.

            If :math:`r_\mathrm{eff} = 1` m, this is essentially equal to
            :math:`\partial V/ \partial r`.
        note : str, *optional*
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, *optional*
            Set this input as active on creation.
        preview : bool, *optional*
            If True, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : ~a5py.data.efield.EfieldRadialPotential
            Input variant created from the given parameters.

        Notes
        -----
        The electric field is evaluated from the gradient of the 1D potential
        and the gradient of the square of the normalized poloidal flux:

        .. math::

            \mathbf{E} = \frac{\partial V}{\partial \rho}
                         \nabla \rho.
        """
        with utils.validate_variables() as v:
            rhogrid = v.validate("rhogrid", rhogrid, (-1,), "1")

        nrho = rhogrid.size
        with utils.validate_variables() as v:
            dvdrho = v.validate("dvdrho", dvdrho, (nrho,), "V")

        utils.validate_abscissa(rhogrid, "rhogrid")
        leaf = EfieldPotential1D(note=note)
        leaf._stage(
            rhogrid=rhogrid, dvdrho=dvdrho,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="efield",
            )
        return leaf
