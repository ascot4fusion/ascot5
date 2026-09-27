"""Defines 2D potential electric field input class and the corresponding factory
method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, Spline2D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in E_1DS.h."""

    _fields_ = [
        ("potential", Spline2D),
        ]

init_fun(
    "EfieldPotential2D_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
)

init_fun("EfieldPotential2D_free", ctypes.POINTER(Struct))

@Leaf.register
class EfieldPotential2D(InputVariant):
    """Axisymmetric electric field evaluated from 2D potential."""

    @property
    def rgrid(self) -> unyt.unyt_array:
        """Radial grid in :math:`R` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "m", "potential")
        assert self._file is not None
        return self._file.read("rgrid")

    @property
    def zgrid(self) -> unyt.unyt_array:
        """Axial grid in :math:`z` in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "m", "potential")
        assert self._file is not None
        return self._file.read("zgrid")

    @property
    def potential(self) -> unyt.unyt_array:
        """Electric field potential as a function of :math:`r` and :math:`z`."""
        if self._cdata is not None:
            return self._cdata.readonly_interp("potential", "V")
        assert self._file is not None
        return self._file.read("potential")


    def _stage(
        self, rgrid: unyt.unyt_array, zgrid: unyt.unyt_array, potential: unyt.unyt_array,
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.EfieldPotential2D_init(
            ctypes.byref(self._cdata),
            rgrid.size,
            zgrid.size,
            rgrid[[0, -1]].v,
            zgrid[[0, -1]].v,
            potential.v,
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("rgrid", self.rgrid)
        self._file.write("zgrid", self.zgrid)
        self._file.write("potential", self.potential)

    def export(self) -> dict[str, unyt.unyt_array]:
        data = {
            "rgrid": self.rgrid,
            "zgrid": self.zgrid,
            "potential": self.potential,
        }
        return data

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.EfieldPotential2D_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_efieldpotential2d(
            self,
            rgrid: utils.ArrayLike,
            zgrid: utils.ArrayLike,
            potential: utils.ArrayLike,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> EfieldPotential2D:
        r"""Create radial electric field input that is evaluated from the
        gradient of a 1D potential.

        This input was designed to use NEOTRANSP output.

        Parameters
        ----------
        rgrid : array_like (nr,)
            Radial grid in :math:`R` in which the data is tabulated.
        zgrid : array_like (nz,)
            Axial grid in :math:`z` in which the data is tabulated.
        potential : array_like (nr, nz)
            Electric field potential as a function of :math:`r` and :math:`z`.
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
        The electric field is evaluated from the gradient of the 2D potential:

        .. math::

            \mathbf{E} = \frac{\partial V}{\partial \r} \hat{\mathbf{r}}
                       + \frac{\partial V}{\partial z} \hat{\mathbf{z}}.
        """
        with utils.validate_variables() as v:
            rgrid = v.validate("rgrid", rgrid, (-1,), "m")
            zgrid = v.validate("zgrid", zgrid, (-1,), "m")

        nr, nz = rgrid.size, zgrid.size
        with utils.validate_variables() as v:
            potential = v.validate("potential", potential, (nr,nz), "V")

        utils.validate_abscissa(rgrid, "rgrid")
        utils.validate_abscissa(zgrid, "zgrid")
        leaf = EfieldPotential2D(note=note)
        leaf._stage(
            rgrid=rgrid, zgrid=zgrid, potential=potential,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="efield",
            )
        return leaf
