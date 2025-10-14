"""Defines 1D potential electric field input class and the corresponding factory
method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, Spline1D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in E_1DS.h."""

    _fields_ = [
        ('reff', ctypes.c_double),
        ('dv', Spline1D),
        ]


@Leaf.register
class EfieldPotential1D(InputVariant):
    """Radial electric field evaluated from the gradient of a 1D potential."""

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.dv.x_min,
                self._struct_.dv.x_max,
                self._struct_.dv.n_x
                ) * unyt.dimensionless
        if self._format == Format.HDF5:
            nx, x0, x1 = self._read_hdf5("nrho", "rhomin", "rhomax")
            return np.linspace(x0, x1, nx) * unyt.dimensionless

    @property
    def dvdrho(self) -> unyt.unyt_array:
        """Derivative of the electric potential with respect to minor radius."""
        if self._staged:
            return self._from_struct_("dV", units="V/m")
        if self._format == Format.HDF5:
            return self._read_hdf5("dvdrho")

    @property
    def reff(self) -> unyt.unyt_array:
        """Effective minor radius."""
        if self._staged:
            return self._from_struct_("reff", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("reff")

    def export(self):
        data = {
            "reff":self.reff,
            "rhogrid":self.rhogrid,
            "dvdrho":self.dvdrho,
        }
        return data

    def stage(self):
        init = LIBASCOT.EfieldPotential1D_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.rhogrid.size,
                self.rhogrid[0].v,
                self.rhogrid[-1].v,
                self.reff[0].v,
                self.dvdrho,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._dvdrho
            self._staged = True

    def unstage(self):
        free = LIBASCOT.EfieldPotential1D_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._dvdrho = self.dvdrho
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_efieldpotential1d(
            self,
            rhogrid: utils.ArrayLike,
            dvdrho: utils.ArrayLike,
            reff: Optional[float]=None,
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
        reff : float, *optional*
            Effective minor radius.

            This is defined as
            :math:`r_\mathrm{eff} = \partial r/ \partial \rho`, and it is used
            to convert :math:`\partial V/ \partial \rho` to
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
        parameters = _variants.parse_parameters(
            rhogrid, dvdrho, reff,
        )
        default_rhogrid = np.linspace(0., 1., 3)
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rhogrid", "dvdrho", "reff",],
            units=["1", "V/m", "m",],
            shape=[(nrho,), (nrho,), (1,)],
            dtype="f8",
            default=[np.array([0., 0.5, 1.]), np.zeros((3,)), 1.],
        )
        meta = _variants.new_metadata("EfieldRadialPotential", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
