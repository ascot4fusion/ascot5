"""Defines :class:`EfieldRadialPotential` electric field input class and the
corresponding factory method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ..cstructs import interp1D_data
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class EfieldRadialPotential(InputVariant):
    """Radial electric field evaluated from the gradient of a 1D potential."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in E_1DS.h."""
        _pack_ = 1
        _fields_ = [
            ('reff', ctypes.c_double),
            ('dV', interp1D_data),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="EfieldRadialPotential",
            struct=EfieldRadialPotential.Struct(),
            )
        self._reff: unyt.unyt_array
        self._rhogrid: unyt.unyt_array
        self._dvdrho: unyt.unyt_array

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.dV.x_min,
                self._struct_.dV.x_max,
                self._struct_.dV.n_x
                ) * unyt.dimensionless
        if self._format == Format.HDF5:
            nx, x0, x1 = self._read_hdf5("nrho", "rhomin", "rhomax")
            return np.linspace(x0, x1, nx) * unyt.dimensionless
        return self._rhogrid.copy()

    @property
    def dvdrho(self) -> unyt.unyt_array:
        """Derivative of the electric potential with respect to minor radius."""
        if self._staged:
            return self._from_struct_("dV", units="V/m")
        if self._format == Format.HDF5:
            return self._read_hdf5("dvdrho")
        return self._dvdrho.copy()

    @property
    def reff(self) -> unyt.unyt_array:
        """Effective minor radius."""
        if self._staged:
            return self._from_struct_("reff", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("reff")
        return self._reff.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        for grid in ["rhogrid"]:
            name = grid.replace("grid", "")
            data["n" + name] = data[grid].size
            data[name + "min"] = data[grid][0]
            data[name + "max"] = data[grid][-1]
            del data[grid]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "reff":self.reff,
            "rhogrid":self.rhogrid,
            "dvdrho":self.dvdrho,
        }
        return data

    def stage(self):
        init = LIBASCOT.E_1DS_init
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
        free = LIBASCOT.E_1DS_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._dvdrho = self.dvdrho
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateEfieldRadialPotentialMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create EfieldRadialPotential input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_efieldradialpotential(
            self,
            rhogrid: utils.ArrayLike | None = None,
            dvdrho: utils.ArrayLike | None = None,
            reff: float | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> EfieldRadialPotential:
        r"""Create radial electric field input that is evaluated from the
        gradient of a 1D potential.

        Parameters
        ----------
        rhogrid : array_like (nrho,)
            Radial grid in rho in which the data is tabulated.
        dvdrho : array_like (nrho,)
            Derivative of electric potential with respect to minor radius.

            If :math:`r_\mathrm{eff} = 1` m, this is essentially equal to
            :math:`\partial V/ \partial r`.
        reff : float
            Effective minor radius.

            This is defined as
            :math:`r_\mathrm{eff} = \partial r/ \partial \rho`, and it is used
            to convert :math:`\partial V/ \partial \rho` to
            :math:`\partial V/ \partial r`.
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
        inputdata : ~a5py.data.efield.EfieldRadialPotential
            Freshly minted input data object.

        Notes
        -----
        The electric field is evaluated from the gradient of the 1D potential
        and the gradient of the square of the normalized poloidal flux:

        .. math::

            \mathbf{E} = \frac{\partial V}{\partial \rho}
                         \nabla \rho.
        """
        parameters = variants.parse_parameters(
            rhogrid, dvdrho, reff,
        )
        default_rhogrid = np.linspace(0., 1., 3)
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        variants.validate_required_parameters(
            parameters,
            names=["rhogrid", "dvdrho", "reff",],
            units=["1", "V/m", "m",],
            shape=[(nrho,), (nrho,), (1,)],
            dtype="f8",
            default=[np.array([0., 0.5, 1.]), np.zeros((3,)), 1.],
        )
        meta = variants.new_metadata("EfieldRadialPotential", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
