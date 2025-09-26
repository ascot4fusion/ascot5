"""Defines MhdDynamic MHD eigenmode input class and the corresponding factory
method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, interp2D_data, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin



# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in mhdnonstat.h."""

    _fields_ = [
        ('n_modes', ctypes.c_int32),
        ('rho_min', ctypes.c_double),
        ('rho_max', ctypes.c_double),
        ('nmode', ctypes.POINTER(ctypes.c_int32)),
        ('mmode', ctypes.POINTER(ctypes.c_int32)),
        ('amplitude_nm', ctypes.POINTER(ctypes.c_double)),
        ('omega_nm', ctypes.POINTER(ctypes.c_double)),
        ('phase_nm', ctypes.POINTER(ctypes.c_double)),
        ('alpha_nm', ctypes.POINTER(interp2D_data)),
        ('phi_nm', ctypes.POINTER(interp2D_data)),
        ]


@Leaf.register
class MhdDynamic(InputVariant):
    """Electric field in Cartesian basis for testing purposes."""

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="MhdDynamic",
            struct=MhdDynamic.Struct(),
            )
        self._rhogrid: unyt.unyt_array
        self._timegrid: unyt.unyt_array
        self._toroidalnumber: unyt.unyt_array
        self._poloidalnumber: unyt.unyt_array
        self._magneticprofile: unyt.unyt_array
        self._electricprofile: unyt.unyt_array
        self._amplitude: unyt.unyt_array
        self._frequency: unyt.unyt_array
        self._phase: unyt.unyt_array

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.alpha_nm[0].x_min,
                self._struct_.alpha_nm[0].x_max,
                self._struct_.alpha_nm[0].n_x
                ) * unyt.dimensionless
        if self._format == Format.HDF5:
            nrho, rho0, rho1 = self._read_hdf5("nrho", "rhomin", "rhomax")
            return np.linspace(rho0, rho1, nrho)
        return self._rhogrid.copy()

    @property
    def timegrid(self) -> unyt.unyt_array:
        """Time grid in in which the data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.alpha_nm[0].y_min,
                self._struct_.alpha_nm[0].y_max,
                self._struct_.alpha_nm[0].n_y
                ) * unyt.s
        if self._format == Format.HDF5:
            ny, y0, y1 = self._read_hdf5("ntime", "timemin", "timemax")
            return np.linspace(y0, y1, ny)
        return self._timegrid.copy()

    @property
    def toroidalnumber(self):
        r"""Toroidal number :math:`n`."""
        if self._staged:
            nmode = self._from_struct_("n_modes", shape=())
            return self._from_struct_("nmode", shape=(nmode,))
        if self._format == Format.HDF5:
            return self._read_hdf5("toroidalnumber")
        return self._toroidalnumber

    @property
    def poloidalnumber(self):
        r"""Poloidal number :math:`m`."""
        if self._staged:
            nmode = self._from_struct_("n_modes", shape=())
            return self._from_struct_("mmode", shape=(nmode,))
        if self._format == Format.HDF5:
            return self._read_hdf5("poloidalnumber")
        return self._poloidalnumber

    @property
    def magneticprofile(self):
        r"""Magnetic eigenmode profile :math:`\alpha`."""
        if self._staged:
            nmode = self.toroidalnumber.size
            data = self._from_struct_("alpha_nm", idx=0)
            for i in range(1, nmode):
                data = np.stack(
                    (data, self._from_struct_("alpha_nm", idx=i)), axis=2
                    )
            if nmode == 1:
                data = np.expand_dims(data, axis=2)
            return data * unyt.m
        if self._format == Format.HDF5:
            return self._read_hdf5("magneticprofile")
        return self._magneticprofile.copy()

    @property
    def electricprofile(self):
        r"""Electric eigenmode profile :math:`\tilde{\Phi}`."""
        if self._staged:
            nmode = self.toroidalnumber.size
            data = self._from_struct_("phi_nm", idx=0)
            for i in range(1, nmode):
                data = np.stack(
                    (data, self._from_struct_("phi_nm", idx=i)), axis=2
                    )
            if nmode == 1:
                data = np.expand_dims(data, axis=2)
            return data * unyt.V
        if self._format == Format.HDF5:
            return self._read_hdf5("electricprofile")
        return self._electricprofile.copy()

    @property
    def amplitude(self):
        r"""Mode amplitude :math:`\lambda`."""
        if self._staged:
            nmode = self.toroidalnumber.size
            return self._from_struct_("amplitude_nm", shape=(nmode,))
        if self._format == Format.HDF5:
            return self._read_hdf5("amplitude")
        return self._amplitude

    @property
    def frequency(self):
        r"""Mode frequency :math:`\omega` [rad/s]."""
        if self._staged:
            nmode = self.toroidalnumber.size
            return self._from_struct_("omega_nm", shape=(nmode,))
        if self._format == Format.HDF5:
            return self._read_hdf5("frequency")
        return self._frequency

    @property
    def phase(self):
        r"""Mode phase :math:`\varphi` [rad]."""
        if self._staged:
            nmode = self.toroidalnumber.size
            return self._from_struct_("phase_nm", shape=(nmode,))
        if self._format == Format.HDF5:
            return self._read_hdf5("phase")
        return self._phase

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        for grid in ["rhogrid", "timegrid"]:
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
            "rhogrid":self.rhogrid,
            "timegrid":self.timegrid,
            "toroidalnumber":self.toroidalnumber,
            "poloidalnumber":self.poloidalnumber,
            "magneticprofile":self.magneticprofile,
            "electricprofile":self.electricprofile,
            "amplitude":self.amplitude,
            "frequency":self.frequency,
            "phase":self.phase,
        }
        return data

    def stage(self):
        init = LIBASCOT.mhd_nonstat_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ndpointer(ctypes.c_int32),
            ndpointer(ctypes.c_int32),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.toroidalnumber.size,
                self.rhogrid.size,
                self.timegrid.size,
                self.rhogrid[0].v,
                self.rhogrid[-1].v,
                self.timegrid[0].v,
                self.timegrid[-1].v,
                self.toroidalnumber,
                self.poloidalnumber,
                self.amplitude.v,
                self.frequency.v,
                self.phase.v,
                self.magneticprofile.v,
                self.electricprofile.v,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._magneticprofile
                del self._electricprofile
            self._staged = True

    def unstage(self):
        free = LIBASCOT.mhd_nonstat_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._magneticprofile = self.magneticprofile
                self._electricprofile = self.electricprofile
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create MhdDynamic input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_mhddynamic(
            self,
            rhogrid: utils.ArrayLike | None = None,
            timegrid: utils.ArrayLike | None = None,
            toroidalnumber: utils.ArrayLike | None = None,
            poloidalnumber: utils.ArrayLike | None = None,
            magneticprofile: utils.ArrayLike | None = None,
            electricprofile: utils.ArrayLike | None = None,
            amplitude: Optional[utils.ArrayLike]=None,
            frequency: Optional[utils.ArrayLike]=None,
            phase: Optional[utils.ArrayLike]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> MhdDynamic:
        r"""Create MHD eigenmode input where the modes evolve in time.

        This input is otherwise equivalent to :class:`MhdStatic` except that the
        profiles evolve in time. This slows simulation significantly if the number of modes is large.

        Parameters
        ----------
        rhogrid : array_like (nrho,)
            Radial grid in rho in which the data is tabulated.
        timegrid : array_like (ntime,)
            Time grid in which the data is tabulated.
        toroidalnumber : array_like (nmode,)
            Toroidal mode number :math:`n`.
        poloidalnumber : array_like (nmode,)
            Poloidal mode number :math:`m`.
        magneticprofile : array_like (nrho,ntime,nmode)
            Magnetic eigenmode profile :math:`\alpha`.
        electricprofile : array_like (nrho,ntime,nmode)
            Electric eigenmode profile :math:`\tilde{\Phi}`.
        amplitude : array_like (nmode,), optional
            Mode amplitude :math:`\lambda`.

            This is a scalar value that scales both :math:`\alpha` and
            :math:`\tilde{\Phi}`. This is equal to one by default.
        frequency : array_like (nmode,), optional
            Mode frequency :math:`\omega` [rad/s].

            By default this is set to zero.
        phase : array_like (nmode,), optional
            Mode phase :math:`\varphi` [rad].

            By default this is set to zero.
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        preview : bool, *optional*
            If True, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : ~a5py.data.mhd.MhdDynamic
            Input variant created from the given parameters.
        """
        parameters = _variants.parse_parameters(
            rhogrid, timegrid, toroidalnumber, poloidalnumber, magneticprofile,
            electricprofile, amplitude, frequency, phase,
        )
        default_rhogrid = np.linspace(0., 1., 3)
        default_timegrid = np.linspace(0., 1., 4)
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        ntime = (default_timegrid.size if parameters["timegrid"] is None
              else parameters["timegrid"].size)
        nmode = (2 if parameters["toroidalnumber"] is None
              else parameters["toroidalnumber"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rhogrid", "timegrid", "toroidalnumber", "poloidalnumber",
                   "magneticprofile", "electricprofile",],
            units=["1", "s", "1", "1", "m", "V",],
            shape=[(nrho,), (ntime,), (nmode,), (nmode,), (nrho,ntime,nmode),
                   (nrho,ntime, nmode),],
            dtype=["f8", "f8", "i4", "i4", "f8", "f8",],
            default=[default_rhogrid, default_timegrid, np.array([1, 2]),
                     np.array([3,3]), np.ones((nrho,ntime,nmode)),
                     np.ones((nrho,ntime,nmode)),],
        )
        _variants.validate_optional_parameters(
            parameters,
            names=["amplitude", "frequency", "phase",],
            units=["1", "rad/s", "rad",],
            shape=[(nmode,), (nmode,), (nmode,),],
            dtype=["f8", "f8", "f8",],
            default=[np.ones(nmode), np.zeros(nmode), np.zeros(nmode),],
        )
        meta = _variants.new_metadata("MhdDynamic", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
