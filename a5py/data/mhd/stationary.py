"""Defines `MhdStationary` stationary MHD eigenmode input class and the
corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ..cstructs import interp1D_data
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class MhdStationary(InputVariant):
    """Stationary MHD eigenmode input."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in mhd_stat.h."""
        _pack_ = 1
        _fields_ = [
            ('n_modes', ctypes.c_int32),
            ('PADDING_0', ctypes.c_ubyte * 4),
            ('rho_min', ctypes.c_double),
            ('rho_max', ctypes.c_double),
            ('nmode', ctypes.POINTER(ctypes.c_int32)),
            ('mmode', ctypes.POINTER(ctypes.c_int32)),
            ('amplitude_nm', ctypes.POINTER(ctypes.c_double)),
            ('omega_nm', ctypes.POINTER(ctypes.c_double)),
            ('phase_nm', ctypes.POINTER(ctypes.c_double)),
            ('alpha_nm', ctypes.POINTER(interp1D_data)),
            ('phi_nm', ctypes.POINTER(interp1D_data)),
        ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="MhdStationary",
            struct=MhdStationary.Struct(),
            )
        self._rhogrid: unyt.unyt_array
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
                    (data, self._from_struct_("alpha_nm", idx=i)), axis=1
                    )
            if nmode == 1:
                data = np.expand_dims(data, axis=1)
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
                    (data, self._from_struct_("phi_nm", idx=i)), axis=1
                    )
            if nmode == 1:
                data = np.expand_dims(data, axis=1)
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
            "rhogrid":self.rhogrid,
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
        init = LIBASCOT.mhd_stat_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_int32,
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
                self.rhogrid[0].v,
                self.rhogrid[-1].v,
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
        free = LIBASCOT.mhd_stat_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._magneticprofile = self.magneticprofile
                self._electricprofile = self.electricprofile
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateMhdStationaryMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create `MhdStationary` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_mhdstationary(
            self,
            rhogrid: utils.ArrayLike | None = None,
            toroidalnumber: utils.ArrayLike | None = None,
            poloidalnumber: utils.ArrayLike | None = None,
            magneticprofile: utils.ArrayLike | None = None,
            electricprofile: utils.ArrayLike | None = None,
            amplitude: Optional[utils.ArrayLike] = None,
            frequency: Optional[utils.ArrayLike] = None,
            phase: Optional[utils.ArrayLike] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> MhdStationary:
        r"""Create MHD eigenmode input where the modes don't evolve in time.

        This input assumes that :math:`\alpha` and :math:`\tilde{\Phi}` have
        only radial dependence (but the modes still rotate if :math:`\omega`
        is non-zero). These are much faster to evaluate than the time-dependent
        modes if the number of modes is large.

        Parameters
        ----------
        rhogrid : array_like (nrho,)
            Radial grid in rho in which the data is tabulated.
        toroidalnumber : array_like (nmode,)
            Toroidal mode number :math:`n`.
        poloidalnumber : array_like (nmode,)
            Poloidal mode number :math:`m`.
        magneticprofile : array_like (nrho,nmode)
            Magnetic eigenmode profile :math:`\alpha`.
        electricprofile : array_like (nrho,nmode)
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
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        inputdata : ~a5py.data.mhd.MhdStationary
            Freshly minted input data object.

        Notes
        -----
        This input can be used to include any perturbations of type

        .. math::

            \alpha       &= \sum_{nm} \lambda_{nm} \alpha_{nm}
                \cos\left(n\zeta-m\theta-\omega_{nm}t + \varphi\right),\\
            \tilde{\Phi} &= \sum_{nm} \lambda_{nm} \Phi_{nm}
                \cos\left(n\zeta-m\theta-\omega_{nm}t + \varphi\right),

        where :math:`\zeta` and :math:`\theta` are the toroidal and poloidal
        Boozer coordinates, to the EM-field as

        .. math::

            \mathbf{B} &= \mathbf{B}_\mathrm{bkg}
                + \nabla\cross\alpha\mathbf{\mathbf{B}},\\
            \mathbf{E} &= \mathbf{E}_\mathrm{bkg} - \nabla \tilde{\Phi}
                - \frac{\partial \alpha \mathbf{E}}{\partial t},

        where :math:`\mathbf{B}_\mathrm{bkg}` and
        :math:`\mathbf{E}_\mathrm{bkg}` are the background magnetic and electric
        fields.

        For rapidly rotating modes, the electrons are able to balance any
        electric field parallel to the field lines. The condition
        :math:`E_\parallel=0` makes the magnetic and electric perturbations
        co-dependent as:

        .. math::

            \omega_{nm}\alpha_{nm} = \frac{nq - m}{I+gq}\Phi_{nm},

        where :math:`q` is safety factor, :math:`g=RB_\mathrm{phi}`, and
        :math:`I(\psi)` is toroidal current function (see
        ~a5py.data.boozer.Boozer for details). However, this module *does not
        enforce* this condition automatically and it is up to the user to ensure
        this whenever appropriate.
        """
        parameters = variants.parse_parameters(
            rhogrid, toroidalnumber, poloidalnumber, magneticprofile,
            electricprofile, amplitude, frequency, phase,
        )
        default_rhogrid = np.linspace(0., 1., 3)
        nrho = (default_rhogrid.size if parameters["rhogrid"] is None
              else parameters["rhogrid"].size)
        nmode = (2 if parameters["toroidalnumber"] is None
              else parameters["toroidalnumber"].size)
        variants.validate_required_parameters(
            parameters,
            names=["rhogrid", "toroidalnumber", "poloidalnumber",
                   "magneticprofile", "electricprofile",],
            units=["1", "1", "1", "m", "V",],
            shape=[(nrho,), (nmode,), (nmode,), (nrho,nmode), (nrho,nmode),],
            dtype=["f8", "i4", "i4", "f8", "f8",],
            default=[default_rhogrid, np.array([1, 2]), np.array([3,3]),
                     np.ones((nrho,nmode)), np.ones((nrho,nmode)),],
        )
        variants.validate_optional_parameters(
            parameters,
            names=["amplitude", "frequency", "phase",],
            units=["1", "rad/s", "rad",],
            shape=[(nmode,), (nmode,), (nmode,),],
            dtype=["f8", "f8", "f8",],
            default=[np.ones(nmode), np.zeros(nmode), np.zeros(nmode),],
        )
        meta = variants.new_metadata("MhdStationary", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
