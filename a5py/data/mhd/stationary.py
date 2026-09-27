"""Defines `MhdStationary` stationary MHD eigenmode input class and the
corresponding factory method.
"""

import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, Spline1D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in mhd_stat.h."""

    _fields_ = [
        ("n", ctypes.c_size_t),
        ("nmode", ctypes.POINTER(ctypes.c_int32)),
        ("mmode", ctypes.POINTER(ctypes.c_int32)),
        ("amplitude", ctypes.POINTER(ctypes.c_double)),
        ("omega", ctypes.POINTER(ctypes.c_double)),
        ("phase", ctypes.POINTER(ctypes.c_double)),
        ("alpha", ctypes.POINTER(Spline1D)),
        ("phi", ctypes.POINTER(Spline1D)),
    ]


init_fun(
    "MhdStationary_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
    ctypes.c_size_t,
    ndpointer(ctypes.c_int32),
    ndpointer(ctypes.c_int32),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
)

init_fun("MhdStationary_free", ctypes.POINTER(Struct))


@Leaf.register
class MhdStationary(InputVariant):
    """Stationary MHD eigenmode input."""

    @property
    def rhogrid(self) -> unyt.unyt_array:
        """Radial grid in rho in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "1", "alpha", 0)
        assert self._file is not None
        return self._file.read("rhogrid")

    @property
    def number_of_modes(self):
        r"""Number of eigenmodes."""
        if self._cdata is not None:
            return self._cdata.n
        assert self._file is not None
        return self._file.read("toroidalnumber").size

    @property
    def toroidalnumber(self):
        r"""Toroidal number :math:`n`."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("nmode", (self.number_of_modes,), "1")
        assert self._file is not None
        return self._file.read("toroidalnumber")

    @property
    def poloidalnumber(self):
        r"""Poloidal number :math:`m`."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("mmode", (self.number_of_modes,), "1")
        assert self._file is not None
        return self._file.read("poloidalnumber")

    @property
    def magneticprofile(self):
        r"""Magnetic eigenmode profile :math:`\alpha`."""
        if self._cdata is not None:
            nmode = self.toroidalnumber.size
            data = self._cdata.readonly_interp("alpha", "m", idx=0)
            for i in range(1, nmode):
                data = np.stack((data, self._cdata.readonly_interp("alpha", "m", idx=i)), axis=1)
            if nmode == 1:
                data = np.expand_dims(data, axis=1)
            return data
        assert self._file is not None
        return self._file.read("magneticprofile")

    @property
    def electricprofile(self):
        r"""Electric eigenmode profile :math:`\tilde{\Phi}`."""
        if self._cdata is not None:
            nmode = self.toroidalnumber.size
            data = self._cdata.readonly_interp("phi", "V", idx=0)
            for i in range(1, nmode):
                data = np.stack((data, self._cdata.readonly_interp("phi", "V", idx=i)), axis=1)
            if nmode == 1:
                data = np.expand_dims(data, axis=1)
            return data
        assert self._file is not None
        return self._file.read("electricprofile")

    @property
    def amplitude(self):
        r"""Mode amplitude :math:`\lambda`."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("amplitude", (self.number_of_modes,))
        assert self._file is not None
        return self._file.read("amplitude")

    @property
    def frequency(self):
        r"""Mode frequency :math:`\omega` [rad/s]."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("omega", (self.number_of_modes,))
        assert self._file is not None
        return self._file.read("frequency")

    @property
    def phase(self):
        r"""Mode phase :math:`\varphi` [rad]."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("phase", (self.number_of_modes,))
        assert self._file is not None
        return self._file.read("phase")

    # pylint: disable=too-many-arguments
    def _stage(
        self,
        rhogrid: unyt.unyt_array,
        poloidalnumber: unyt.unyt_array,
        toroidalnumber: unyt.unyt_array,
        magneticprofile: unyt.unyt_array,
        electricprofile: unyt.unyt_array,
        amplitude: unyt.unyt_array,
        frequency: unyt.unyt_array,
        phase: unyt.unyt_array,
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.MhdStationary_init(
            ctypes.byref(self._cdata),
            poloidalnumber.size,
            rhogrid.size,
            toroidalnumber,
            poloidalnumber,
            rhogrid[[0, -1]],
            amplitude,
            frequency,
            phase,
            magneticprofile.T.ravel(),
            electricprofile.T.ravel(),
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "rhogrid",
            "poloidalnumber",
            "toroidalnumber",
            "magneticprofile",
            "electricprofile",
            "amplitude",
            "frequency",
            "phase",
        ]:
            self._file.write(field, getattr(self, field))

    def export(self) -> dict[str, unyt.unyt_array]:
        fields = [
            "rhogrid",
            "poloidalnumber",
            "toroidalnumber",
            "magneticprofile",
            "electricprofile",
            "amplitude",
            "frequency",
            "phase",
        ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            rhogrid=self.rhogrid,
            poloidalnumber=self.poloidalnumber,
            toroidalnumber=self.toroidalnumber,
            magneticprofile=self.magneticprofile,
            electricprofile=self.electricprofile,
            amplitude=self.amplitude,
            frequency=self.frequency,
            phase=self.phase,
        )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.MhdStationary_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create `MhdStationary` input."""

    # pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_mhdstationary(
        self,
        rhogrid: utils.ArrayLike,
        toroidalnumber: utils.ArrayLike,
        poloidalnumber: utils.ArrayLike,
        magneticprofile: utils.ArrayLike,
        electricprofile: utils.ArrayLike,
        amplitude: Optional[utils.ArrayLike] = None,
        frequency: Optional[utils.ArrayLike] = None,
        phase: Optional[utils.ArrayLike] = None,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
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
        preview : bool, *optional*
            If True, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : ~a5py.data.mhd.MhdStationary
            Input variant created from the given parameters.

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
        with utils.validate_variables() as v:
            rhogrid = v.validate("rhogrid", rhogrid, (-1,), "1")
            toroidalnumber = v.validate("toroidalnumber", toroidalnumber, (-1,), "1", "i4")
            poloidalnumber = v.validate("poloidalnumber", poloidalnumber, (-1,), "1", "i4")

        if toroidalnumber.size != poloidalnumber.size:
            raise ValueError(
                "There must be equal number of toroidal and poloidal modes "
                "(toroidalnumber and poloidalnumber must have the same size)."
            )

        nrho, nmode = rhogrid.size, toroidalnumber.size
        with utils.validate_variables() as v:
            phase = v.validate("phase", phase, (nmode,), "rad")
            amplitude = v.validate("amplitude", amplitude, (nmode,), "1")
            frequency = v.validate("frequency", frequency, (nmode,), "rad/s")
            magneticprofile = v.validate(
                "magneticprofile", magneticprofile, (nrho, nmode), "m"
            )
            electricprofile = v.validate(
                "electricprofile", electricprofile, (nrho, nmode), "V"
            )

        utils.validate_abscissa(rhogrid, "rhogrid")
        leaf = MhdStationary(note=note)
        leaf._stage(
            rhogrid=rhogrid,
            toroidalnumber=toroidalnumber,
            poloidalnumber=poloidalnumber,
            magneticprofile=magneticprofile,
            electricprofile=electricprofile,
            amplitude=amplitude,
            frequency=frequency,
            phase=phase,
        )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf,
            activate=activate,
            save=save,
            category="mhd",
        )
        return leaf
