"""Defines MhdDynamic MHD eigenmode input class and the corresponding factory
method.
"""

import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, Spline2D, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in mhdnonstat.h."""

    _fields_ = [
        ("n", ctypes.c_size_t),
        ("nmode", ctypes.POINTER(ctypes.c_int32)),
        ("mmode", ctypes.POINTER(ctypes.c_int32)),
        ("amplitude", ctypes.POINTER(ctypes.c_double)),
        ("omega", ctypes.POINTER(ctypes.c_double)),
        ("phase", ctypes.POINTER(ctypes.c_double)),
        ("alpha", ctypes.POINTER(Spline2D)),
        ("phi", ctypes.POINTER(Spline2D)),
    ]


init_fun(
    "MhdDynamic_init",
    ctypes.POINTER(Struct),
    ctypes.c_size_t,
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
    ndpointer(ctypes.c_double),
)

init_fun("MhdDynamic_free", ctypes.POINTER(Struct))


@Leaf.register
class MhdDynamic(InputVariant):
    """Time-dependent MHD eigenmode input."""

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
    def timegrid(self) -> unyt.unyt_array:
        """Time grid in in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "s", "alpha", 0)
        assert self._file is not None
        return self._file.read("timegrid")

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
                data = np.stack((data, self._cdata.readonly_interp("alpha", "m", idx=i)), axis=2)
            if nmode == 1:
                data = np.expand_dims(data, axis=2)
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
                data = np.stack((data, self._cdata.readonly_interp("phi", "V", idx=i)), axis=2)
            if nmode == 1:
                data = np.expand_dims(data, axis=2)
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
        timegrid: unyt.unyt_array,
        poloidalnumber: unyt.unyt_array,
        toroidalnumber: unyt.unyt_array,
        magneticprofile: unyt.unyt_array,
        electricprofile: unyt.unyt_array,
        amplitude: unyt.unyt_array,
        frequency: unyt.unyt_array,
        phase: unyt.unyt_array,
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.MhdDynamic_init(
            ctypes.byref(self._cdata),
            poloidalnumber.size,
            rhogrid.size,
            timegrid.size,
            toroidalnumber,
            poloidalnumber,
            rhogrid[[0, -1]],
            timegrid[[0, -1]],
            amplitude,
            frequency,
            phase,
            np.ascontiguousarray(magneticprofile.transpose((2,0,1)), dtype="f8"),
            np.ascontiguousarray(electricprofile.transpose((2,0,1)), dtype="f8"),
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "rhogrid",
            "timegrid",
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
            "timegrid",
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
            timegrid=self.timegrid,
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
        LIBASCOT.MhdDynamic_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Mixin class used by `Data` to create MhdDynamic input."""

    # pylint: disable=protected-access, too-many-arguments
    def create_mhddynamic(
        self,
        rhogrid: utils.ArrayLike,
        timegrid: utils.ArrayLike,
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
    ) -> MhdDynamic:
        r"""Create MHD eigenmode input where the modes evolve in time.

        This input is otherwise equivalent to :class:`MhdStationary` except that the
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
        with utils.validate_variables() as v:
            rhogrid = v.validate("rhogrid", rhogrid, (-1,), "1")
            timegrid = v.validate("timegrid", timegrid, (-1,), "s")
            toroidalnumber = v.validate("toroidalnumber", toroidalnumber, (-1,), "1", "i4")
            poloidalnumber = v.validate("poloidalnumber", poloidalnumber, (-1,), "1", "i4")

        if toroidalnumber.size != poloidalnumber.size:
            raise ValueError(
                "There must be equal number of toroidal and poloidal modes "
                "(toroidalnumber and poloidalnumber must have the same size)."
            )

        nrho, ntime, nmode = rhogrid.size, timegrid.size, toroidalnumber.size
        with utils.validate_variables() as v:
            phase = v.validate("phase", phase, (nmode,), "rad")
            amplitude = v.validate("amplitude", amplitude, (nmode,), "1")
            frequency = v.validate(
                "frequency", frequency, (nmode,), "rad/s"
            )
            magneticprofile = v.validate(
                "magneticprofile", magneticprofile, (nrho, ntime, nmode), "m"
            )
            electricprofile = v.validate(
                "electricprofile", electricprofile, (nrho, ntime, nmode), "V"
            )

        utils.validate_abscissa(rhogrid, "rhogrid")
        leaf = MhdDynamic(note=note)
        leaf._stage(
            rhogrid=rhogrid,
            timegrid=timegrid,
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
