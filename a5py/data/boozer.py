"""Defines :class:`Boozer` Boozer coordinate mapping input class and the
corresponding factory method.
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

_NPADDING = 4
"""How many indices are used to "pad" the Boozer poloidal angle data to extend
it beoynd [0,2pi].

The Boozer poloidal angle is interpolated with splines. Since it is a periodic
quantity, it would make sense to use the periodic boundary condition (in the
axis corresponding to the geometrical poloidal angle) BUT that boundary quantity
is for continuous quantities (whereas the Boozer poloidal angle is cyclic).

Therefore we must use the natural boundary condition. Since this would make the
splines inaccurate near the boundary, we instead add some extra values to both
ends (where we never actually interpolate as they are outside the range [0,2pi])
so that the interpolation happens far enough that the boundary condition doesn't
have an effect. This value is the number of points that we add on both ends.

Seeing how long this explanation is, there should be a better way to do this.
"""


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in boozer.h."""

    _fields_ = [
        ("nrz", ctypes.c_size_t),
        ("rlim", ctypes.POINTER(ctypes.c_double)),
        ("zlim", ctypes.POINTER(ctypes.c_double)),
        ("theta", Spline2D),
        ("nu", Spline2D),
    ]


init_fun(
    "Boozer_init",
    ctypes.POINTER(Struct),
    *(4 * [ctypes.c_size_t]),
    ctypes.c_int32,
    *(5 * [ndpointer(ctypes.c_double)]),
)

init_fun("Boozer_free", ctypes.POINTER(Struct))


@Leaf.register
class Boozer(InputVariant):
    """Mapping between cylindrical and Boozer coordinates."""

    @property
    def psigrid(self) -> unyt.unyt_array:
        """Radial grid in psi in which the data is tabulated."""
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "1", "nu")
        assert self._file is not None
        return self._file.read("psigrid")

    @property
    def separatrix(self) -> unyt.unyt_array:
        """Separatrix :math:`(R,z)` coordinates."""
        if self._cdata is not None:
            r = self._cdata.readonly_carray("rlim", (self._cdata.nrz,), "m")
            z = self._cdata.readonly_carray("zlim", (self._cdata.nrz,), "m")
            return np.stack((r, z), axis=1).T
        assert self._file is not None
        return self._file.read("separatrix")

    @property
    def boozertoroidal(self):
        """Boozer toroidal coordinates tabulated as a function of psi and
        the poloidal Boozer angle."""
        if self._cdata is not None:
            return self._cdata.readonly_interp("nu", "rad")
        assert self._file is not None
        return self._file.read("boozertoroidal")

    @property
    def boozerpoloidal(self):
        """Boozer poloidal coordinates tabulated as a function of psi and
        the geometric poloidal angle."""
        if self._cdata is not None:
            data = self._cdata.readonly_interp("theta", "rad")
            return data[:, _NPADDING:-_NPADDING]
        assert self._file is not None
        return self._file.read("boozerpoloidal")[:, _NPADDING:-_NPADDING]

    def _stage(
        self,
        psigrid: unyt.unyt_array,
        boozerpoloidal: unyt.unyt_array,
        boozertoroidal: unyt.unyt_array,
        separatrix: unyt.unyt_array,
    ) -> None:
        self._cdata = Struct()
        if LIBASCOT.Boozer_init(
            ctypes.byref(self._cdata),
            psigrid.size,
            boozertoroidal.shape[1],
            boozerpoloidal.shape[1],
            separatrix.shape[1],
            _NPADDING,
            psigrid[[0, -1]],
            boozertoroidal,
            boozerpoloidal,
            separatrix[:, 0],
            separatrix[:, 1],
        ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("psigrid", self.psigrid)
        self._file.write("separatrix", self.separatrix)
        self._file.write("boozerpoloidal", self.boozerpoloidal)
        self._file.write("boozertoroidal", self.boozertoroidal)

    def export(self) -> dict[str, unyt.unyt_array]:
        data = {
            "psigrid": self.psigrid,
            "separatrix": self.separatrix,
            "boozerpoloidal": self.boozerpoloidal,
            "boozertoroidal": self.boozertoroidal,
        }
        return data

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.Boozer_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateBoozerMixin(TreeMixin):
    """Mixin class used by :class:`Data` to create :class:`Boozer` input."""

    # pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_boozer(
        self,
        psigrid: utils.ArrayLike,
        boozerpoloidal: utils.ArrayLike,
        boozertoroidal: utils.ArrayLike,
        separatrix: utils.ArrayLike,
        note: Optional[str] = None,
        activate: bool = False,
        preview: bool = False,
        save: Optional[bool] = None,
    ) -> Boozer:
        r"""Create an input that implements a mapping between the cylindrical
        and the Boozer coordinates.

        The mapping between the cylindrical and Boozer coordinates is required
        in simulations where the input with MHD eigenfunctions is used.

        Note that due to limitations of numerical method, the coordinate mapping
        can be ill-defined close to the magnetic axis and near the separatrix.
        Because of this, it makes sense to limit ``psigrid`` to only where the
        mapping is valid. Markers are **not** aborted if they exit this region
        as the code just assumes that there's no MHD present (but they are
        aborted if the MHD data does not cover the whole ``psigrid``).

        Parameters
        ----------
        psigrid : float
            The uniform grid in psi in which the Boozer poloidal and toroidal
            coordinates are tabulated.
        boozerpoloidal : array_like (npsi, nthetag)
            Boozer poloidal coordinates tabulated as a function of psi and
            the geometric poloidal angle.
        boozertoroidal : array_like (npsi, nthetab)
            Boozer toroidal coordinates tabulated as a function of psi and
            the poloidal Boozer angle.

            To be specific, this is the difference between the Boozer toroidal
            angle and the geometric toroidal angle, which depends only on
            psi and the Boozer poloidal angle.
        separatrix : array_like (n,2)
            Separatrix :math:`(R,z)` coordinates where the first and last points
            coincide.

            Boozer coordinates are not defined outside the separatrix, so with
            this we can separate the actual plasma region from e.g. the private
            plasma region.
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
        inputdata : ~a5py.data.wall.BoozerMap
            Freshly minted input data object.

        Notes
        -----
        During the simulation, the marker cylindrical coordinates are mapped to
        straight-field line coordinates if the MHD perturbations are enabled.
        This mapping is implemented only for stationary tokamak fields and we
        further assume that the field is axisymmetric (but the code allows
        mapping to be used in non-axisymmetric fields as well).

        Our choice of the coordinate system are the Boozer coordinates
        :math:`(\psi,\theta,\zeta)`[1]_, where :math:`\psi` is the normalized
        poloidal flux, :math:`\theta` is the Boozer poloidal angle (which points
        in same direction as the geometrical poloidal angle
        :math:`\theta_\mathrm{geo}` i.e. counter-clockwise when looking at the
        same direction as positive :math:`\hat{\phi}`), and
        :math:`\zeta = \phi - \nu`, where :math:`\nu=\nu(\psi,\theta)`, is the
        Boozer toroidal angle (with the same positive direction as the
        cylindrical toroidal angle). Both Boozer angular coordinates have the
        periodicity of :math:`2\pi`.

        To faciliate the mapping in run-time, we precalculate
        :math:`\theta(\psi,\theta_\mathrm{geo})` and :math:`\nu(\psi,\theta)`
        in an uniform grid and use the tabulated values together with the
        cubic-spline interpolation to perform the mapping. This is the purpose
        of this input.

        [1] For brevity, here we use :math:`\theta` for the Boozer poloidal
            angle and :math:`\phi` for the Boozer toroidal angle. Normally these
            symbols refer to the geometrical poloidal and toroidal angle,
            respectively.
        """

        with utils.validate_variables() as v:
            psigrid = v.validate("psigrid", psigrid, (-1,), "1")
            separatrix = v.validate("separatrix", separatrix, (2, 4), "m")

        npsi = psigrid.size
        with utils.validate_variables() as v:
            boozerpoloidal = v.validate(
                "boozerpoloidal", boozerpoloidal, (npsi, -1), "rad"
            )
            boozertoroidal = v.validate(
                "boozertoroidal", boozertoroidal, (npsi, -1), "rad"
            )

        utils.validate_abscissa(psigrid, "psigrid")

        # Extending boozerpoloidal data, see _PADDING for why we do it
        nthetag = boozerpoloidal.shape[1]
        data = np.copy(boozerpoloidal).T
        boozerpoloidal = np.concatenate(
            (data, data[-1, :] + data[1 : _NPADDING + 1, :])
        )
        boozerpoloidal = np.concatenate(
            (data[int(nthetag - _NPADDING - 1) : -1, :] - data[-1, :], boozerpoloidal)
        )
        boozerpoloidal = boozerpoloidal.T

        leaf = Boozer(note=note)
        leaf._stage(
            psigrid,
            boozerpoloidal,
            boozertoroidal,
            separatrix,
        )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf,
            activate=activate,
            save=save,
            category="boozer",
        )
        return leaf
