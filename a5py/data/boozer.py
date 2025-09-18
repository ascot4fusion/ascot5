"""Defines :class:`BoozerMap` Boozer coordinate mapping input class and the
corresponding factory method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from .access import _variants, InputVariant, Format, TreeCreateClassMixin
from .cstructs import interp2D_data
from .. import utils
from ..libascot import LIBASCOT
from ..exceptions import AscotIOException


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


class BoozerMap(InputVariant):
    """Mapping between cylindrical and Boozer coordinates."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in boozer.h."""
        #_pack_ = 1
        _fields_ = [
            ('psi_min', ctypes.c_double),
            ('psi_max', ctypes.c_double),
            ('rs', ctypes.POINTER(ctypes.c_double)),
            ('zs', ctypes.POINTER(ctypes.c_double)),
            ('nrzs', ctypes.c_int32),
            ('PADDING_0', ctypes.c_ubyte * 4),
            ('nu_psitheta', interp2D_data),
            ('theta_psithetageom', interp2D_data),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="BoozerMap",
            struct=BoozerMap.Struct(),
            )
        self._nthetag: unyt.unyt_array
        self._nthetab: unyt.unyt_array
        self._psigrid: unyt.unyt_array
        self._separatrix: unyt.unyt_array
        self._boozerpoloidal: unyt.unyt_array
        self._boozertoroidal: unyt.unyt_array

    @property
    def nthetag(self) -> int:
        r"""Number of geometric poloidal angle grid values in the data."""
        if self._staged:
            return self._struct_.theta_psithetageom.n_y
        if self._format == Format.HDF5:
            return self._read_hdf5("nthetag")
        return self._nthetag

    @property
    def nthetab(self) -> int:
        r"""Number of Boozer poloidal angle grid values in the data."""
        if self._staged:
            return self._struct_.nu_psitheta.n_y
        if self._format == Format.HDF5:
            return self._read_hdf5("nthetab")
        return self._nthetab

    @property
    def psigrid(self) -> unyt.unyt_array:
        """Radial grid in psi in which the data is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.nu_psitheta.x_min,
                self._struct_.nu_psitheta.x_max,
                self._struct_.nu_psitheta.n_x
                ) * unyt.dimensionless
        if self._format == Format.HDF5:
            nx, x0, x1 = self._read_hdf5("npsi", "psimin", "psimax")
            return np.linspace(x0, x1, nx)
        return self._psigrid.copy()

    @property
    def separatrix(self) -> unyt.unyt_array:
        """Separatrix :math:`(R,z)` coordinates."""
        if self._staged:
            return np.stack(
                (self._struct_.rs[: self._struct_.nrzs],
                 self._struct_.zs[: self._struct_.nrzs],), axis=1,
            ) * unyt.m
        if self._format == Format.HDF5:
            return np.stack(self._read_hdf5("rs", "zs"))
        return self._separatrix.copy()

    @property
    def boozertoroidal(self):
        """Boozer toroidal coordinates tabulated as a function of psi and
        the poloidal Boozer angle."""
        if self._staged:
            return self._from_struct_("nu_psitheta", units="rad")
        if self._format == Format.HDF5:
            return self._read_hdf5("boozertoroidal")
        return self._boozertoroidal.copy()

    @property
    def boozerpoloidal(self):
        """Boozer poloidal coordinates tabulated as a function of psi and
        the geometric poloidal angle."""
        if self._staged:
            data = self._from_struct_("theta_psithetageom", units="rad")
            return data[:,_NPADDING:-_NPADDING]
        if self._format == Format.HDF5:
            return self._read_hdf5("boozerpoloidal")[:,_NPADDING:-_NPADDING]
        return self._boozerpoloidal.copy()[:,_NPADDING:-_NPADDING]

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "nthetag": self.nthetag,
            "nthetab": self.nthetab,
            "psigrid": self.psigrid,
            "separatrix": self.separatrix,
            "boozerpoloidal": self.boozerpoloidal,
            "boozertoroidal": self.boozertoroidal,
        }
        return data

    def stage(self):
        init = LIBASCOT.boozer_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_int32,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ctypes.c_int32,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.exyz,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._separatrix
                del self._boozerpoloidal
                del self._boozertoroidal
            self._staged = True

    def unstage(self):
        free = LIBASCOT.boozer_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._exyz = self.exyz
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateBoozerMixin(TreeCreateClassMixin):
    """Mixin class used by :class:`Data` to create :class:`BoozerMap` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_boozermap(
            self,
            psigrid: utils.ArrayLike | None = None,
            nthetag: int | None = None,
            nthetab: int | None = None,
            boozerpoloidal: utils.ArrayLike | None = None,
            boozertoroidal: utils.ArrayLike | None = None,
            separatrix: utils.ArrayLike | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> BoozerMap:
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
        nthetag : int
            Number of geometric poloidal angle (angle between a point and the
            outer mid plane) grid values.
        nthetab : int
            Number of poloidal Boozer coordinate grid values.
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
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

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
        parameters = _variants.parse_parameters(
            psigrid, nthetag, nthetab, boozerpoloidal, boozertoroidal,
            separatrix,
        )
        default_psigrid = np.linspace(0.0, 1.0, 10)
        default_sep = np.stack(
            (np.array([0.0, 1.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0, 1.0]))
            )
        npsi = (default_psigrid.size if parameters["psigrid"] is None
             else parameters["psigrid"].size)
        nsep = (default_sep.shape[1] if parameters["separatrix"] is None
                else parameters["separatrix"].shape[0])
        nthetag = (6 if parameters["nthetag"] is None
                   else int(parameters["nthetag"]))
        nthetab = (12 if parameters["nthetab"] is None
                   else int(parameters["nthetab"]))
        _variants.validate_required_parameters(
            parameters,
            names=["psigrid", "boozerpoloidal", "boozertoroidal", "separatrix",
                   "nthetag", "nthetab",],
            units=["1", "rad", "rad", "m", "1", "1",],
            shape=[(npsi,), (npsi, nthetag), (npsi, nthetab), (nsep, 2),
                   (), ()],
            dtype=["f8", "f8", "f8", "f8", "i4", "i4"],
            default=[default_psigrid, np.zeros((npsi, nthetag)),
                     np.zeros((npsi, nthetab)), default_sep, nthetag, nthetab,],
        )

        # Extending boozerpoloidal data, see _PADDING for why we do it
        data = np.copy(parameters["boozerpoloidal"]).T
        parameters["boozerpoloidal"] = np.concatenate(
            (data, data[-1,:] + data[1:_NPADDING+1,:]) )
        parameters["boozerpoloidal"] = np.concatenate(
            (data[int(nthetag-_NPADDING-1):-1,:] - data[-1,:],
             parameters["boozerpoloidal"]) )
        parameters["boozerpoloidal"] = parameters["boozerpoloidal"].T

        meta = _variants.new_metadata("BoozerMap", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
