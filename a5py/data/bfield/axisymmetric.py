"""Defines :class:`bfield2d` axisymmetric magnetic field input class and the
corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ..cstructs import interp2D_data
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class Bfield2D(InputVariant):
    """Axisymmetric tokamak field interpolated with splines."""

    # pylint: disable=too-few-public-methods, too-many-instance-attributes
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in B_2DS.h."""
        _pack_ = 1
        _fields_ = [
            ('psi0', ctypes.c_double),
            ('psi1', ctypes.c_double),
            ('axis_r', ctypes.c_double),
            ('axis_z', ctypes.c_double),
            ('psi', interp2D_data),
            ('B_r', interp2D_data),
            ('B_phi', interp2D_data),
            ('B_z', interp2D_data),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Bfield2D",
            struct=Bfield2D.Struct(),
            )
        self._rgrid: unyt.unyt_array
        self._zgrid: unyt.unyt_array
        self._axisrz: unyt.unyt_array
        self._psilimits: unyt.unyt_array
        self._psi: unyt.unyt_array
        self._bphi: unyt.unyt_array
        self._br: unyt.unyt_array
        self._bz: unyt.unyt_array

    @property
    def rgrid(self):
        """The uniform grid in R in which B and psi are tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.psi.x_min,
                self._struct_.psi.x_max,
                self._struct_.psi.n_x
                ) * unyt.m
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nr", "rmin", "rmax")
            return np.linspace(r0, r1, nr)
        return self._rgrid.copy()

    @property
    def zgrid(self):
        """The uniform grid in z in which B and psi are tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.psi.y_min,
                self._struct_.psi.y_max,
                self._struct_.psi.n_y
                ) * unyt.m
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        return self._zgrid.copy()

    @property
    def axisrz(self):
        """Magnetic axis R and z coordinates."""
        if self._staged:
            return unyt.unyt_array((
                self._from_struct_("axis_r", shape=(1,), units="m"),
                self._from_struct_("axis_z", shape=(1,), units="m")
            )).T
        if self._format == Format.HDF5:
            return unyt.unyt_array((
                self._read_hdf5("axisr"), self._read_hdf5("axisz")
                ))
        return self._axisrz.copy()

    @property
    def psilimits(self):
        """Poloidal flux values on the magnetic axis and on the separatrix."""
        if self._staged:
            return unyt.unyt_array((
                self._from_struct_("psi0", shape=(1,), units="Wb/m"),
                self._from_struct_("psi1", shape=(1,), units="Wb/m")
            )).T
        if self._format == Format.HDF5:
            return unyt.unyt_array((
                self._read_hdf5("psi0"), self._read_hdf5("psi1")
                ))
        return self._psilimits.copy()

    @property
    def psi(self):
        """Poloidal flux values on the (R,z) grid."""
        if self._staged:
            return self._from_struct_("psi", units="Wb/m")
        if self._format == Format.HDF5:
            return self._read_hdf5("psi")
        return self._psi.copy()

    @property
    def bphi(self):
        """Toroidal component of the magnetic field on the (R,z) grid."""
        if self._staged:
            return self._from_struct_("B_phi", units="T")
        if self._format == Format.HDF5:
            return self._read_hdf5("bphi")
        return self._bphi.copy()

    @property
    def br(self):
        """Magnetic field R component (excl. equil. comp.) on the (R,z) grid.
        """
        if self._staged:
            return self._from_struct_("B_r", units="T")
        if self._format == Format.HDF5:
            return self._read_hdf5("br")
        return self._br.copy()

    @property
    def bz(self):
        """Magnetic field z component (excl. equil. comp.) on the (R,z) grid.
        """
        if self._staged:
            return self._from_struct_("B_z", units="T")
        if self._format == Format.HDF5:
            return self._read_hdf5("bz")
        return self._bz.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        for grid in ["rgrid", "zgrid"]:
            name = grid.replace("grid", "")
            data["n" + name] = data[grid].size
            data[name + "min"] = data[grid][0]
            data[name + "max"] = data[grid][-1]
            del data[grid]

        data["axisr"], data["axisz"], data["psi0"], data["psi1"] = (
            data["axisrz"][0], data["axisrz"][1],
            data["psilimits"][0], data["psilimits"][1],
        )
        del data["axisrz"]
        del data["psilimits"]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "rgrid":self.rgrid,
            "zgrid":self.zgrid,
            "axisrz":self.axisrz,
            "psilimits":self.psilimits,
            "psi":self.psi,
            "bphi":self.bphi,
            "br":self.br,
            "bz":self.bz,
        }
        return data

    def stage(self):
        init = LIBASCOT.B_2DS_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.rgrid.size,
                self.rgrid[0].v,
                self.rgrid[-1].v,
                self.zgrid.size,
                self.zgrid[0].v,
                self.zgrid[-1].v,
                self.axisrz[0].v,
                self.axisrz[1].v,
                self.psilimits[0].v,
                self.psilimits[1].v,
                self.psi.v,
                self.br.v,
                self.bphi.v,
                self.bz.v,
                ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._psi
                del self._br
                del self._bphi
                del self._bz
            self._staged = True


    def unstage(self):
        free = LIBASCOT.B_2DS_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._br = self.br
                self._bz = self.bz
                self._psi = self.psi
                self._bphi = self.bphi
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateBfield2DMixin(TreeCreateClassMixin):
    """Mixin class used by :class:`Data` to create :class:`Bfield2D` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_bfield2d(
            self,
            rgrid: utils.ArrayLike | None = None,
            zgrid: utils.ArrayLike | None = None,
            axisrz: Tuple[float, float] | None = None,
            psilimits: Tuple[float, float] | None = None,
            psi: utils.ArrayLike | None = None,
            bphi: utils.ArrayLike | None = None,
            br: Optional[utils.ArrayLike] = None,
            bz: Optional[utils.ArrayLike] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> Bfield2D:
        r"""Create axisymmetric tokamak field that is interpolated with splines.

        This method creates a field that is suitable for simulations where
        toroidal ripple or other 3D perturbations are not relevant [1]_. Using
        this input instead of it's 3D counterpart
        (:class:`~a5py.data.bfield.Bfield3D`) makes the simulations *much*
        faster, so using this input when applicable is strongly recommended.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in R in which B and psi are tabulated.
        zgrid : array_like (nz,)
            The uniform grid in z in which B and psi are tabulated.
        axisrz : tuple[float, float]
            Magnetic axis R and z coordinates.
        psilimits : tuple[float, float]
            Poloidal flux values on the magnetic axis and on the separatrix.
        psi : array_like (nr, nz)
            Poloidal flux values on the (R,z) grid.
        bphi : array_like (nr, nz)
            Toroidal component of the magnetic field on the (R,z) grid.
        br : array_like (nr, nz), optional
            Magnetic field R component (excl. equil. comp.) on the (R,z) grid.

            In most cases this should not be set as the poloidal component of
            the field is calculated from the poloidal flux.
        bz : array_like (nr, nz), optional
            Magnetic field z component (excl. equil. comp.) on the (R,z) grid.

            In most cases this should not be set as the poloidal component of
            the field is calculated from the poloidal flux.
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
        inputdata : ~a5py.data.bfield.Bfield2D
            Freshly minted input data object.

        Notes
        -----
        The input consists of poloidal flux psi (COCOS3) and the toroidal
        magnetic field defined in uniform :math:`(R,z)` grid. During the
        simulation, the values are interpolated using bicubic splines:

        .. math::

            f(x, y) = \sum_{i=0}^{3} \sum_{j=0}^{3} \
                a_{ij} (x - x_0)^i (y - y_0)^j

        where the coefficients :math:`a_{ij}` are calculated at the start of the
        simulation.

        The magnetic field vector is evaluated as the sum of the explicitly
        given values of :math:`\mathbf{B}` and the contribution from the
        gradient of :math:`\psi`:

        .. math::

            \mathbf{B} = \frac{1}{R}\nabla\psi\times\hat{\mathbf{e}}_\phi \
                + B_\phi\hat{\mathbf{e}}_\phi

        This means that by default `br` and `bz` should not contain any data as
        the equilibrium component of :math:`B_\theta` is calculated from
        the poloidal flux [2]_. This makes the field inherently divergence-free.

        .. [1] Except resonant magnetic perturbations which can be included
           via the dedicated MHD module and used together with this field.

        .. [2] In special cases where `psi` has poor quality, `br` and `bz` can
           be used instead. Doing so requires that `psi` is scaled down to
           insignificant value so that :math:`B_\theta` can be provided
           explicitly via `br` and `bz`. Remember to scale `psilimits` as well.
           Caution is advised as this invalidates the divergence-free
           quality of the field.
        """
        parameters = variants.parse_parameters(
            rgrid, zgrid, axisrz, psilimits, psi, bphi, br, bz,
        )
        default_rgrid = np.linspace(1., 2., 45)
        default_zgrid = np.linspace(-1., 1., 90)
        nr = (default_rgrid.size if parameters["rgrid"] is None
              else parameters["rgrid"].size)
        nz = (default_zgrid.size if parameters["zgrid"] is None
              else parameters["zgrid"].size)
        variants.validate_required_parameters(
            parameters,
            names=["rgrid", "zgrid", "axisrz", "psilimits", "psi", "bphi"],
            units=["m", "m", "m", "Wb/m", "Wb/m", "T"],
            shape=[(nr,), (nz,), (2,), (2,), (nr, nz), (nr, nz)],
            dtype="f8",
            default=[
                default_rgrid, default_zgrid, (1.5, 0.), (-1., 1.),
                np.zeros((nr, nz)), np.ones((nr, nz)),
                ],
        )
        variants.validate_optional_parameters(
            parameters,
            ["br", "bz"], ("T", "T"), (nr, nz), "f8",
            [np.zeros((nr, nz)), np.zeros((nr, nz))],
        )
        for abscissa in ["rgrid", "zgrid"]:
            utils.check_abscissa(parameters[abscissa], abscissa)

        meta = variants.new_metadata("Bfield2D", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
