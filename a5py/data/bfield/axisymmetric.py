"""Defines :class:`bfield2d` axisymmetric magnetic field input class and the
corresponding factory method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, interp2D_data, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods, too-many-instance-attributes
class Struct(DataStruct):
    """Python wrapper for the struct in B_2DS.h."""

    _fields_ = [
        ("psi0", ctypes.c_double),
        ("psi1", ctypes.c_double),
        ("axis_r", ctypes.c_double),
        ("axis_z", ctypes.c_double),
        ("psi", interp2D_data),
        ("B_r", interp2D_data),
        ("B_phi", interp2D_data),
        ("B_z", interp2D_data),
        ]

init_fun(
    "B_2DS_init",
    ctypes.POINTER(Struct),
    ctypes.c_int32,
    *(2*[ctypes.c_double]),
    ctypes.c_int32,
    *(6*[ctypes.c_double]),
    *(4*[ndpointer(ctypes.c_double)]),
    restype=ctypes.c_int32,
    )

init_fun("B_2DS_free", ctypes.POINTER(Struct))


@Leaf.register
class Bfield2D(InputVariant):
    """Axisymmetric tokamak field interpolated with splines."""

    _cdata: Optional[Struct]

    @property
    def rgrid(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`R` in which :math:`\mathbf{B}` and
        :math:`psi` are tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "m", "psi")
        assert self._file is not None
        return self._file.read("rgrid")

    @property
    def zgrid(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`z` in which :math:`\mathbf{B}` and
        :math:`psi` are tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "m", "psi")
        assert self._file is not None
        return self._file.read("zgrid")

    @property
    def axisrz(self) -> unyt.unyt_array:
        r"""Magnetic axis :math:`(R, z)` coordinates."""
        if self._cdata is not None:
            return unyt.unyt_array((
                self._cdata.readonly_carray("axis_r", (), "m"),
                self._cdata.readonly_carray("axis_z", (), "m"),
                ))
        assert self._file is not None
        return self._file.read("axisrz")

    @property
    def psilimits(self) -> unyt.unyt_array:
        """Poloidal flux values on the magnetic axis and on the separatrix."""
        if self._cdata is not None:
            return unyt.unyt_array((
                self._cdata.readonly_carray("psi0", (), "Wb/rad"),
                self._cdata.readonly_carray("psi1", (), "Wb/rad"),
                ))
        assert self._file is not None
        return self._file.read("psilimits")

    @property
    def psi(self) -> unyt.unyt_array:
        """Tabulated values of poloidal flux."""
        if self._cdata is not None:
            return self._cdata.readonly_interp("psi", "Wb/rad")
        assert self._file is not None
        return self._file.read("psi")

    @property
    def bphi(self) -> unyt.unyt_array:
        """Tabulated values of  toroidal component of the magnetic field
        including any perturbations.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("B_phi", "T")
        assert self._file is not None
        return self._file.read("bphi")

    @property
    def br(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`R` component of the non-equilibrium
        magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("B_r", "T")
        assert self._file is not None
        return self._file.read("br")

    @property
    def bz(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`z` component of the non-equilibrium
        magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("B_z", "T")
        assert self._file is not None
        return self._file.read("bz")

    #pylint: disable=too-many-arguments
    def _stage(
            self, rgrid: unyt.unyt_array,
            zgrid: unyt.unyt_array,
            axisrz: unyt.unyt_array,
            psilimits: unyt.unyt_array,
            psi: unyt.unyt_array,
            bphi: unyt.unyt_array,
            br: unyt.unyt_array,
            bz: unyt.unyt_array,
            ) -> None:
        self._cdata = Struct()
        if LIBASCOT.B_2DS_init(
            ctypes.byref(self._cdata), rgrid.size, rgrid[0].v, rgrid[-1].v,
            zgrid.size, zgrid[0].v, zgrid[-1].v, axisrz[0].v, axisrz[1].v,
            psilimits[0].v, psilimits[1].v, psi.v, br.v, bphi.v, bz.v,
            ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        for field in [
            "rgrid", "zgrid", "axisrz", "psilimits", "psi", "bphi", "br", "bz",
        ]:
            assert self._file is not None
            self._file.write(field, getattr(self, field))

    def export(self) -> dict[str, unyt.unyt_array | np.ndarray | int]:
        fields = [
            "rgrid", "zgrid", "axisrz", "psilimits", "psi", "bphi", "br", "bz",
        ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            rgrid=self.rgrid, zgrid=self.zgrid, axisrz=self.axisrz,
            psilimits=self.psilimits, psi=self.psi, bphi=self.bphi, br=self.br,
            bz=self.bz,
            )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.B_2DS_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_bfield2d(
            self,
            rgrid: utils.ArrayLike,
            zgrid: utils.ArrayLike,
            psi: utils.ArrayLike,
            bphi: utils.ArrayLike,
            br: Optional[utils.ArrayLike]=None,
            bz: Optional[utils.ArrayLike]=None,
            axisrz: Optional[utils.ArrayLike]=None,
            psilimits: Optional[utils.ArrayLike]= None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> Bfield2D:
        r"""Create axisymmetric tokamak field that is interpolated with splines.

        This method creates a field that is suitable for simulations where
        toroidal ripple or other 3D perturbations are not relevant [1]_. Using
        this input instead of it's 3D counterpart
        (:class:`.Bfield3D`) makes the simulations *much* faster, so using this
        input whenever applicable is strongly recommended.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in :math:`R` in which :math:`\mathbf{B}` and
            :math:`psi` are tabulated.
        zgrid : array_like (nz,)
            The uniform grid in :math:`z` in which :math:`\mathbf{B}` and
            :math:`psi` are tabulated.
        psi : array_like (nr,nz)
            Tabulated values of poloidal flux.
        bphi : array_like (nr,nz)
            Tabulated values of  toroidal component of the magnetic field.
        br : array_like (nr,nz), *optional*
            Tabulated values of :math:`R` component of the non-equilibrium
            magnetic field.

            In most cases this should not be set, see the notes below.
        bz : array_like (nr,nz), *optional*
            Tabulated values of :math:`z` component of the non-equilibrium
            magnetic field.

            In most cases this should not be set, see the notes below.
        axisrz : array_like (2,), *optional*
            Magnetic axis :math:`(R,z)` coordinates.

            The recommendation is not to set this variable in which case the
            location of the axis is interpolated internally. If the
            interpolation fails (which can happen when e.g.the data covers the
            PF coil positions), this can be used to set the correct location.
        psilimits : array_like (2,), *optional*
            Poloidal flux values on the magnetic axis and on the separatrix.

            By default the inner limit is interpolated internally using the
            ``axisrz`` position. The outer limit is assumed to be
            :math:`\psi=0`.
        note : str, *optional*
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, *optional*
            Set this input as active on creation.
        preview : bool, *optional*
            If ``True``, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        data : :class:`.Bfield2D`
            Input variant created from the given parameters.

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

        This means that by default ``br`` and ``bz`` should not contain any data
        as the equilibrium component of :math:`B_\theta` is calculated from
        the poloidal flux [2]_. This makes the field inherently divergence-free.

        .. [1] Except resonant magnetic perturbations which can be included
           via the dedicated MHD module and used together with this field.

        .. [2] In special cases where `psi` has poor quality, ``br`` and ``bz``
           can be used instead. Doing so requires that ``psi`` is scaled down to
           insignificant value so that :math:`B_\theta` can be provided
           explicitly via ``br`` and ``bz``. Remember to scale ``psilimits`` as
           well (if provided explicitly). Caution is advised as this invalidates
           the divergence-free quality of the field.
        """
        with utils.validate_variables() as v:
            rgrid = v.validate("rgrid", rgrid, (-1,), "m")
            zgrid = v.validate("zgrid", zgrid, (-1,), "m")

        nr, nz = rgrid.size, zgrid.size
        with utils.validate_variables() as v:
            psi = v.validate("psi", psi, (nr, nz), "Wb/rad")
            bphi = v.validate("bphi", bphi, (nr, nz), "T")
            br = v.validate(
                "br", br, (nr, nz), "T", default=np.full((nr, nz), 0.0)
                )
            bz = v.validate(
                "bz", bz, (nr, nz), "T", default=np.full((nr, nz), 0.0)
                )

        interpolated_psi = unyt.unyt_array(0.0, "Wb/rad")
        interpolated_axisrz = unyt.unyt_array([6.2, 0], "m")
        with utils.validate_variables() as v:
            axisrz = v.validate(
                "axisrz", axisrz, (2,), "m", default=interpolated_axisrz.v
                )
            psilimits = v.validate(
                "psilimits", psilimits, (2,), "Wb/rad",
                default=np.array([interpolated_psi.v, 0.0]),
                )

        for abscissa in ["rgrid", "zgrid"]:
            utils.check_abscissa(locals()[abscissa], abscissa)

        leaf = Bfield2D(note=note)
        leaf._stage(
            rgrid=rgrid, zgrid=zgrid, axisrz=axisrz, psilimits=psilimits,
            psi=psi, bphi=bphi, br=br, bz=bz,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="bfield",
            )
        return leaf
