"""Defines the perturbed magnetic field input class and the corresponding
factory method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import (
    LIBASCOT, DataStruct, interp2D_data, interp3D_data, init_fun,
    )
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods, too-many-instance-attributes
class Struct(DataStruct):
    """Python wrapper for the struct in B_3DS.h."""

    _fields_ = [
        ("axisrz", ctypes.c_double * 2),
        ("psilimits", ctypes.c_double * 2),
        ("psi", interp2D_data),
        ("br", interp3D_data),
        ("bz", interp3D_data),
        ("bphi", interp3D_data),
        ]


init_fun(
    "BfieldSpline3D_init",
    ctypes.POINTER(Struct),
    *(5*[ctypes.c_int32]),
    *(11*[ndpointer(ctypes.c_double)]),
    restype=ctypes.c_int32,
    )

init_fun("BfieldSpline3D_free", ctypes.POINTER(Struct))

@Leaf.register
class BfieldSpline3D(InputVariant):
    """Perturbed tokamak field interpolated with splines."""

    _cdata: Optional[Struct]

    @property
    def rgrid(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`R` in which :math:`\mathbf{B}` is
        tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "m", "bphi")
        assert self._file is not None
        return self._file.read("rgrid")

    @property
    def phigrid(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`\phi` in which :math:`psi` is tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "rad", "bphi").to("deg")
        assert self._file is not None
        return self._file.read("phigrid")

    @property
    def zgrid(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`z` in which :math:`\mathbf{B}` is
        tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("z", "m", "bphi")
        assert self._file is not None
        return self._file.read("zgrid")

    @property
    def rgridpsi(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`R` in which :math:`psi` is tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "m", "psi")
        assert self._file is not None
        return self._file.read("rgridpsi")

    @property
    def zgridpsi(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`z` in which :math:`psi` is tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "m", "psi")
        assert self._file is not None
        return self._file.read("zgridpsi")

    @property
    def axisrz(self) -> unyt.unyt_array:
        r"""Magnetic axis :math:`(R, z)` coordinates."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("axisrz", (2,), "m")
        assert self._file is not None
        return self._file.read("axisrz")

    @property
    def psilimits(self) -> unyt.unyt_array:
        """Poloidal flux values on the magnetic axis and on the separatrix."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("psilimits", (2,), "Wb/rad")
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
        """Tabulated values of toroidal component of the magnetic field."""
        if self._cdata is not None:
            return self._cdata.readonly_interp("bphi", "T")
        assert self._file is not None
        return self._file.read("bphi")

    @property
    def br(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`R` component of the non-equilibrium
        magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("br", "T")
        assert self._file is not None
        return self._file.read("br")

    @property
    def bz(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`z` component of the non-equilibrium
        magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("bz", "T")
        assert self._file is not None
        return self._file.read("bz")

    #pylint: disable=too-many-arguments
    def _stage(
            self, rgridpsi: unyt.unyt_array,
            zgridpsi: unyt.unyt_array,
            rgrid: unyt.unyt_array,
            phigrid: unyt.unyt_array,
            zgrid: unyt.unyt_array,
            axisrz: unyt.unyt_array,
            psilimits: unyt.unyt_array,
            psi: unyt.unyt_array,
            bphi: unyt.unyt_array,
            br: unyt.unyt_array,
            bz: unyt.unyt_array,
            ) -> None:
        self._cdata = Struct()
        if LIBASCOT.BfieldSpline3D_init(
            ctypes.byref(self._cdata), rgridpsi.size, zgridpsi.size, rgrid.size,
            zgrid.size, phigrid.size, rgridpsi[[0, -1]].v, zgridpsi[[0, -1]].v,
            rgrid[[0, -1]].v, zgrid[[0, -1]].v, phigrid[[0, -1]].to("rad").v,
            axisrz.v, psilimits.v, psi.v, br.v, bz.v, bphi.v,
            ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        for field in [
            "rgridpsi", "zgridpsi", "rgrid", "phigrid", "zgrid", "axisrz",
            "psilimits", "psi", "bphi", "br", "bz",
            ]:
            assert self._file is not None
            self._file.write(field, getattr(self, field))

    def export(self) -> dict[str, unyt.unyt_array | np.ndarray | int]:
        fields = [
            "rgridpsi", "zgridpsi", "rgrid", "phigrid", "zgrid", "axisrz",
            "psilimits", "psi", "bphi", "br", "bz",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            rgrid=self.rgrid, phigrid=self.phigrid, zgrid=self.zgrid,
            rgridpsi=self.rgridpsi, zgridpsi=self.zgridpsi, axisrz=self.axisrz,
            psilimits=self.psilimits, psi=self.psi, bphi=self.bphi, br=self.br,
            bz=self.bz,
            )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.BfieldSpline3D_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_bfieldspline3d(
            self,
            rgrid: utils.ArrayLike,
            phigrid: utils.ArrayLike,
            zgrid: utils.ArrayLike,
            psi: utils.ArrayLike,
            bphi: utils.ArrayLike,
            br: utils.ArrayLike,
            bz: utils.ArrayLike,
            rgridpsi: Optional[utils.ArrayLike]=None,
            zgridpsi: Optional[utils.ArrayLike]=None,
            axisrz: Optional[utils.ArrayLike]=None,
            psilimits: Optional[utils.ArrayLike]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> BfieldSpline3D:
        r"""Create spline-interpolated non-axisymmetric tokamak field.

        This method creates a tokamak magnetic field that can include arbitrary
        3D perturbations such as TF ripple or RMP coils. Using 3D magnetic field
        increases memory consumption and simulation time significantly.

        The toroidal angle is assumed to be periodic but neither the grid nor
        the data should contain the value at :math:`\phi_\mathrm{max}` to avoid
        unnecessary duplication of data. In other words, assuming

        .. math::

            A(\varphi=\varphi_\mathrm{min}) == A(\varphi=\varphi_\mathrm{max})

        then value :math:`A(\varphi=\varphi_\mathrm{max}` should not be included
        and the toroidal grid should be given by

        .. math::

            \varphi_i = \varphi_\mathrm{min}, \varphi_\mathrm{min}
            + \Delta \varphi, \dots, \varphi_\mathrm{min}
            + (n-1) \Delta \varphi,

        where
        :math:`\varphi_\mathrm{max} = \varphi_\mathrm{min} + n\Delta \varphi`
        is not included.

        It is possible to use different :math:`(R,z)` grids for :math:`\psi` and
        magnetic field components by giving a separate :math:`(R,z)` grid for
        :math:`\psi`. 3D data can be memory intensive which necessitates sparser
        grid for :math:`\mathbf{B}` components, but :math:`\psi` can still be
        evaluated on a dense grid.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in :math:`R` in which :math:`\mathbf{B}` (and
            :math:`psi`) are tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in :math:`\phi` in which :math:`\mathbf{B}` is
            tabulated.
        zgrid : array_like (nz,)
            The uniform grid in :math:`z` in which :math:`\mathbf{B}` (and
            :math:`psi`) are tabulated.
        psi : array_like (nr,nz) *or* (nrpsi,nzpsi)
            Tabulated values of poloidal flux.
        bphi : array_like (nr,nphi,nz)
            Tabulated values of  toroidal component of the magnetic field
            including any perturbations.
        br : array_like (nr,nphi,nz)
            Tabulated values of :math:`R` component of the non-equilibrium
            magnetic field.

            This should not include the axisymmetric contribution from the
            equilibrium.
        bz : array_like (nr,nphi,nz)
            Tabulated values of :math:`z` component of the non-equilibrium
            magnetic field.

            This should not include the axisymmetric contribution from the
            equilibrium.
        rgridpsi : array_like (nrpsi,), *optional*
            If provided, the uniform grid in :math:`R` in which :math:`psi` is
            tabulated.
        zgridpsi : array_like (nzpsi,), *optional*
            If provided, the uniform grid in :math:`z` in which :math:`psi` is
            tabulated.
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
        data : :class:`.BfieldSpline3D`
            Input variant created from the given parameters.

        Notes
        -----
        The input consists of poloidal flux psi (COCOS3) and the (perturbed)
        magnetic field defined in uniform :math:`(R,phi,z)` grid. During the
        simulation, the values are interpolated using bicubic splines for
        :math:`\psi`:

        .. math::

            f(x, y) = \sum_{i=0}^{3} \sum_{j=0}^{3} \
                a_{ij} (x - x_0)^i (y - y_0)^j

        and tricubic splines for the magnetic field components. The coefficients
        :math:`a_{ij}` are calculated at the start of the simulation.

        The magnetic field vector is evaluated as the sum of the explicitly
        given values of :math:`\mathbf{B}` and the contribution from the
        gradient of :math:`\psi`:

        .. math::

            \mathbf{B} = \frac{1}{R}\nabla\psi\times\hat{\mathbf{e}}_\phi \
                + B_\phi\hat{\mathbf{e}}_\phi

        This means that by default `br` and `bz` should contain only the
        perturbation component as the equilibrium component of :math:`B_\theta`
        is calculated from the poloidal flux. Since the 3D magnetic field is
        interpolated directly, the resulting field is not divergence free.
        """
        with utils.validate_variables() as v:
            rgrid = v.validate("rgrid", rgrid, (-1,), "m")
            zgrid = v.validate("zgrid", zgrid, (-1,), "m")
            phigrid = v.validate("phigrid", phigrid, (-1,), "deg")

        with utils.validate_variables() as v:
            rgridpsi = v.validate(
                "rgridpsi", rgridpsi, (-1,), "m", default=rgrid.v
                )
            zgridpsi = v.validate(
                "zgridpsi", zgridpsi, (-1,), "m", default=zgrid.v
                )

        assert isinstance(rgridpsi, unyt.unyt_array)
        assert isinstance(zgridpsi, unyt.unyt_array)
        nrpsi, nzpsi = rgridpsi.size, zgridpsi.size
        nr, nphi, nz = rgrid.size, phigrid.size, zgrid.size
        with utils.validate_variables() as v:
            psi = v.validate("psi", psi, (nrpsi, nzpsi), "Wb/rad")
            bphi = v.validate("bphi", bphi, (nr, nphi, nz), "T")
            br = v.validate("br", br, (nr, nphi, nz), "T")
            bz = v.validate("bz", bz, (nr, nphi, nz), "T")

        for abscissa in ["rgrid", "phigrid", "zgrid", "rgridpsi", "zgridpsi"]:
            periodic = abscissa in ["phigrid"]
            utils.validate_abscissa(
                locals()[abscissa], abscissa, periodic=periodic
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

        leaf = BfieldSpline3D(note=note)
        leaf._stage(
            rgrid=rgrid, phigrid=phigrid, zgrid=zgrid, rgridpsi=rgridpsi,
            zgridpsi=zgridpsi, axisrz=axisrz, psilimits=psilimits,
            psi=psi, bphi=bphi, br=br, bz=bz,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="bfield",
            )
        return leaf
