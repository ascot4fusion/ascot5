"""Defines stellarator magnetic field input class and the corresponding factory
method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import (
    LIBASCOT, DataStruct, linint1D_data, interp3D_data, init_fun,
    )
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods, too-many-instance-attributes
class Struct(DataStruct):
    """Python wrapper for the struct in B_STS.h."""

    _fields_ = [
        ("psilimits", ctypes.c_double * 2),
        ("axisr", linint1D_data),
        ("axisz", linint1D_data),
        ("psi", interp3D_data),
        ("br", interp3D_data),
        ("bz", interp3D_data),
        ("bphi", interp3D_data),
        ]


init_fun(
    "BfieldStellarator_init",
    ctypes.POINTER(Struct),
    *(7*[ctypes.c_int32]),
    *(14*[ndpointer(ctypes.c_double)]),
    restype=ctypes.c_int32,
    )

init_fun("BfieldStellarator_free", ctypes.POINTER(Struct))

@Leaf.register
class BfieldStellarator(InputVariant):
    """Stellarator magnetic field interpolated with splines."""

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
    def phigridpsi(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`\phi` in which :math:`psi` is tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "rad", "psi").to("deg")
        assert self._file is not None
        return self._file.read("phigridpsi")

    @property
    def zgridpsi(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`z` in which :math:`psi` is tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("z", "m", "psi")
        assert self._file is not None
        return self._file.read("zgridpsi")

    @property
    def axisgrid(self) -> unyt.unyt_array:
        r"""The uniform grid in :math:`\phi` in which axis :math:`(R,z)`
        coordinates are tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("x", "rad", "axisr").to("deg")
        assert self._file is not None
        return self._file.read("axisgrid")

    @property
    def axisrz(self) -> unyt.unyt_array:
        r"""Tabulated magnetic axis :math:`(R,z)` coordinates."""
        if self._cdata is not None:
            return unyt.unyt_array((
                self._cdata.readonly_interp("axisr", "m"),
                self._cdata.readonly_interp("axisz", "m"),
                )).T
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
        r"""Tabulated values of :math:`R` component of the magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("br", "T")
        assert self._file is not None
        return self._file.read("br")

    @property
    def bz(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`z` component of the magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_interp("bz", "T")
        assert self._file is not None
        return self._file.read("bz")

     #pylint: disable=too-many-arguments
    def _stage(
            self, rgridpsi: unyt.unyt_array,
            phigridpsi: unyt.unyt_array,
            zgridpsi: unyt.unyt_array,
            rgrid: unyt.unyt_array,
            phigrid: unyt.unyt_array,
            zgrid: unyt.unyt_array,
            axisgrid: unyt.unyt_array,
            axisrz: unyt.unyt_array,
            psilimits: unyt.unyt_array,
            psi: unyt.unyt_array,
            bphi: unyt.unyt_array,
            br: unyt.unyt_array,
            bz: unyt.unyt_array,
            ) -> None:
        self._cdata = Struct()
        if LIBASCOT.BfieldStellarator_init(
            ctypes.byref(self._cdata), rgridpsi.size, zgridpsi.size,
            phigridpsi.size, rgrid.size, zgrid.size, phigrid.size, axisgrid.size,
            rgridpsi[[0, -1]].v, zgridpsi[[0, -1]].v,
            phigridpsi[[0, -1]].to("rad").v, rgrid[[0, -1]].v, zgrid[[0, -1]].v,
            phigrid[[0, -1]].to("rad").v, axisgrid[[0, -1]].to("rad").v,
            axisrz[:,0].v, axisrz[:,1].v, psilimits.v, psi.v, br.v, bz.v,
            bphi.v,
            ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        for field in [
            "rgridpsi", "phigridpsi", "zgridpsi", "rgrid", "phigrid", "zgrid",
            "axisgrid", "axisrz",  "psilimits", "psi", "bphi", "br", "bz",
            ]:
            assert self._file is not None
            self._file.write(field, getattr(self, field))

    def export(self) -> dict[str, unyt.unyt_array | np.ndarray | int]:
        fields = [
            "rgridpsi", "phigridpsi", "zgridpsi", "rgrid", "phigrid", "zgrid",
            "axisgrid", "axisrz",  "psilimits", "psi", "bphi", "br", "bz",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            rgrid=self.rgrid, phigrid=self.phigrid, zgrid=self.zgrid,
            rgridpsi=self.rgridpsi, phigridpsi=self.phigridpsi,
            zgridpsi=self.zgridpsi, axisgrid=self.axisgrid, axisrz=self.axisrz,
            psilimits=self.psilimits, psi=self.psi, bphi=self.bphi, br=self.br,
            bz=self.bz,
            )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.BfieldStellarator_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_bfieldstellarator(
            self,
            rgrid: utils.ArrayLike,
            phigrid: utils.ArrayLike ,
            zgrid: utils.ArrayLike,
            psi: utils.ArrayLike,
            bphi: utils.ArrayLike,
            br: utils.ArrayLike,
            bz: utils.ArrayLike,
            rgridpsi: Optional[utils.ArrayLike]=None,
            phigridpsi: Optional[utils.ArrayLike]=None,
            zgridpsi: Optional[utils.ArrayLike]=None,
            axisgrid: Optional[utils.ArrayLike]=None,
            axisrz: Optional[utils.ArrayLike]=None,
            psilimits: Optional[utils.ArrayLike]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> BfieldStellarator:
        r"""Create spline-interpolated stellarators field.

        This method creates a field that extends
        :class:`~a5py.data.bfield.BfieldSpline3D` to stellarators by providing support
        for 3D :math:`\psi` and 3D magnetic axis. Furthermore, the magnetic
        field is evaluated from :math:`\mathbf{B}` alone, leaving :math:`\psi`
        to act only as radial coordinate. This means that :math:`\mathbf{B}`
        should be given with great resolution to ensure low divergence, and that
        the value of :math:`\psi` does not matter that much outside the plasma
        region.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in :math:`R` in which :math:`\mathbf{B}` (and
            :math:`psi`) are tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in :math:`\phi` in which :math:`\mathbf{B}` (and
            :math:`psi`) are tabulated
        zgrid : array_like (nz,)
            The uniform grid in :math:`z` in which :math:`\mathbf{B}` (and
            :math:`psi`) are tabulated.
        psi : array_like (nr,nphi,nz) *or* (nrpsi,nphipsi,nzpsi)
            Tabulated values of poloidal flux.
        br : array_like (nr,nphi,nz)
            Tabulated values of :math:`R` component of the magnetic field.
        bphi : array_like (nr,nphi,nz)
            Tabulated values of  toroidal component of the magnetic field.
        bz : array_like (nr,nphi,nz)
            Tabulated values of :math:`z` component of the magnetic field.
        rgridpsi : array_like (nrpsi,), *optional*
            If provided, the uniform grid in :math:`R` in which :math:`psi` is
            tabulated.
        phigridpsi : array_like (nphipsi,), *optional*
            If provided, the uniform grid in :math:`\phi` in which :math:`psi`
            is tabulated.
        zgridpsi : array_like (nzpsi,), *optional*
            If provided, the uniform grid in :math:`z` in which :math:`psi` is
            tabulated.
        axisgrid : array_like (naxis,), *optional*
            The uniform grid in :math:`\phi` in which axis coordinates are
            tabulated.

            If not given, ``phigrid`` is used.
        axisrz : array_like (naxis,2), *optional*
            Magnetic axis :math:`(R,z)` coordinates as a function of
            :math:`\phi`.

            The recommendation is not to set this variable in which case the
            location of the axis is interpolated internally. If the
            interpolation fails this can be used to set the correct location.
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
        data : :class:`.BfieldStellarator`
            Input variant created from the given parameters.
        """
        with utils.validate_variables() as v:
            rgrid = v.validate("rgrid", rgrid, (-1,), "m")
            zgrid = v.validate("zgrid", zgrid, (-1,), "m")
            phigrid = v.validate("phigrid", phigrid, (-1,), "deg")

        with utils.validate_variables() as v:
            rgridpsi = v.validate(
                "rgridpsi", rgridpsi, (-1,), "m", default=rgrid.v
                )
            phigridpsi = v.validate(
                "phigridpsi", phigridpsi, (-1,), "deg", default=phigrid.v
                )
            zgridpsi = v.validate(
                "zgridpsi", zgridpsi, (-1,), "m", default=zgrid.v
                )

        assert isinstance(rgridpsi, unyt.unyt_array)
        assert isinstance(zgridpsi, unyt.unyt_array)
        assert isinstance(phigridpsi, unyt.unyt_array)
        nr, nphi, nz = rgrid.size, phigrid.size, zgrid.size
        nrpsi, nphipsi, nzpsi = rgridpsi.size, phigridpsi.size, zgridpsi.size
        with utils.validate_variables() as v:
            psi = v.validate("psi", psi, (nrpsi, nphipsi, nzpsi), "Wb/rad")
            bphi = v.validate("bphi", bphi, (nr, nphi, nz), "T")
            br = v.validate("br", br, (nr, nphi, nz), "T")
            bz = v.validate("bz", bz, (nr, nphi, nz), "T")

        with utils.validate_variables() as v:
            axisgrid = v.validate(
                "axisgrid", axisgrid, (-1,), "deg", default=phigrid.v
                )

        assert isinstance(axisgrid, unyt.unyt_array)
        interpolated_psi = unyt.unyt_array(0.0, "Wb/rad")
        interpolated_axisrz = unyt.unyt_array(np.ones(axisgrid.shape), "m")
        with utils.validate_variables() as v:
            axisrz = v.validate(
                "axisrz", axisrz, (axisgrid.size, 2), "m",
                default=interpolated_axisrz.v
                )
            psilimits = v.validate(
                "psilimits", psilimits, (2,), "Wb/rad",
                default=np.array([interpolated_psi.v, 0.0])
                )

        for abscissa in [
            "rgrid", "phigrid", "zgrid", "rgridpsi", "phigridpsi", "zgridpsi"
            ]:
            periodic = abscissa in ["phigrid", "phigridpsi"]
            utils.validate_abscissa(
                locals()[abscissa], abscissa, periodic=periodic
                )

        leaf = BfieldStellarator(note=note)
        leaf._stage(
            rgrid=rgrid, phigrid=phigrid, zgrid=zgrid, rgridpsi=rgridpsi,
            phigridpsi=phigridpsi, zgridpsi=zgridpsi, axisgrid=axisgrid,
            axisrz=axisrz, psilimits=psilimits, psi=psi, bphi=bphi, br=br,
            bz=bz,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="bfield",
            )
        return leaf
