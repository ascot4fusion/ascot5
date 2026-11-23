"""Defines the perturbed magnetic field (linear interpolation) input class and
the corresponding factory method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import (
    LIBASCOT, DataStruct, Linear2D, Linear3D, init_fun,
    )
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods, too-many-instance-attributes
class Struct(DataStruct):
    """Python wrapper for the BfieldLinear3D struct in bfield_data.h."""

    _fields_ = [
        ("axisrz", ctypes.c_double * 2),
        ("psilimits", ctypes.c_double * 2),
        ("psi", Linear2D),
        ("br", Linear3D),
        ("bz", Linear3D),
        ("bphi", Linear3D),
        ]

    _buffers: dict[str, np.ndarray]
    """Tabulated data the interpolants point to.

    The data is owned by Python and it is not copied on the C side, so these
    arrays must be kept alive as long as the struct is.
    """

    def readonly_buffer(self, name: str, units: str) -> unyt.unyt_array:
        """Return a read-only view on the tabulated data."""
        arr = self._buffers[name].view()
        arr.flags.writeable = False
        return unyt.unyt_array(arr, units)


init_fun(
    "BfieldLinear3D_init",
    ctypes.POINTER(Struct),
    *(5*[ctypes.c_size_t]),
    *(11*[ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]),
    restype=ctypes.c_int32,
    )

init_fun("BfieldLinear3D_free", ctypes.POINTER(Struct))

@Leaf.register
class BfieldLinear3D(InputVariant):
    """Perturbed tokamak field interpolated linearly."""

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
        r"""The uniform grid in :math:`\phi` in which :math:`\mathbf{B}` is
        tabulated.
        """
        if self._cdata is not None:
            return self._cdata.readonly_grid("y", "rad", "bphi", convert="deg")
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
            return self._cdata.readonly_buffer("psi", "Wb/rad")
        assert self._file is not None
        return self._file.read("psi")

    @property
    def bphi(self) -> unyt.unyt_array:
        """Tabulated values of toroidal component of the magnetic field."""
        if self._cdata is not None:
            return self._cdata.readonly_buffer("bphi", "T")
        assert self._file is not None
        return self._file.read("bphi")

    @property
    def br(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`R` component of the non-equilibrium
        magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_buffer("br", "T")
        assert self._file is not None
        return self._file.read("br")

    @property
    def bz(self) -> unyt.unyt_array:
        r"""Tabulated values of :math:`z` component of the non-equilibrium
        magnetic field.
        """
        if self._cdata is not None:
            return self._cdata.readonly_buffer("bz", "T")
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
        # The interpolants use the data in place, so it is copied here to make
        # sure it is C-contiguous float64 and not shared with the caller. The
        # struct keeps the copies alive; they are released with the struct.
        buffers = {
            name: np.array(data.v, dtype=np.float64, order="C", copy=True)
            for name, data in [("psi", psi), ("br", br), ("bz", bz),
                               ("bphi", bphi)]
            }
        for buffer in buffers.values():
            buffer.flags.writeable = False
        self._cdata = Struct()
        self._cdata._buffers = buffers
        if LIBASCOT.BfieldLinear3D_init(
            ctypes.byref(self._cdata), rgridpsi.size, zgridpsi.size, rgrid.size,
            zgrid.size, phigrid.size-1, rgridpsi[[0, -1]].v, zgridpsi[[0, -1]].v,
            rgrid[[0, -1]].v, zgrid[[0, -1]].v,
            phigrid[[0, -1]].to("rad").v,
            axisrz.v, psilimits.v, buffers["psi"], buffers["br"],
            buffers["bz"], buffers["bphi"],
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
        LIBASCOT.BfieldLinear3D_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_bfieldlinear3d(
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
            ) -> BfieldLinear3D:
        r"""Create linearly interpolated non-axisymmetric tokamak field.

        The input created by this method is equivalent to
        :class:`BfieldSpline3D` with the exception that here the magnetic field
        is interpolated linearly. Linear interpolation is faster but at the same
        time less accurate and less smooth. The poloidal flux is still
        interpolated with splines.

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
        data : :class:`.BfieldLinear3D`
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
            zgridpsi = v.validate(
                "zgridpsi", zgridpsi, (-1,), "m", default=zgrid.v
                )

        assert isinstance(rgridpsi, unyt.unyt_array)
        assert isinstance(zgridpsi, unyt.unyt_array)
        nrpsi, nzpsi = rgridpsi.size, zgridpsi.size
        nr, nphi, nz = rgrid.size, phigrid.size-1, zgrid.size
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

        leaf = BfieldLinear3D(note=note)
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
