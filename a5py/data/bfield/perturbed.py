"""Defines `Bfield3D` perturbed magnetic field input class and the corresponding
factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import _variants, InputVariant, Format, TreeCreateClassMixin
from ..cstructs import interp2D_data, interp3D_data
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class Bfield3D(InputVariant):
    """Perturbed tokamak field interpolated with splines."""

    # pylint: disable=too-few-public-methods, too-many-instance-attributes
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in B_3DS.h."""
        _pack_ = 1
        _fields_ = [
            ('psi0', ctypes.c_double),
            ('psi1', ctypes.c_double),
            ('axis_r', ctypes.c_double),
            ('axis_z', ctypes.c_double),
            ('psi', interp2D_data),
            ('B_r', interp3D_data),
            ('B_phi', interp3D_data),
            ('B_z', interp3D_data),
            ]

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Bfield3D",
            struct=Bfield3D.Struct(),
            )
        self._rgrid: unyt.unyt_array
        self._phigrid: unyt.unyt_array
        self._zgrid: unyt.unyt_array
        self._axisrz: unyt.unyt_array
        self._psilimits: unyt.unyt_array
        self._psi: unyt.unyt_array
        self._bphi: unyt.unyt_array
        self._br: unyt.unyt_array
        self._bz: unyt.unyt_array
        self._rgridpsi: unyt.unyt_array
        self._zgridpsi: unyt.unyt_array

    @property
    def rgrid(self):
        """The uniform grid in R in which B is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.B_phi.x_min,
                self._struct_.B_phi.x_max,
                self._struct_.B_phi.n_x
                ) * unyt.m
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nr", "rmin", "rmax")
            return np.linspace(r0, r1, nr)
        return self._rgrid.copy()

    @property
    def phigrid(self):
        """The uniform grid in phi in which B is tabulated."""
        if self._staged:
            return (np.linspace(
                self._struct_.B_phi.y_min,
                self._struct_.B_phi.y_max,
                self._struct_.B_phi.n_y
                ) * unyt.rad).to("deg")
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nphi", "phimin", "phimax")
            return np.linspace(r0, r1, nr)
        return self._phigrid.copy()

    @property
    def zgrid(self):
        """The uniform grid in z in which B is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.B_phi.z_min,
                self._struct_.B_phi.z_max,
                self._struct_.B_phi.n_z
                ) * unyt.m
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        return self._zgrid.copy()

    @property
    def rgridpsi(self):
        """The uniform grid in R in which psi is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.psi.x_min,
                self._struct_.psi.x_max,
                self._struct_.psi.n_x
                ) * unyt.m
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nr", "rmin", "rmax")
            return np.linspace(r0, r1, nr)
        return self._rgridpsi.copy()

    @property
    def zgridpsi(self):
        """The uniform grid in z in which psi is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.psi.y_min,
                self._struct_.psi.y_max,
                self._struct_.psi.n_y
                ) * unyt.m
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        return self._zgridpsi.copy()

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
        """Toroidal component of the magnetic field on the (R,phi,z) grid."""
        if self._staged:
            return self._from_struct_("B_phi", units="T")
        if self._format == Format.HDF5:
            return self._read_hdf5("bphi")
        return self._bphi.copy()

    @property
    def br(self):
        """Magnetic field R component (excl. equil. comp.) on
        the (R,phi,z) grid.
        """
        if self._staged:
            return self._from_struct_("B_r", units="T")
        if self._format == Format.HDF5:
            return self._read_hdf5("br")
        return self._br.copy()

    @property
    def bz(self):
        """Magnetic field z component (excl. equil. comp.) on
        the (R,phi,z) grid.
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
        for grid in ["rgrid", "phigrid", "zgrid", "rgridpsi", "zgridpsi"]:
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
            "phigrid":self.phigrid,
            "zgrid":self.zgrid,
            "axisrz":self.axisrz,
            "psilimits":self.psilimits,
            "psi":self.psi,
            "bphi":self.bphi,
            "br":self.br,
            "bz":self.bz,
            "rgridpsi":self.rgridpsi,
            "zgridpsi":self.zgridpsi,
        }
        return data

    def stage(self):
        init = LIBASCOT.B_3DS_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
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
                self.rgridpsi.size,
                self.rgridpsi[0].v,
                self.rgridpsi[-1].v,
                self.zgridpsi.size,
                self.zgridpsi[0].v,
                self.zgridpsi[-1].v,
                self.rgrid.size,
                self.rgrid[0].v,
                self.rgrid[-1].v,
                self.phigrid.size,
                self.phigrid[0].to("rad").v,
                self.phigrid[-1].to("rad").v,
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
        free = LIBASCOT.B_3DS_free
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
class CreateBfield3DMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create `Bfield3D` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_bfield3d(
            self,
            rgrid: utils.ArrayLike | None = None,
            phigrid: utils.ArrayLike | None = None,
            zgrid: utils.ArrayLike | None = None,
            axisrz: Tuple[float, float] | None = None,
            psilimits: Tuple[float, float] | None = None,
            psi: utils.ArrayLike | None = None,
            bphi: utils.ArrayLike | None = None,
            br: utils.ArrayLike | None = None,
            bz: utils.ArrayLike | None = None,
            rgridpsi: Optional[utils.ArrayLike] = None,
            zgridpsi: Optional[utils.ArrayLike] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> Bfield3D:
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
            The uniform grid in R in which B (and psi) are tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in phi in which B is tabulated.
        zgrid : array_like (nz,)
            The uniform grid in z in which B (and psi) are tabulated.
        axisrz : tuple[float, float]
            Magnetic axis R and z coordinates.
        psilimits : tuple[float, float]
            Poloidal flux values on the magnetic axis and on the separatrix.
        psi : array_like (nr, nz) or (nrpsi, nzpsi)
            Poloidal flux values on the (R,z) grid.
        br : array_like (nr, nphi, nz)
            Radial component of the magnetic field on the (R,phi,z) grid.
        bphi : array_like (nr, nphi, nz)
            Toroidal component of the magnetic field on the (R,phi,z) grid.
        bz : array_like (nr, nphi, nz)
            Axial component of the magnetic field on the (R,phi,z) grid.
        rgridpsi : float, optional
            If provided, the uniform grid in R in which psi is tabulated.
        zgridpsi : float, optional
            If provided, the uniform grid in z in which psi is tabulated.
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
        inputdata : ~a5py.data.bfield.Bfield3D
            Freshly minted input data object.

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
        parameters = _variants.parse_parameters(
            rgrid, phigrid, zgrid, axisrz, psilimits, psi, bphi, br, bz,
            rgridpsi, zgridpsi,
        )
        default_rgrid, default_phigrid, default_zgrid = (
            np.linspace(1, 2, 45),
            np.linspace(0, 360, 31)[:-1],
            np.linspace(-1, 1, 90),
        )
        nr = (default_rgrid.size if parameters["rgrid"] is None
              else parameters["rgrid"].size)
        nphi = (default_phigrid.size if parameters["phigrid"] is None
              else parameters["phigrid"].size)
        nz = (default_zgrid.size if parameters["zgrid"] is None
              else parameters["zgrid"].size)
        if parameters["rgridpsi"] is None:
            parameters["rgridpsi"] = parameters["rgrid"]
        if parameters["zgridpsi"] is None:
            parameters["zgridpsi"] = parameters["zgrid"]
        nrpsi = (default_rgrid.size if parameters["rgridpsi"] is None
              else parameters["rgridpsi"].size)
        nzpsi = (default_zgrid.size if parameters["zgridpsi"] is None
              else parameters["zgridpsi"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rgrid", "phigrid", "zgrid", "axisrz", "psilimits", "psi",
                   "bphi", "br", "bz"],
            units=["m", "deg", "m", "m", "Wb/m", "Wb/m", "T", "T", "T"],
            shape=[(nr,), (nphi,), (nz,), (2,), (2,), (nrpsi, nzpsi),
                   (nr, nphi, nz), (nr, nphi, nz), (nr, nphi, nz)],
            dtype="f8",
            default=[
                default_rgrid, default_phigrid, default_zgrid, (1.5, 0.),
                (0., 1.), np.zeros((nr, nz)), np.ones((nr, nphi, nz)),
                np.ones((nr, nphi, nz)), np.ones((nr, nphi, nz)),
                ],
        )
        _variants.validate_optional_parameters(
            parameters,
            ["rgridpsi", "zgridpsi"], ("m", "m"), [(nrpsi,), (nzpsi,)], "f8",
            [parameters["rgrid"].v, parameters["zgrid"].v],
        )
        for abscissa in ["rgrid", "phigrid", "zgrid", "rgridpsi", "zgridpsi"]:
            periodic = abscissa in ["phigrid"]
            utils.check_abscissa(
                parameters[abscissa], abscissa, periodic=periodic
                )

        meta = _variants.new_metadata("Bfield3D", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
