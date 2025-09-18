"""Defines `BfieldStellarator`stellarator magnetic field input class and the
corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import _variants, InputVariant, Format, TreeCreateClassMixin
from ..cstructs import linint1D_data, interp3D_data
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class BfieldStellarator(InputVariant):
    """Stellarator magnetic field interpolated with splines."""

    # pylint: disable=too-few-public-methods, too-many-instance-attributes
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in B_STS.h."""
        _pack_ = 1
        _fields_ = [
            ('psi0', ctypes.c_double),
            ('psi1', ctypes.c_double),
            ('axis_r', linint1D_data),
            ('axis_z', linint1D_data),
            ('psi', interp3D_data),
            ('B_r', interp3D_data),
            ('B_phi', interp3D_data),
            ('B_z', interp3D_data),
            ]

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="BfieldStellarator",
            struct=BfieldStellarator.Struct(),
            )
        self._rgrid: unyt.unyt_array
        self._phigrid: unyt.unyt_array
        self._zgrid: unyt.unyt_array
        self._axisgrid: unyt.unyt_array
        self._axisrz: unyt.unyt_array
        self._psilimits: unyt.unyt_array
        self._psi: unyt.unyt_array
        self._bphi: unyt.unyt_array
        self._br: unyt.unyt_array
        self._bz: unyt.unyt_array
        self._rgridpsi: unyt.unyt_array
        self._phigridpsi: unyt.unyt_array
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
    def phigridpsi(self):
        """The uniform grid in phi in which psi is tabulated."""
        if self._staged:
            return (np.linspace(
                self._struct_.psi.y_min,
                self._struct_.psi.y_max,
                self._struct_.psi.n_y
                ) * unyt.rad).to("deg")
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nphi", "phimin", "phimax")
            return np.linspace(r0, r1, nr)
        return self._phigridpsi.copy()

    @property
    def zgridpsi(self):
        """The uniform grid in z in which psi is tabulated."""
        if self._staged:
            return np.linspace(
                self._struct_.psi.z_min,
                self._struct_.psi.z_max,
                self._struct_.psi.n_z
                ) * unyt.m
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        return self._zgridpsi.copy()

    @property
    def axisgrid(self):
        """The uniform grid in phi in which axis coordinates are tabulated."""
        if self._staged:
            return (np.linspace(
                self._struct_.axis_r.x_min,
                self._struct_.axis_r.x_max,
                self._struct_.axis_r.n_x
                ) * unyt.rad).to("deg")
        if self._format == Format.HDF5:
            nx, x0, x1 = self._read_hdf5("naxis", "axismin", "axismax")
            return np.linspace(x0, x1, nx)
        return self._axisgrid.copy()

    @property
    def axisrz(self):
        """Magnetic axis R and z coordinates."""
        if self._staged:
            return unyt.unyt_array((
                self._from_struct_("axis_r", units="m"),
                self._from_struct_("axis_z", units="m")
            )).T
        if self._format == Format.HDF5:
            return unyt.unyt_array((
                self._read_hdf5("axisr"), self._read_hdf5("axisz")
                )).T
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

        for grid in ["rgrid", "phigrid", "zgrid", "axisgrid",
                     "rgridpsi", "phigridpsi", "zgridpsi"]:
            name = grid.replace("grid", "")
            data["n" + name] = data[grid].size
            data[name + "min"] = data[grid][0]
            data[name + "max"] = data[grid][-1]
            del data[grid]

        data["axisr"], data["axisz"], data["psi0"], data["psi1"] = (
            data["axisrz"][:,0], data["axisrz"][:,1],
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
            "axisgrid":self.axisgrid,
            "axisrz":self.axisrz,
            "psilimits":self.psilimits,
            "psi":self.psi,
            "bphi":self.bphi,
            "br":self.br,
            "bz":self.bz,
            "rgridpsi":self.rgridpsi,
            "phigridpsi":self.phigridpsi,
            "zgridpsi":self.zgridpsi,
        }
        return data

    def stage(self):
        init = LIBASCOT.B_STS_init
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
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
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
                self.phigridpsi.size,
                self.phigridpsi[0].to("rad").v,
                self.phigridpsi[-1].to("rad").v,
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
                self.axisgrid.size,
                self.axisgrid[0].to("rad").v,
                self.axisgrid[-1].to("rad").v,
                self.axisrz[:,0].v,
                self.axisrz[:,1].v,
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
                del self._axisrz
            self._staged = True


    def unstage(self):
        free = LIBASCOT.B_STS_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._br = self.br
                self._bz = self.bz
                self._psi = self.psi
                self._bphi = self.bphi
                self._axisrz = self.axisrz
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateBfieldStellaratorMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create BfieldStellarator input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_bfieldstellarator(
            self,
            rgrid: utils.ArrayLike | None = None,
            phigrid: utils.ArrayLike | None = None,
            zgrid: utils.ArrayLike | None = None,
            axisgrid: utils.ArrayLike | None = None,
            axisrz: utils.ArrayLike | None = None,
            psilimits: Tuple[float, float] | None = None,
            psi: utils.ArrayLike | None = None,
            bphi: utils.ArrayLike | None = None,
            br: utils.ArrayLike | None = None,
            bz: utils.ArrayLike | None = None,
            rgridpsi: Optional[utils.ArrayLike] = None,
            phigridpsi: Optional[utils.ArrayLike] = None,
            zgridpsi: Optional[utils.ArrayLike] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> BfieldStellarator:
        r"""Create spline-interpolated stellarators field.

        This method creates a field that extends
        :class:`~a5py.data.bfield.Bfield3D` to stellarators by providing support
        for 3D :math:`\psi` and 3D magnetic axis. Furthermore, the magnetic
        field is evaluated from :math:`\mathbf{B}` alone, leaving :math:`\psi`
        to act only as radial coordinate. This means that :math:`\mathbf{B}`
        should be given with great resolution to ensure low divergence, and that
        the value of :math:`\psi` does not matter that much outside the plasma
        region.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in R in which B (and psi) are tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in phi in which B (and psi) are tabulated.
        zgrid : array_like (nz,)
            The uniform grid in z in which B (and psi) are tabulated.
        axisgrid : array_like (naxis,), optional
            The uniform grid in phi in which axis coordinates are tabulated.
        axisrz : array_like (naxis,2)
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
        rgridpsi : array_like (nrpsi,), optional
            If provided, the uniform grid in R in which psi is tabulated.
        phigridpsi : array_like (nphipsi,), optional
            If provided, the uniform grid in phi in which psi is tabulated.
        zgridpsi : array_like (nzpsi,), optional
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
        inputdata : ~a5py.data.bfield.BfieldStellarator
            Freshly minted input data object.
        """
        parameters = _variants.parse_parameters(
            rgrid, phigrid, zgrid, axisgrid, axisrz, psilimits, psi,
            bphi, br, bz, rgridpsi, phigridpsi, zgridpsi,
        )
        default_rgrid, default_phigrid, default_zgrid, default_axisgrid = (
            np.linspace(1, 2, 45),
            np.linspace(0, 360, 31)[:-1],
            np.linspace(-1, 1, 90),
            np.linspace(0, 360, 11)[:-1],
        )
        nr = (default_rgrid.size if parameters["rgrid"] is None
              else parameters["rgrid"].size)
        nphi = (default_phigrid.size if parameters["phigrid"] is None
              else parameters["phigrid"].size)
        nz = (default_zgrid.size if parameters["zgrid"] is None
              else parameters["zgrid"].size)
        naxis = (default_axisgrid.size if parameters["axisgrid"] is None
              else parameters["axisgrid"].size)
        if parameters["rgridpsi"] is None:
            parameters["rgridpsi"] = parameters["rgrid"]
        if parameters["phigridpsi"] is None:
            parameters["phigridpsi"] = parameters["phigrid"]
        if parameters["zgridpsi"] is None:
            parameters["zgridpsi"] = parameters["zgrid"]
        nrpsi = (default_rgrid.size if parameters["rgridpsi"] is None
              else parameters["rgridpsi"].size)
        nphipsi = (default_phigrid.size if parameters["phigridpsi"] is None
              else parameters["phigridpsi"].size)
        nzpsi = (default_zgrid.size if parameters["zgridpsi"] is None
              else parameters["zgridpsi"].size)
        _variants.validate_required_parameters(
            parameters,
            names=["rgrid", "phigrid", "zgrid", "axisgrid", "axisrz",
                   "psilimits", "psi", "bphi", "br", "bz"],
            units=["m", "deg", "m", "deg", "m", "Wb/m", "Wb/m", "T", "T", "T"],
            shape=[(nr,), (nphi,), (nz,), (naxis,), (naxis, 2), (2,),
                   (nrpsi, nphipsi, nzpsi), (nr, nphi, nz), (nr, nphi, nz),
                   (nr, nphi, nz)],
            dtype="f8",
            default=[
                default_rgrid, default_phigrid, default_zgrid, default_axisgrid,
                np.ones((naxis, 2)), (0., 1.),
                np.zeros((nrpsi, nphipsi, nzpsi)), np.ones((nr, nphi, nz)),
                np.ones((nr, nphi, nz)), np.ones((nr, nphi, nz)),
                ],
        )
        _variants.validate_optional_parameters(
            parameters,
            ["rgridpsi", "phigridpsi", "zgridpsi"],
            ("m", "deg", "m"),
            [(nrpsi,), (nphipsi), (nzpsi,)],
            "f8",
            [parameters["rgrid"].v, parameters["phigrid"].v,
             parameters["zgrid"].v],
        )
        for abscissa in ["rgrid", "phigrid", "zgrid", "axisgrid",
                         "rgridpsi", "phigridpsi", "zgridpsi"]:
            periodic = abscissa in ["phigrid", "axisgrid", "phigridpsi"]
            utils.check_abscissa(
                parameters[abscissa], abscissa, periodic=periodic
                )

        meta = _variants.new_metadata("BfieldStellarator", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
