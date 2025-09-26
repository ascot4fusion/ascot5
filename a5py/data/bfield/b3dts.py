"""Defines B_3DTS magnetic field input class and the corresponding factory
method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np

from a5py.data.access import InputVariant, Status, DataStruct
from ... import utils
from a5py.libascot import LIBASCOT
from ...exceptions import AscotIOException


class B_3DTS(InputVariant):
    """Tokamak with time-dependent 3D perturbation."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in B_3DTS.h."""
        _pack_ = 1
        _fields_ = [
        ]

    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="B_3DTS",
            struct=B_3DTS.Struct(),
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
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nr", "rmin", "rmax")
            return np.linspace(r0, r1, nr)
        if self._format == Format.CSTRUCT:
            return np.linspace(
                self._struct_.r_min, self._struct_.r_max, self._struct_.n_r
                ) * unyt.m
        return self._rgrid.copy()

    @property
    def phigrid(self):
        """The uniform grid in phi in which B is tabulated."""
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nphi", "phimin", "phimax")
            return np.linspace(r0, r1, nr)
        if self._format == Format.CSTRUCT:
            return (np.linspace(
                self._struct_.phi_min, self._struct_.phi_max,
                self._struct_.n_phi
                ) * unyt.rad).to("deg")
        return self._phigrid.copy()

    @property
    def zgrid(self):
        """The uniform grid in z in which B is tabulated."""
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        if self._format == Format.CSTRUCT:
            return np.linspace(
                self._struct_.z_min, self._struct_.z_max, self._struct_.n_z
                ) * unyt.m
        return self._zgrid.copy()

    @property
    def rgridpsi(self):
        """The uniform grid in R in which psi is tabulated."""
        if self._format == Format.HDF5:
            nr, r0, r1 = self._read_hdf5("nr", "rmin", "rmax")
            return np.linspace(r0, r1, nr)
        if self._format == Format.CSTRUCT:
            return np.linspace(
                self._struct_.r_min, self._struct_.r_max, self._struct_.n_r
                ) * unyt.m
        return self._rgridpsi.copy()

    @property
    def zgridpsi(self):
        """The uniform grid in z in which psi is tabulated."""
        if self._format == Format.HDF5:
            nz, z0, z1 = self._read_hdf5("nz", "zmin", "zmax")
            return np.linspace(z0, z1, nz)
        if self._format == Format.CSTRUCT:
            return np.linspace(
                self._struct_.z_min, self._struct_.z_max, self._struct_.n_z
                ) * unyt.m
        return self._zgridpsi.copy()

    @property
    def axisrz(self):
        """Magnetic axis R and z coordinates."""
        if self._format == Format.HDF5:
            return unyt.unyt_array((
                self._read_hdf5("axisr"), self._read_hdf5("axisz")
                ))
        if self._format == Format.CSTRUCT:
            return unyt.unyt_array((
                self._from_struct_("axisr", shape=(1,), units="m"),
                self._from_struct_("axisz", shape=(1,), units="m")
            ))
        return self._axisrz.copy()

    @property
    def psilimits(self):
        """Poloidal flux values on the magnetic axis and on the separatrix."""
        if self._format == Format.HDF5:
            return unyt.unyt_array((
                self._read_hdf5("psi0"), self._read_hdf5("psi1")
                ))
        if self._format == Format.CSTRUCT:
            return unyt.unyt_array((
                self._from_struct_("psi0", shape=(1,), units="Wb/m"),
                self._from_struct_("psi1", shape=(1,), units="Wb/m")
            ))
        return self._psilimits.copy()

    @property
    def psi(self):
        """Poloidal flux values on the (R,z) grid."""
        if self._format == Format.HDF5:
            return self._read_hdf5("psi")
        if self._format == Format.CSTRUCT:
            shape = (self._rgridpsi.size, self._zgridpsi.size)
            return self._from_struct_("psi", shape=shape, units="Wb/m")
        return self._psi.copy()

    @property
    def bphi(self):
        """Toroidal component of the magnetic field on the (R,phi,z) grid."""
        if self._format == Format.HDF5:
            return self._read_hdf5("bphi")
        if self._format == Format.CSTRUCT:
            shape = (self._rgrid.size, self._phigrid.size, self._zgrid.size)
            return self._from_struct_("bphi", shape=shape, units="T")
        return self._bphi.copy()

    @property
    def br(self):
        """Magnetic field R component (excl. equil. comp.) on
        the (R,phi,z) grid.
        """
        if self._format == Format.HDF5:
            return self._read_hdf5("br")
        if self._format == Format.CSTRUCT:
            shape = (self._rgrid.size, self._phigrid.size, self._zgrid.size)
            return self._from_struct_("br", shape=shape, units="T")
        return self._br.copy()

    @property
    def bz(self):
        """Magnetic field z component (excl. equil. comp.) on
        the (R,phi,z) grid.
        """
        if self._format == Format.HDF5:
            return self._read_hdf5("bz")
        if self._format == Format.CSTRUCT:
            shape = (self._rgrid.size, self._phigrid.size, self._zgrid.size)
            return self._from_struct_("bz", shape=shape, units="T")
        return self._bz.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        data["axisr"], data["axisz"] = data["axisrz"][0], data["axisrz"][1]
        data["nr"], data["nz"], data["nphi"], data["nrpsi"], data["nzpsi"] = (
            data["rgrid"].size, data["zgrid"].size, data["phigrid"].size,
            data["rgridpsi"].size, data["zgridpsi"].size
        )
        data["rmin"], data["zmin"], data["phimin"] = (
            data["rgrid"][0], data["zgrid"][0], data["phigrid"][0]
        )
        data["rminpsi"], data["zminpsi"] = (
            data["rgridpsi"][0], data["zgridpsi"][0]
        )
        data["rmax"], data["zmax"], data["phimax"] = (
            data["rgrid"][-1], data["zgrid"][-1], data["phigrid"][-1]
        )
        data["rmaxpsi"], data["zmaxpsi"] = (
            data["rgridpsi"][-1], data["zgridpsi"][-1]
        )
        data["psi0"], data["psi1"] = data["psilimits"][0], data["psilimits"][1]
        del data["rgrid"]
        del data["phigrid"]
        del data["zgrid"]
        del data["rgridpsi"]
        del data["zgridpsi"]
        del data["axisrz"]
        del data["psilimits"]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        """Return a dictionary with sufficient data to duplicate this instance.

        Returns
        -------
        data : dict[str, np.ndarray or unyt.unyt_array]
            Data that can be passed to create_b3ds to duplicate this instance.
        """
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

# pylint: disable=too-few-public-methods
class CreateB3DTSMixin():
    """Mixin class used by `Data` to create B_3DTS input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_b3dts(
            self,
            rgrid: utils.ArrayLike | None = None,
            phigrid: utils.ArrayLike | None = None,
            zgrid: utils.ArrayLike | None = None,
            tgrid: utils.ArrayLike | None = None,
            axisrz: Tuple[float, float] | None = None,
            psilimits: Tuple[float, float] | None = None,
            psi: utils.ArrayLike | None = None,
            bphi: utils.ArrayLike | None = None,
            br: utils.ArrayLike = None,
            bz: utils.ArrayLike = None,
            rgridpsi: Optional[utils.ArrayLike] | None = None,
            zgridpsi: Optional[utils.ArrayLike] | None = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> B_3DTS:
        r"""Create spline-interpolated non-axisymmetric tokamak field where the
        perturbation is time-dependent.

        This field is extension from :class:`B_3DS`, working in the same way but
        the perturbation is time-dependent.

        Note that the equilibrium remains time independent.

        Parameters
        ----------
        rgrid : array_like (nr,)
            The uniform grid in R in which B (and psi) are tabulated.
        phigrid : array_like (nphi,)
            The uniform grid in phi in which B is tabulated.
        zgrid : array_like (nz,)
            The uniform grid in z in which B (and psi) are tabulated.
        tgrid : array_like (nt,)
            The uniform grid in time in which B is tabulated.
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
        inputdata : B_3DTS
            Freshly minted input data object.
        """

