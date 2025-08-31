"""Defines `BfieldAnalytical` analytical magnetic field input class and
the corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import mpi4py as MPI
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import variants, InputVariant, Format, TreeCreateClassMixin
from ... import utils
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException
from ...physlib import aeq


class BfieldAnalytical(InputVariant):
    """Analytical tokamak field which can be either axisymmetric or 3D."""

    # pylint: disable=too-few-public-methods, too-many-instance-attributes
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in B_GS.h."""
        _pack_ = 1
        _fields_ = [
            ('R0', ctypes.c_double),
            ('z0', ctypes.c_double),
            ('raxis', ctypes.c_double),
            ('zaxis', ctypes.c_double),
            ('B_phi0', ctypes.c_double),
            ('psi0', ctypes.c_double),
            ('psi1', ctypes.c_double),
            ('psi_mult', ctypes.c_double),
            ('psi_coeff', ctypes.c_double * 13),
            ('Nripple', ctypes.c_int32),
            ('PADDING_0', ctypes.c_ubyte * 4),
            ('a0', ctypes.c_double),
            ('alpha0', ctypes.c_double),
            ('delta0', ctypes.c_double),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="BfieldAnalytical",
            struct=BfieldAnalytical.Struct(),
            )
        self._rmajor: unyt.unyt_array
        self._axisr: unyt.unyt_array
        self._axisz: unyt.unyt_array
        self._axisb: unyt.unyt_array
        self._psiscaling: unyt.unyt_array
        self._coefficients: np.ndarray
        self._psilimits: unyt.unyt_array
        self._nripple: int
        self._rminor: unyt.unyt_array
        self._ripplepenetration: unyt.unyt_array
        self._ripplescaling: unyt.unyt_array

    @property
    def rmajor(self) -> unyt.unyt_array:
        """Major radius of the tokamak."""
        if self._staged:
            return self._from_struct_("R0", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("r0")
        return self._rmajor.copy()

    @property
    def axisr(self) -> unyt.unyt_array:
        """Magnetic axis R coordinate."""
        if self._staged:
            return self._from_struct_("raxis", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("axisr")
        return self._axisr.copy()

    @property
    def axisz(self) -> unyt.unyt_array:
        """Magnetic axis z coordinate."""
        if self._staged:
            return self._from_struct_("zaxis", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("z0")
        return self._axisz.copy()

    @property
    def axisb(self) -> unyt.unyt_array:
        """Toroidal field strength at the major radius location."""
        if self._staged:
            return self._from_struct_("B_phi0", shape=(1,), units="T")
        if self._format == Format.HDF5:
            return self._read_hdf5("bphi0")
        return self._axisb.copy()

    @property
    def psiscaling(self) -> unyt.unyt_array:
        """Scaling factor for the poloidal flux."""
        if self._staged:
            return self._from_struct_("psi_mult", shape=(1,), units="Wb/m")
        if self._format == Format.HDF5:
            return self._read_hdf5("psimult")
        return self._psiscaling.copy()

    @property
    def coefficients(self) -> unyt.unyt_array:
        """Coefficients defining psi: [c0, c1, ..., c11, A]."""
        if self._staged:
            return self._from_struct_("psi_coeff", shape=(13,), units="1")
        if self._format == Format.HDF5:
            return self._read_hdf5("coefficients")
        return self._coefficients.copy()

    @property
    def psilimits(self) -> unyt.unyt_array:
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
    def nripple(self) -> unyt.unyt_array:
        """Number of TF coils."""
        if self._staged:
            return self._from_struct_("Nripple", shape=(1,), units="1")
        if self._format == Format.HDF5:
            return self._read_hdf5("nripple")
        return self._nripple.copy()

    @property
    def rminor(self) -> unyt.unyt_array:
        """Minor radius."""
        if self._staged:
            return self._from_struct_("a0", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("a0")
        return self._rminor.copy()

    @property
    def ripplepenetration(self) -> unyt.unyt_array:
        """Ripple penetration."""
        if self._staged:
            return self._from_struct_("alpha0", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("alpha0")
        return self._ripplepenetration.copy()

    @property
    def ripplescaling(self) -> unyt.unyt_array:
        """Ripple scaling parameter."""
        if self._staged:
            return self._from_struct_("delta0", shape=(1,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("delta0")
        return self._ripplescaling.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        for name, hdf5name in (
            ("rmajor", "r0"),
            ("axisz", "z0"),
            ("axisb", "bphi0"),
            ("psiscaling", "psimult"),
            ("rminor", "a0"),
            ("ripplepenetration", "alpha0"),
            ("ripplescaling", "delta0"),
        ):
            data[hdf5name] = data[name]
            del data[name]
        data["psi0"] = data["psilimits"][0]
        data["psi1"] = data["psilimits"][1]
        del data["psilimits"]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "rmajor":self.rmajor,
            "axisz":self.axisz,
            "axisr":self.axisr,
            "axisb":self.axisb,
            "psiscaling":self.psiscaling,
            "coefficients":self.coefficients,
            "psilimits":self.psilimits,
            "nripple":self.nripple,
            "rminor":self.rminor,
            "ripplepenetration":self.ripplepenetration,
            "ripplescaling":self.ripplescaling,
        }
        return data

    def stage(self):
        init = LIBASCOT.B_GS_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ndpointer(ctypes.c_double),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.rmajor[0].v,
                self.axisz[0].v,
                self.axisr[0].v,
                self.axisz[0].v,
                self.axisb[0].v,
                self.psilimits[0].v,
                self.psilimits[1].v,
                self.psiscaling[0].v,
                self.coefficients.v,
                self.nripple[0].v,
                self.rminor[0].v,
                self.ripplepenetration[0].v,
                self.ripplescaling[0].v,
                ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._coefficients
            self._staged = True

    def unstage(self):
        free = LIBASCOT.B_GS_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._coefficients = self.coefficients

            free(ctypes.byref(self._struct_))
            self._staged = False

    def bcast(self, isroot=0):
        """Broadcast data."""
        s = self._struct_
        buf = (ctypes.c_char * ctypes.sizeof(BfieldAnalytical.Struct)).from_buffer(s)
        MPI.COMM_WORLD.Bcast([buf, MPI.BYTE], root=0)


# pylint: disable=too-few-public-methods
class CreateBfieldAnalyticalMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create `BfieldAnalytical` input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_bfieldanalytical(
            self,
            rmajor: float | None = None,
            axisz: float | None = None,
            axisb: float | None = None,
            psiscaling: float | None = None,
            coefficients : utils.ArrayLike | None = None,
            axisr: Optional[float] = None,
            psilimits: Optional[Tuple[float,float]] = None,
            nripple: Optional[int] = None,
            rminor: Optional[float] = None,
            ripplepenetration: Optional[float] = None,
            ripplescaling: Optional[float] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> BfieldAnalytical:
        r"""Create an analytic tokamak field.

        Parameters
        ----------
        rmajor : float
            Major radius of the tokamak.

            Note that this is not the same as the magnetic axis R coordinate.
        axisz : float
            Magnetic axis z coordinate.
        axisb : float
            Toroidal field strength at the major radius location.
        psiscaling : float
            Scaling factor for the poloidal flux.

            This scaling directly determines the magnitude of the plasma
            current and the poloidal field.
        coefficients : array_like (13,1)
            Coefficients defining psi: [c0, c1, ..., c11, A].

            These can be evaluated from the equilibrium properties such as
            triangularity and elongation.
        raxis : float, optional
            R coordinate of the magnetic axis.

            If not provided, this value will be interpolated internally.
        psilimits : tuple[float, float], optional
            Poloidal flux values on the magnetic axis and on the separatrix.
        nripple : float, optional
            Number of TF coils.

            If this value is set and it is larger than zero, the resulting field
            will contain the toroidal field ripple.
        rminor : float, optional
            Minor radius.

            Used only for computing ripple strength.
        ripplepenetration : float, optional
            Ripple penetration.

            Increasing this parameter decreases ripple strength as
            :math:`\propto\frac{r}{a}^{\alpha}`, where :math:`a` is the minor
            radius and :math:`\alpha` is this parameter.
        ripplescaling : float, optional
            Ripple scaling parameter.

            Corresponds to the strength of the ripple at the separatrix.

        Returns
        -------
        inputdata : ~a5py.data.bfield.BfieldAnalytical
            Freshly minted input data object.

        Notes
        -----

        This field can either be axisymmetric or it can include a ripple-like
        perturbation. The field is very fast to interpolate and it is not memory
        intensive. However, it is still intended mostly for testing purposes as
        it is difficult to recreate experimental equilibria with it.

        The axisymmetric field is divergence-free whereas the ripple field is
        not because it incorporates only the toroidal component of the
        perturbation.

        This field implements an analytical solution to the Grad-Shafranov
        equation [1]_. In this model, the poloidal flux psi is calculated as

        .. math::

            \psi(R,Z) =
            &\psi_c [(1-A) (r^4/8)                                         \\
            &+ A (r^2\log(r)/2)                                            \\
            &+ c_0                                                         \\
            &+ c_1 (r^2)                                                   \\
            &+ c_2 (r^2\log(r) - z^2)                                      \\
            &+ c_3 (r^4 - 4r^2z^2)                                         \\
            &+ c_4 (3r^4\log(r) - 9r^2z^2 - 12r^2\log(r)z^2 + 2z^4)        \\
            &+ c_5 (r^6 - 12r^4z^2 + 8r^2z^4)                              \\
            &+ c_6 (8z^6 - 140r^2z^4 - 120r^2\log(r)z^4 + 180r^4\log(r)z^2
                                    + 75r^4 z^2 - 15r^6\log(r))            \\
            &+ c_7  (z)                                                    \\
            &+ c_8  (z r^2)                                                \\
            &+ c_9  (z^3 - 3zr^2\log(r))                                   \\
            &+ c_{10} (3zr^4 - 4z^3r^2)                                    \\
            &+ c_{12} (8z^5 - 45zr^4 - 80z^3r^2\log(r) + 60zr^4\log(r)) ]

        where :math:`c_i` and :math:`A` are pre-calculated coefficients which
        can be chosen so that realistic equilibria resembling different machines
        are produced. The equilibrium can be non-symmetric with respect to
        magnetic plane, and can have zero, one, or two X-points. :math:`\psi_c`
        is a scaling constant, and :math:`r = R/R_0` :math:`z = Z/R_0`.
        From :math:`\psi` the poloidal magnetic field components is evaluated
        from:

        .. math::

            B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z},\\
            B_z &= \frac{1}{R}\frac{\partial\psi}{\partial R},

        and toroidal field is evaluated as

        .. math::

            B_\phi = \frac{B_0 R_0}{R}.

        This input also includes the possibility to have an analytical model
        for toroidal field ripple, which is used if ripple period :math:`N>0`.
        The rippled toroidal field is

        .. math::

            \tilde{B}_\phi = B_\phi( 1 + \delta\cos(N\phi) ),

        where

        .. math::

            \delta = \delta_0 \frac{r}{a}^{\alpha}  e^{-\theta^2},

        and :math:`\theta` is the (geometrical) poloidal angle and :math:`r` is
        the distance from the magnetic axis. The ripple strength,
        :math:`\delta_0`, minor radius, :math:`a`, and ripple penetration
        :math:`\alpha`, can be tuned to adjust the ripple.

        .. [1] A.J. Cerfon, J.P. Freidberg. "One size fits all" analytic
           solutions to the Grad-Shafranov equation. Physics of Plasmas (2010).
           https://doi.org/10.1063/1.3328818
        """
        parameters = variants.parse_parameters(
            rmajor, axisz, axisb, psiscaling, coefficients, axisr, psilimits,
            nripple, rminor, ripplepenetration, ripplescaling,
        )
        default_coefs = np.array([
            2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
            6.200e-03, -1.205e-03, -3.701e-05,  0,
            0,          0,          0,          0,
            -0.155])
        variants.validate_required_parameters(
            parameters,
            names=["rmajor", "axisz", "axisb", "psiscaling", "coefficients"],
            units=["m", "m", "T", "Wb/m", "1"],
            shape=[(1,), (1,), (1,), (1,), (13,)],
            dtype="f8",
            default=[1., 0., 1., 1., default_coefs],
        )

        axisxy = aeq.find_axis(
            parameters["coefficients"],
            )
        psi0 = aeq.psi(
            axisxy[0], axisxy[1], parameters["coefficients"][:-1],
            parameters["coefficients"][-1],
            ) * parameters["psiscaling"]
        psi1  = 0
        raxis = axisxy[0]*parameters["rmajor"]
        padding = 1e-8 * unyt.Wb/unyt.m
        psi0 = psi0 - padding if psi0 < psi1 else psi0 + padding

        variants.validate_optional_parameters(
            parameters,
            ["axisr", "psilimits", "nripple", "rminor", "ripplepenetration",
             "ripplescaling"],
            ("m", "Wb/m", "1", "m", "1", "1"),
            [(1,), (2,), (1,), (1,), (1,), (1,)],
            ["f8", "f8", "i8", "f8", "f8", "f8"],
            [raxis.v, np.array([psi0[0], 0.]), 0, 1., 1., 1.],
        )

        meta = variants.new_metadata("BfieldAnalytical", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
