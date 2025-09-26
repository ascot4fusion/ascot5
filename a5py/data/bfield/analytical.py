"""Defines analytical magnetic field input class and the corresponding factory
method.
"""
import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.physlib import aeq
from a5py.libascot import LIBASCOT, DataStruct, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods, too-many-instance-attributes
class Struct(DataStruct):
    """Python wrapper for the struct in B_GS.h."""

    _fields_ = [
        ("R0", ctypes.c_double),
        ("z0", ctypes.c_double),
        ("raxis", ctypes.c_double),
        ("zaxis", ctypes.c_double),
        ("B_phi0", ctypes.c_double),
        ("psi0", ctypes.c_double),
        ("psi1", ctypes.c_double),
        ("psi_mult", ctypes.c_double),
        ("psi_coeff", ctypes.c_double * 13),
        ("Nripple", ctypes.c_int32),
        ("a0", ctypes.c_double),
        ("alpha0", ctypes.c_double),
        ("delta0", ctypes.c_double),
        ]


init_fun(
    "B_GS_init",
    ctypes.POINTER(Struct),
    *(8*[ctypes.c_double]),
    ndpointer(ctypes.c_double),
    ctypes.c_int32,
    *(3*[ctypes.c_double]),
    restype=ctypes.c_int32,
    )

init_fun("B_GS_free", ctypes.POINTER(Struct))


@Leaf.register
class BfieldAnalytical(InputVariant):
    """Analytical tokamak field which can be either axisymmetric or 3D."""

    _cdata: Optional[Struct]

    @property
    def rmajor(self) -> unyt.unyt_array:
        """Major radius of the tokamak."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("R0", (), "m")
        assert self._file is not None
        return self._file.read("rmajor")

    @property
    def axisrz(self) -> unyt.unyt_array:
        r"""Magnetic axis :math:`R` and :math:`z` coordinates."""
        if self._cdata is not None:
            r = self._cdata.readonly_carray("raxis", (), "m")
            z = self._cdata.readonly_carray("zaxis", (), "m")
            return unyt.unyt_array([r, z])
        assert self._file is not None
        return self._file.read("axisrz")

    @property
    def axisb(self) -> unyt.unyt_array:
        """Toroidal field strength at the major radius location."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("B_phi0", (), "T")
        assert self._file is not None
        return self._file.read("axisb")

    @property
    def psiscaling(self) -> unyt.unyt_array:
        """Scaling factor for the poloidal flux."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("psi_mult", (), "Wb/rad")
        assert self._file is not None
        return self._file.read("psiscaling")

    @property
    def coefficients(self) -> np.ndarray:
        """Coefficients defining psi: [c0, c1, ..., c11, A]."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("psi_coeff", (13,))
        assert self._file is not None
        return self._file.read("coefficients")

    @property
    def psilimits(self) -> unyt.unyt_array:
        """Poloidal flux values on the magnetic axis and on the separatrix."""
        if self._cdata is not None:
            return unyt.unyt_array((
                self._cdata.readonly_carray("psi0", (), "Wb/rad"),
                self._cdata.readonly_carray("psi1", (), "Wb/rad")
            ))
        assert self._file is not None
        return self._file.read("psilimits")

    @property
    def nripple(self) -> int:
        """Number of TF coils."""
        if self._cdata is not None:
            return int(self._cdata.readonly_carray("Nripple", ()))
        assert self._file is not None
        return int(self._file.read("nripple"))

    @property
    def rminor(self) -> unyt.unyt_array:
        """Minor radius."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("a0", (), "m")
        assert self._file is not None
        return self._file.read("rminor")

    @property
    def ripplepenetration(self) -> unyt.unyt_array:
        """Ripple penetration."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("alpha0", (), "m")
        assert self._file is not None
        return self._file.read("ripplepenetration")

    @property
    def ripplescaling(self) -> unyt.unyt_array:
        """Ripple scaling parameter."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("delta0", (), "m")
        assert self._file is not None
        return self._file.read("ripplescaling")

    #pylint: disable=too-many-arguments
    def _stage(
            self, rmajor: unyt.unyt_array,
            axisrz: unyt.unyt_array,
            axisb: unyt.unyt_array,
            psilimits: unyt.unyt_array,
            coefficients: np.ndarray,
            nripple: int,
            rminor: unyt.unyt_array,
            ripplepenetration: unyt.unyt_array,
            ripplescaling: unyt.unyt_array,
            psiscaling: unyt.unyt_array,
            ) -> None:
        self._cdata = Struct()
        if LIBASCOT.B_GS_init(
            ctypes.byref(self._cdata), rmajor.v, axisrz[1].v, axisrz[0].v,
            axisrz[1].v, axisb.v, psilimits[0].v, psilimits[1].v, psiscaling.v,
            coefficients, nripple, rminor.v, ripplepenetration.v, ripplescaling.v,
            ):
            self._cdata = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "rmajor", "axisrz", "axisb", "psiscaling", "coefficients",
            "psilimits", "rminor", "ripplepenetration", "ripplescaling",
            ]:
            self._file.write(field, getattr(self, field))
        self._file.write(
            "nripple", np.array(self.nripple, dtype="i4")
            )

    def export(self) -> dict[str, unyt.unyt_array | np.ndarray | int]:
        fields = [
            "rmajor", "axisrz", "axisb", "psiscaling", "coefficients",
            "psilimits", "nripple", "rminor", "ripplepenetration",
            "ripplescaling",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            rmajor=self.rmajor, axisrz=self.axisrz, axisb=self.axisb,
            psilimits=self.psilimits, coefficients=self.coefficients,
            nripple=self.nripple, rminor=self.rminor,
            ripplepenetration=self.ripplepenetration,
            ripplescaling=self.ripplescaling, psiscaling=self.psiscaling,
            )

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.B_GS_free(ctypes.byref(self._cdata))
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_bfieldanalytical(
            self,
            rmajor: utils.ArrayLike,
            axisb: utils.ArrayLike,
            psiscaling: utils.ArrayLike,
            coefficients: utils.ArrayLike,
            psilimits: Optional[utils.ArrayLike]=None,
            axisrz: Optional[utils.ArrayLike]=None,
            nripple: Optional[int]=None,
            rminor: Optional[utils.ArrayLike]=None,
            ripplepenetration: Optional[utils.ArrayLike]=None,
            ripplescaling: Optional[utils.ArrayLike]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> BfieldAnalytical:
        r"""Create an analytic tokamak field.

        Parameters
        ----------
        rmajor : float
            Major radius of the tokamak.

            Note that this is not the same as the magnetic axis R coordinate.
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
        axisrz : float, *optional*
            :math:`R` and :math:`z` coordinates of the magnetic axis.

            If not provided, this :math:`R` will be interpolated internally and
            :math:`z=0`.
        psilimits : tuple[float, float], *optional*
            Poloidal flux values on the magnetic axis and on the separatrix.
        nripple : float, *optional*
            Number of TF coils.

            If this value is set and it is larger than zero, the resulting field
            will contain the toroidal field ripple.
        rminor : float, *optional*
            Minor radius.

            Used only for computing ripple strength.
        ripplepenetration : float, *optional*
            Ripple penetration.

            Increasing this parameter decreases ripple strength as
            :math:`\propto\frac{r}{a}^{\alpha}`, where :math:`a` is the minor
            radius and :math:`\alpha` is this parameter.
        ripplescaling : float, *optional*
            Ripple scaling parameter.

            Corresponds to the strength of the ripple at the separatrix.
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
        data : :class:`.BfieldAnalytical`
            Input variant created from the given parameters.

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
        with utils.validate_variables() as v:
            axisb = v.validate("axisb", axisb, (), "T")
            rmajor = v.validate("rmajor", rmajor, (), "m")
            psiscaling = v.validate("psiscaling", psiscaling, (), "Wb/rad")
            coefficients = v.validate("coefficients", coefficients, (13,))

        axisxy = aeq.find_axis(coefficients)
        axisr = axisxy[0] * rmajor.v
        padding = 1e-8 * unyt.Wb/unyt.rad
        psi0 = aeq.psi(axisxy[0], axisxy[1], coefficients[:-1], coefficients[-1]) * psiscaling
        psi1 = 0
        psi0 = psi0 - padding if psi0 < psi1 else psi0 + padding
        with utils.validate_variables() as v:
            axisrz = v.validate("axisrz", axisrz, (2,), "m", default=(axisr, 0.))
            rminor = v.validate("rminor", rminor, (), "m", default=1.0)
            nripple = v.validate("nripple", nripple, (), dtype="i4", default=18)
            psilimits = v.validate("psilimits", psilimits, (2,), "Wb/rad",
                                   default=[psi0.v, 0.0])
            ripplepenetration = v.validate(
                "ripplepenetration", ripplepenetration, (), "1", default=0.01
                )
            ripplescaling = v.validate(
                "ripplescaling", ripplescaling, (), "1", default=1.0
                )

        assert nripple is not None
        leaf = BfieldAnalytical(note=note)
        leaf._stage(
            rmajor=rmajor, axisrz=axisrz, axisb=axisb, psilimits=psilimits,
            coefficients=coefficients, nripple=nripple, rminor=rminor,
            ripplepenetration=ripplepenetration, ripplescaling=ripplescaling,
            psiscaling=psiscaling,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="bfield",
            )
        return leaf
