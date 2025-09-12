"""Analytic tokamak equilibrium.

Solution is based on this article:  http://dx.doi.org/10.1063/1.3328818

Note that thorough this module we assume :math:`x=R/R_0` and :math:`y=z`.
"""
from typing import Optional

import scipy.optimize
import numpy as np
from numpy import log

from a5py.plotting.plotting import openfigureifnoaxes


psi_i = (
    lambda *_: 1.0,
    lambda x, _: x**2,
    lambda x, y: x**2*log(x) - y**2,
    lambda x, y: x**4 - 4*x**2*y**2,
    lambda x, y: 3.0  * x**4 * log(x) - 9.0  * x**2 * y**2
        - 12.0 * x**2 * log(x) * y**2 + 2.0  * y**4,
    lambda x, y: x**6 - 12*x**4*y**2 + 8*x**2*y**4,
    lambda x, y: 8*y**6 - 140*x**2*y**4 - 120*x**2*log(x)*y**4
        + 180*x**4*log(x)*y**2 + 75*x**4*y**2 - 15*x**6*log(x),
    lambda _, y: y,
    lambda x, y: y*x**2,
    lambda x, y: y**3 - 3*y*x**2*log(x),
    lambda x, y: 3*y*x**4 - 4*y**3*x**2,
    lambda x, y: 8*y**5 - 45*y*x**4 - 80*y**3*x**2*log(x) + 60*y*x**4*log(x),
    )
"""The basis functions in Eq. (26) which are defined in Eqs. (8) and (27)."""


psix_i = (
    lambda *_: 0.0,
    lambda x, _: 2*x,
    lambda x, _: 2*x*log(x)+x,
    lambda x, y: 4*x**3 - 8*x*y**2,
    lambda x, y: 12*x**3*log(x) + 3*x**3-30*x*y**2 - 24*x*log(x)*y**2,
    lambda x, y: 6*x**5 - 48*x**3*y**2 + 16*x*y**4,
    lambda x, y: -400*x*y**4 - 240*x*log(x)*y**4 + 720*x**3*log(x)*y**2
        + 480*x**3*y**2 - 90*x**5*log(x) - 15*x**5,
    lambda *_: 0.0,
    lambda x, y: 2*y*x,
    lambda x, y: -6*y*x*log(x) - 3*y*x,
    lambda x, y: 12*y*x**3 - 8*y**3*x,
    lambda x, y: -120*y*x**3 - 160*y**3*x*log(x) - 80*y**3*x
        + 240*y*x**3*log(x),
    )
""":math:`x`-derivate of the basis functions in Eq. (26)."""


psixx_i = (
    lambda *_: 0.0,
    lambda *_: 2.0,
    lambda x, _: 2*log(x) + 3.0,
    lambda x, y: 12*x**2 - 8*y**2,
    lambda x, y: 36*x**2*log(x) + 21*x**2-54*y**2 - 24*log(x)*y**2,
    lambda x, y: 30*x**4 - 144*x**2*y**2 + 16*y**4,
    lambda x, y: -640*y**4 - 240*log(x)*y**4 + 2160*x**2*log(x)*y**2
        + 2160*x**2*y**2 - 450*x**4*log(x) - 165*x**4,
    lambda *_: 0.0,
    lambda _, y: 2*y,
    lambda x, y: -6*y*log(x) - 9*y,
    lambda x, y: 36*y*x**2 - 8*y**3,
    lambda x, y: -120*y*x**2 - 160*y**3*log(x) - 240*y**3 + 720*y*x**2*log(x),
    )
"""Second-order :math:`x`-derivate of the basis functions in Eq. (26)."""


psiy_i = (
    lambda *_: 0.0,
    lambda *_: 0.0,
    lambda _, y: -2*y,
    lambda x, y: -8*x**2*y,
    lambda x, y: -18*x**2*y - 24*x**2*log(x)*y + 8*y**3,
    lambda x, y: -24*x**4*y + 32*x**2*y**3,
    lambda x, y: 48*y**5 - 560*x**2*y**3 - 480*x**2*log(x)*y**3
        + 360*x**4*log(x)*y + 150*x**4*y,
    lambda *_: 1.0,
    lambda x, _: x**2,
    lambda x, y: 3*y**2 - 3*x**2*log(x),
    lambda x, y: 3*x**4 - 12*y**2*x**2,
    lambda x, y: 40*y**4 - 45*x**4 - 240*y**2*x**2*log(x) + 60*x**4*log(x),
    )
""":math:`y`-derivate of the basis functions in Eq. (26)."""


psiyy_i = (
    lambda *_: 0.0,
    lambda *_: 0.0,
    lambda *_: -2.0,
    lambda x, _: -8*x**2,
    lambda x, y: -18*x**2 - 24*x**2*log(x) + 24*y**2,
    lambda x, y: -24*x**4 + 96*x**2*y**2,
    lambda x, y: 240*y**4 - 1680*x**2*y**2 - 1440*x**2*log(x)*y**2
        + 360*x**4*log(x) + 150*x**4,
    lambda *_: 0.0,
    lambda *_: 0.0,
    lambda _, y: 6*y,
    lambda x, y: -24*y*x**2,
    lambda x, y: 160*y**3 - 480*y*x**2*log(x),
    )
"""Second-order :math:`y`-derivate of the basis functions in Eq. (26)."""


psia = lambda x, _, A: 1.0/8*x**4 + A * ( 0.5*x**2*log(x) - 1.0/8*x**4 )
"""The first two terms in Eq. (26) that are separate from the basis functions.
"""


psiax = lambda x, _, A: 0.5*x**3 + A * ( x*log(x) + x/2 - x**3/2 )
"""The :math:`x`-derivative of the first two terms in Eq. (26)."""


psiaxx = lambda x, _, A: 1.5*x**2 + A * ( log(x) - 1.5*x**2 + 1.5 )
"""The second order :math:`x`-derivative of the first two terms in Eq. (26).
"""


psiay = lambda *_: 0.0
"""The :math:`y`-derivative of the first two terms in Eq. (26).
"""


psiayy = lambda *_: 0.0
"""The second order :math:`y`-derivative of the first two terms inEq. (26).
"""


def psi(x, y, c, A):
    """Evaluate the (total) psi function at the given position.

    Parameters
    ----------
    x : float
        The :math:`R/R_0` coordinate.
    y : float
        The :math:`z/R_0` coordinate.
    c : array_like
        The coefficients :math:`c_i`, where :math:`i=0...11`.
    A : float
        The free parameter :math:`A` in the Eq. (26).

    Returns
    -------
    psi : float
        The value of the poloidal flux at the given position.
    """
    val = psia(x, y, A)
    for i in range(12):
        val += c[i] * psi_i[i](x,y)
    return val


def boundaryconditions2coefficients(
    constant_a: float,
    outer_equatorial_point: tuple[float, float],
    inner_equatorial_point: tuple[float, float],
    upper_high_point: tuple[float, float],
    outer_slope: float,
    inner_slope: float,
    curvature_outer_equatorial: float,
    curvature_inner_equatorial: float,
    curvature_upper_high: float,
    lower_x_point: Optional[tuple[float, float]] = None,
    symmetric: bool = True,
    ):
    """Solve the equilibrium coefficients from the given boundary conditions.

    Parameters
    ----------
    constant_a : float
        The free parameter :math:`A` in the paper.
    outer_equatorial_point : float
        :math:`(x/R0, y/R0)`-coordinates of the outer equatorial point at
        the last closed flux surface.
    curvature_outer_equatorial : float
        Curvature at ``outer_equatorial_point``.
    inner_equatorial_point : float
        :math:`(x/R0, y/R0)`-coordinates of the inner equatorial point at
        the last closed flux surface.
    curvature_inner_equatorial : float
        Curvature at ``inner_equatorial_point``.
    upper_high_point : float
        :math:`(x/R0, y/R0)`-coordinates of the highest point at the last closed
        flux surface.
    curvature_upper_high : float
        Curvature at ``upper_high_point``.
    outer_slope : float
        Outer equatorial point slope.
    inner_slope : float
        Inner equatorial point slope.
    lower_x_point : float, optional
        Lower X-point :math:`(x/R0,y/R0)`-coordinates if the equilibrium is not
        a limiter plasma.

        ``None`` for limiter plasma.
    symmetric: bool, optional
        Whether the equilibrium is symmetric (double-null  or limiter plasma).

    Returns
    -------
    c : array_like
        The coefficients :math:`c` in Eq. (26).
    """
    is_limiter_plasma = lower_x_point is None

    if not symmetric and is_limiter_plasma:
        raise ValueError(
            "Limiter plasmas have to be symmetric. Either specify the lower "
            "X-point or make the plasma symmetric."
            )
    elif symmetric and is_limiter_plasma:
        # Limiter plasma
        A, B = np.zeros((7, 7)), np.zeros((7, 1))
        for i in range(7):
            A[0][i] = psi_i[i](*outer_equatorial_point)
            A[1][i] = psi_i[i](*inner_equatorial_point)
            A[2][i] = psi_i[i](*upper_high_point)
            A[3][i] = psix_i[i](*upper_high_point)
            A[4][i] = (
                  psiyy_i[i](*outer_equatorial_point)
                + curvature_outer_equatorial * psix_i[i](*outer_equatorial_point)
                )
            A[5][i] = (
                  psiyy_i[i](*inner_equatorial_point)
                + curvature_inner_equatorial * psix_i[i](*inner_equatorial_point)
                )
            A[6][i] = (
                  psixx_i[i](*upper_high_point)
                + curvature_upper_high * psiy_i[i](*upper_high_point)
                )

        B[0] = psia(*outer_equatorial_point, constant_a)
        B[1] = psia(*inner_equatorial_point, constant_a)
        B[2] = psia(*upper_high_point, constant_a)
        B[3] = psiax(*upper_high_point, constant_a)
        B[4] = (
              psiayy(*outer_equatorial_point, constant_a)
            + curvature_outer_equatorial * psiax(*outer_equatorial_point, constant_a)
            )
        B[5] = (
            psiayy(*inner_equatorial_point, constant_a)
            + curvature_inner_equatorial * psiax(*inner_equatorial_point, constant_a)
            )
        B[6] = (
            psiaxx(*upper_high_point, constant_a)
            + curvature_upper_high * psiay(*upper_high_point, constant_a)
            )
    elif symmetric and not is_limiter_plasma:
        # Double null
        A, B = np.zeros((7, 7)), np.zeros((7, 1))
        for i in range(7):
            A[0][i] = psi_i[i](*outer_equatorial_point)
            A[1][i] = psi_i[i](*inner_equatorial_point)
            A[2][i] = psi_i[i](*lower_x_point)
            A[3][i] = psix_i[i](*lower_x_point)
            A[4][i] = psiy_i[i](*lower_x_point)
            A[5][i] = (
                  psiyy_i[i](*outer_equatorial_point)
                + curvature_outer_equatorial * psix_i[i](*outer_equatorial_point)
                )
            A[6][i] = (
                  psiyy_i[i](*inner_equatorial_point)
                + curvature_inner_equatorial * psix_i[i](*inner_equatorial_point)
                )

        B[0] = psia(*outer_equatorial_point, constant_a)
        B[1] = psia(*inner_equatorial_point, constant_a)
        B[2] = psia(*lower_x_point, constant_a)
        B[3] = psiax(*lower_x_point, constant_a)
        B[4] = psiay(*lower_x_point, constant_a)
        B[5] = (
              psiayy(*outer_equatorial_point, constant_a)
            + curvature_outer_equatorial * psiax(*outer_equatorial_point, constant_a)
            )
        B[6] = (
              psiayy(*inner_equatorial_point, constant_a)
            + curvature_inner_equatorial * psiax(*inner_equatorial_point, constant_a)
            )
    elif not symmetric and not is_limiter_plasma:
        # Bottom divertor
        A, B = np.zeros((12, 12)), np.zeros((12, 1))
        for i in range(12):
            A[0][i] = psi_i[i](*outer_equatorial_point)
            A[1][i] = psi_i[i](*inner_equatorial_point)
            A[2][i] = psi_i[i](*upper_high_point)
            A[3][i] = psi_i[i](*lower_x_point)
            A[4][i] = (
                  outer_slope * psix_i[i](*outer_equatorial_point)
                + psiy_i[i](*outer_equatorial_point)
                )
            A[5][i] = (
                  inner_slope * psix_i[i](*inner_equatorial_point)
                + psiy_i[i](*inner_equatorial_point)
                )
            A[6][i] = psix_i[i](*upper_high_point)
            A[7][i] = psix_i[i](*lower_x_point)
            A[8][i] = psiy_i[i](*lower_x_point)
            A[9][i] = (
                  psiyy_i[i](*outer_equatorial_point)
                + curvature_outer_equatorial * psix_i[i](*outer_equatorial_point)
                )
            A[10][i]= (
                  psiyy_i[i](*inner_equatorial_point)
                + curvature_inner_equatorial * psix_i[i](*inner_equatorial_point)
                )
            A[11][i]= (
                  psixx_i[i](*upper_high_point)
                + curvature_upper_high * psiy_i[i](*upper_high_point)
                )

        B[0] = psia(*outer_equatorial_point, constant_a)
        B[1] = psia(*inner_equatorial_point, constant_a)
        B[2] = psia(*upper_high_point, constant_a)
        B[3] = psia(*lower_x_point, constant_a)
        B[4] = (
              outer_slope * psiax(*outer_equatorial_point, constant_a)
            + psiay(*outer_equatorial_point, constant_a)
            )
        B[5] = (
              outer_slope * psiax(*inner_equatorial_point, constant_a)
            + psiay(*inner_equatorial_point, constant_a)
            )
        B[6] =psiax(*upper_high_point, constant_a)
        B[7] = psiax(*lower_x_point, constant_a)
        B[8] = psiay(*lower_x_point, constant_a)
        B[9] = (
              curvature_outer_equatorial * psiax(*outer_equatorial_point, constant_a)
            + psiayy(*outer_equatorial_point, constant_a)
            )
        B[10] = (
              curvature_inner_equatorial * psiax(*inner_equatorial_point, constant_a)
            + psiayy(*inner_equatorial_point, constant_a)
            )
        B[11] = (
              curvature_upper_high * psiay(*upper_high_point, constant_a)
            + psiaxx(*upper_high_point, constant_a)
            )

    B = np.negative(B)
    c = np.linalg.lstsq(A, B, rcond=None)[0].flatten()
    if symmetric:
        c = np.append(c, [0.0, 0.0, 0.0, 0.0, 0.0])
    return c


def parameters2coefficients(
        A: float,
        epsilon: float,
        kappa: float,
        delta: float,
        xpoint: Optional[tuple[float, float]] = None,
        symmetric: bool = True,
        ):
    """Solve the equilibrium coefficients based on equilibrium parameters.

    Parameters
    ----------
    A : float
        The free parameter in the paper.
    epsilon : float
        Inverse aspect ratio.
    kappa : float
        Elongation.
    delta : float
        Triangularity.
    xpoint : float, optional
        Lower X-point :math:`(x/R0,y/R0)`-coordinates if the equilibrium is not
        a limiter plasma.

        ``None`` for limiter plasma.
    symmetric: bool, optional
        Whether the equilibrium is symmetric (double-null  or limiter plasma).

    Returns
    -------
    c : array_like
        The coefficients :math:`c` in Eq. (26).
    """
    alpha = np.asin(delta)
    return boundaryconditions2coefficients(
        constant_a = A,
        outer_equatorial_point = (1.0 + epsilon, 0.0),
        curvature_outer_equatorial = -(1.0 + alpha)**2 / (epsilon*kappa**2),
        inner_equatorial_point = (1.0 - epsilon, 0.0),
        curvature_inner_equatorial = (1.0 - alpha)**2 / (epsilon*kappa**2),
        upper_high_point = (1.0 - epsilon*delta, kappa*epsilon),
        curvature_upper_high = -kappa / ( epsilon*np.cos(alpha)**2),
        outer_slope = 0.0,
        inner_slope = 0.0,
        lower_x_point = xpoint,
        symmetric=symmetric,
        )


def find_axis(coefficients):
    """Find axis position (i.e. where we have psi minimum).

    Parameters
    ----------
    coefficients : array_like
        The coefficients :math:`c_i`, where :math:`i=0...11`, and :math:`A`.

    Returns
    -------
    r : float
        Axis :math:`R/R_0` coordinate.
    z : float
        Axis :math:`z/R_0` coordinate.
    """
    pfun = lambda x: psi(x[0], x[1], coefficients[:-1], coefficients[-1])
    return scipy.optimize.minimize(
        pfun, np.array([1.0, 0.0]).ravel(), tol=1e-9,
        ).x


@openfigureifnoaxes()
def plot_equilibrium(coefficients, r0, epsilon, kappa, axes=None):
    """Plot the analytical equilibrium.

    Parameters
    ----------
    coefficients : array_like
        The coefficients :math:`c_i`, where :math:`i=0...11`, and :math:`A`.
    r0 : float
        Major axis R coordinate (not the same as the magnetic axis) [m].
    epsilon : float
        Inverse aspect ratio.
    kappa : float
        Elongation.
    """
    rmin, rmax, zmin, zmax = (
        np.amax((1.0 - epsilon - 0.1, 1e-8)),
        1.0 + epsilon + 0.1,
        -1.2*kappa*epsilon - 0.1,
        1.2*kappa*epsilon + 0.1,
    )
    r = np.linspace(rmin, rmax, 200)
    z = np.linspace(zmin, zmax, 200)
    X, Y = np.meshgrid(r, z)
    psival = psi(X, Y, coefficients[:-1], coefficients[-1]).v

    psi_at_separatrix = 0
    axes.pcolormesh(X*r0, Y*r0, psival)
    axes.contour(X*r0, Y*r0, psival, [psi_at_separatrix])
    axes.plot([r0, r0], [zmin*r0, zmax*r0],
              linestyle="--", label="Major axis")
    axes.set_aspect("equal", adjustable="box")
