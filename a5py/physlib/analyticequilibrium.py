"""Analytic tokamak equilibrium.

Solution is based on this article:  http://dx.doi.org/10.1063/1.3328818

Note that the x and y coordinates here does not correspond to R and z directly.
"""
import math as m
import numpy as np
import scipy.optimize
from numpy import log as log
from numpy import power as p

def psi(x,y,i):
    """Evaluate basis functions Eqs. (8) and (27).
    """
    return {
        0 : 1.0,
        1 : p(x,2),
        2 : p(x,2)*log(x) - p(y,2),
        3 : p(x,4) - 4*p(x,2)*p(y,2),
        4 : 3.0  * p(x,4.0) * log(x)
          - 9.0  * p(x,2.0) * p(y,2.0)
          - 12.0 * p(x,2.0) * log(x) * p(y,2.0)
          + 2.0  * p(y,4.0),
        5 : p(x,6) - 12*p(x,4)*p(y,2) + 8*p(x,2)*p(y,4),
        6 : 8*p(y,6) - 140*p(x,2)*p(y,4) - 120*p(x,2)*log(x)*p(y,4)
        + 180*p(x,4)*log(x)*p(y,2) + 75*p(x,4)*p(y,2) - 15*p(x,6)*log(x),
        7 : y,
        8 : y*p(x,2),
        9 : p(y,3)-3*y*p(x,2)*log(x),
        10: 3*y*p(x,4)-4*p(y,3)*p(x,2),
        11: 8*p(y,5)-45*y*p(x,4)-80*p(y,3)*p(x,2)*log(x)+60*y*p(x,4)*log(x),
        }[i]

def psix(x,y,i):
    """Evaluate x-derivative of the basis functions.
    """
    return {
        0 : 0.0,
        1 : 2*x,
        2 : 2*x*log(x)+x,
        3 : 4*x**3-8*x*y**2,
        4 : 12*x**3*log(x)+3*x**3-30*x*y**2-24*x*log(x)*y**2,
        5 : 6*x**5-48*x**3*y**2+16*x*y**4,
        6 : -400*x*y**4-240*x*log(x)*y**4+720*x**3*log(x)*y**2
        +480*x**3*y**2-90*x**5*log(x)-15*x**5,
        7 : 0,
        8 : 2*y*x,
        9 : -6*y*x*log(x)-3*y*x,
        10: 12*y*x**3-8*y**3*x,
        11: -120*y*x**3-160*y**3*x*log(x)-80*y**3*x+240*y*x**3*log(x),
        }[i]


def psixx(x,y,i):
    """Evaluate x^2-derivative of the basis functions.
    """
    return {
        0 : 0,
        1 : 2.0,
        2 : 2*log(x)+3,
        3 : 12*x**2-8*y**2,
        4 : 36*x**2*log(x)+21*x**2-54*y**2-24*log(x)*y**2,
        5 : 30*x**4-144*x**2*y**2+16*y**4,
        6 : -640*y**4-240*log(x)*y**4+2160*x**2*log(x)*y**2+2160*x**2*y**2
        -450*x**4*log(x)-165*x**4,
        7 : 0,
        8 : 2*y,
        9 : -6*y*log(x)-9*y,
        10: 36*y*x**2-8*y**3,
        11: -120*y*x**2-160*y**3*log(x)-240*y**3+720*y*x**2*log(x),
        }[i]

def psiy(x,y,i):
    """Evaluate y-derivative of the basis functions.
    """
    return {
        0 : 0,
        1 : 0,
        2 : -2*y,
        3 : -8*x**2*y,
        4 : -18*x**2*y-24*x**2*log(x)*y+8*y**3,
        5 : -24*x**4*y+32*x**2*y**3,
        6 : 48*y**5-560*x**2*y**3-480*x**2*log(x)*y**3+360*x**4*log(x)*y
        +150*x**4*y,
        7 : 1,
        8 : x**2,
        9 : 3*y**2-3*x**2*log(x),
        10: 3*x**4-12*y**2*x**2,
        11: 40*y**4-45*x**4-240*y**2*x**2*log(x)+60*x**4*log(x) ,
        }[i]

def psiyy(x,y,i):
    """Evaluate y^2-derivative of the basis functions.
    """
    return {
        0 : 0,
        1 : 0,
        2 : -2,
        3 : -8*x**2,
        4 : -18*x**2-24*x**2*log(x)+24*y**2,
        5 : -24*x**4+96*x**2*y**2,
        6 : 240*y**4-1680*x**2*y**2-1440*x**2*log(x)*y**2+360*x**4*log(x)
        +150*x**4,
        7 : 0,
        8 : 0,
        9 : 6*y,
        10: -24*y*x**2,
        11: 160*y**3-480*y*x**2*log(x),
        }[i]

def psipart(x,y,i):
    """This is the (first) part in Eq. (8) that the basis functions don't cover.
    """
    return {
        0 : 1.0/2*p(x,2)*log(x),
        1 : 1.0/8*p(x,4),
        }[i]

def psipartx(x,y,i):
    """Evaluate x-derivative of the first part in Eq. (8).
    """
    return {
        0 : x*log(x)+1.0/2*x,
        1 :  1.0/2*x**3
        }[i]

def psipartxx(x,y,i):
    """Evaluate x^2-derivative of the first part in Eq. (8).
    """
    return {
        0 : log(x)+3.0/2,
        1 : 3.0/2*x**2,
        }[i]

def psiparty(x,y,i):
    """Evaluate y-derivative of the first part in Eq. (8).
    """
    return 0

def psipartyy(x,y,i):
    """Evaluate y^2-derivative of the first part in Eq. (8).
    """
    return 0

def psiX(x, y, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, A):
    """Evaluate the d(total) psi function/dx at the given position.
    """
    return   c0*psix(x,y,0)     \
        +    c1*psix(x,y,1)     \
        +    c2*psix(x,y,2)     \
        +    c3*psix(x,y,3)     \
        +    c4*psix(x,y,4)     \
        +    c5*psix(x,y,5)     \
        +    c6*psix(x,y,6)     \
        +    c7*psix(x,y,7)     \
        +    c8*psix(x,y,8)     \
        +    c9*psix(x,y,9)     \
        +   c10*psix(x,y,10)    \
        +   c11*psix(x,y,11)    \
        +     A*psipartx(x,y,0) \
        + (1-A)*psipartx(x,y,1)

def psiXX(x, y, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, zeta):
    """Evaluate the d^2(total) psi function/dxdx at the given position.
    """
    return c1*2                                                      \
         + c2*(2*log(x) + 3)                                         \
         + c3*(12*x**2 - 8*y**2)                                     \
         + c4*(36*x**2*log(x) + 21*x**2 - 54*y**2 - 24*log(x)*y**2)  \
         + c5*(30*x**4 - 144*x**2*y**2 + 16*y**4)                    \
         + c6*( -640*y**4 - 240*log(x)*y**4 + 2160*x**2*log(x)*y**2  \
               + 2160*x**2*y**2 - 450*x**4*log(x) - 165*x**4)        \
         + cos(zeta)*log(x) + 3.0/2*cos(zeta) + 3.0/2*sin(zeta)*x**2

def psiY(x, y, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, A):
    """Evaluate the d(total) psi function/dy at the given position.
    """
    return     c0*psiy(x,y,0)    \
         +     c1*psiy(x,y,1)    \
         +     c2*psiy(x,y,2)    \
         +     c3*psiy(x,y,3)    \
         +     c4*psiy(x,y,4)    \
         +     c5*psiy(x,y,5)    \
         +     c6*psiy(x,y,6)    \
         +     c7*psiy(x,y,7)    \
         +     c8*psiy(x,y,8)    \
         +     c9*psiy(x,y,9)    \
         +    c10*psiy(x,y,10)   \
         +    c11*psiy(x,y,11)   \
         +     A*psiparty(x,y,0) \
         + (1-A)*psiparty(x,y,1)

def psiYY(x, y, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, zeta):
    """Evaluate the d^2(total) psi function/dydy at the given position.
    """
    return -c2*2                                                     \
          - c3*8*x**2                                                \
          + c4*( -18*x**2 - 24*x**2*log(x) + 24*y**2)                \
          + c5*( -24*x**4 + 96*x**2*y**2)                            \
          + c6*(  240*y**4 - 1680*x**2*y**2 - 1440*x**2*log(x)*y**2
                + 360*x**4*log(x) + 150*x**4)

def psi0(x, y, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, A):
    """Evaluate the (total) psi function at the given position.
    """
    return   c0*psi(x,y,0)     \
        +    c1*psi(x,y,1)     \
        +    c2*psi(x,y,2)     \
        +    c3*psi(x,y,3)     \
        +    c4*psi(x,y,4)     \
        +    c5*psi(x,y,5)     \
        +    c6*psi(x,y,6)     \
        +    c7*psi(x,y,7)     \
        +    c8*psi(x,y,8)     \
        +    c9*psi(x,y,9)     \
        +   c10*psi(x,y,10)    \
        +   c11*psi(x,y,11)    \
        +     A*psipart(x,y,0) \
        + (1-A)*psipart(x,y,1) \

def find_axis(R0, z0, *args):
    """Find axis position (i.e. where we have psi minimum).

    Parameters
    ----------
    R0 : float
        Initial guess for the axis R coordinate [m].
    z0 : float
        Initial guess for the axis z coordinate [m].
    *args
        The coefficients c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,
        and A.

    Returns
    -------
    R : float
        Axis R coordinate [m].
    z : float
        Axis z coordinate [m].
    """
    pfun = lambda x : psi0(x[0], x[1], *args)
    return scipy.optimize.minimize(pfun, np.array([R0/R0, z0/R0]).ravel(),
                                   tol=1e-8).x

def analyticGS(A, R0, epsilon, kappa, delta, Xpointx, Xpointy, sym):
    """Solve the equilibrium coefficients based in simple input.

    This function accepts intuitive tokamak parameters to solve for
    the coefficients.

    Parameters
    ----------
    A : float
        The free parameter in the paper.
    R0 : float
        Major axis R coordinate (not the same as the magnetic axis) [m].
    epsilon : float
        Inverse aspect ratio.
    kappa : float
        Elongation.
    delta : float
        Triangularity.
    Xpointx : float
        Separatrix x-coordinate [m] (if non-symmetric Eq.).
    Xpointy : float
        Separatrix y-coordinate [m] (if non-symmetric Eq.).
    sym : bool
        Make up-down symmetric background.

    Returns
    -------
    c : array_like
        The c coefficients.
    """
    alpha = m.asin(delta)

    C = aGS_allInputs(
        A                    = A,
        outerEqPoint         = [1+epsilon, 0],
        innerEqPoint         = [1-epsilon, 0],
        upperHighpoint       = [1-epsilon*delta, kappa*epsilon],
        lowerXPoint          = [Xpointx, Xpointy],
        outerSlope           = 0,
        innerSlope           = 0,
        curvatureOuterEq     = -(1+alpha)**2/(epsilon*kappa**2),
        curvatureInnerEq     = (1-alpha)**2/(epsilon*kappa**2),
        curvatureHighpointEq = -kappa/(epsilon*(m.cos(alpha))**2),
        sym                  = sym)

    return C

def aGS_allInputs(A, outerEqPoint, innerEqPoint, upperHighpoint, lowerXPoint,
                  outerSlope, innerSlope, curvatureOuterEq, curvatureInnerEq,
                  curvatureHighpointEq, sym):
    """Solve the equilibrium coefficients from the given boundary conditions.

    Parameters
    ----------
    outerEqPoint : float
        xy-coordinates of the outer equatorial point at separatrix.
    innerEqPoint : float
        xy-coordinates of the inner equatorial point at separatrix.
    upperHighPoint : float
        xy-coordinates of the highest separatrix point.
    lowerXpoint : float
        Outer equatorial point slope.
    outerSlope : float
        Outer equatorial point slope.
    innerSlope : float
        Inner equatorial point slope.
    curvatureOuterEq : float
        Curvature at the outboard equatorial plane.
    curvatureInnerEq : float
        Curvature at the inboard equatorial plane.
    curvatureHighPointEq : float
        Curvature at the top point.
    sym : float
        Make up-down symmetric background.

    Returns
    -------
    c : array_like
        The c coefficients.
    """

    # Rename A as Acoeff
    Acoeff = A

    ###########################################################################
    #                                                                         #
    #   Construct the matrix A of the boundary conditions for the funtions    #
    #   which are solutions to the homogeneous equation                       #
    #                                                                         #
    ###########################################################################
    if sym:
        A = np.zeros((7,7))
        for ipsi in range(7):
            A[0][ipsi] = psi(outerEqPoint[0], outerEqPoint[1], ipsi)
            A[1][ipsi] = psi(innerEqPoint[0], innerEqPoint[1], ipsi)
            A[2][ipsi] = psi(upperHighpoint[0], upperHighpoint[1], ipsi)
            A[3][ipsi] = psix(upperHighpoint[0], upperHighpoint[1], ipsi)
            A[4][ipsi] = psiyy(outerEqPoint[0], outerEqPoint[1], ipsi) \
                         + curvatureOuterEq * psix(outerEqPoint[0],
                                                   outerEqPoint[1], ipsi)
            A[5][ipsi] = psiyy(innerEqPoint[0], innerEqPoint[1], ipsi) \
                         + curvatureInnerEq * psix(innerEqPoint[0],
                                                   innerEqPoint[1], ipsi)
            A[6][ipsi] = psixx(upperHighpoint[0], upperHighpoint[1], ipsi) \
                         + curvatureHighpointEq*psiy(upperHighpoint[0],
                                                     upperHighpoint[1], ipsi)
    else:
        A = np.zeros((12,12))
        for ipsi in range(12):
            A[0][ipsi] = psi(outerEqPoint[0], outerEqPoint[1], ipsi)
            A[1][ipsi] = psi(innerEqPoint[0], innerEqPoint[1], ipsi)
            A[2][ipsi] = psi(upperHighpoint[0], upperHighpoint[1], ipsi)
            A[3][ipsi] = psi(lowerXPoint[0], lowerXPoint[1], ipsi)
            A[4][ipsi] = psix(outerEqPoint[0],
                              outerEqPoint[1], ipsi) * outerSlope \
                       + psiy(outerEqPoint[0], outerEqPoint[1], ipsi)
            A[5][ipsi] = psix(innerEqPoint[0],
                              innerEqPoint[1], ipsi) \
                       * innerSlope+psiy(innerEqPoint[0], innerEqPoint[1], ipsi)
            A[6][ipsi] = psix(upperHighpoint[0], upperHighpoint[1], ipsi)
            A[7][ipsi] = psix(lowerXPoint[0], lowerXPoint[1], ipsi)
            A[8][ipsi] = psiy(lowerXPoint[0], lowerXPoint[1], ipsi)
            A[9][ipsi] = psiyy(outerEqPoint[0], outerEqPoint[1], ipsi) \
                         + curvatureOuterEq*psix(outerEqPoint[0],
                                                 outerEqPoint[1], ipsi)
            A[10][ipsi]= psiyy(innerEqPoint[0], innerEqPoint[1], ipsi) \
                         + curvatureInnerEq*psix(innerEqPoint[0],
                                                 innerEqPoint[1], ipsi)
            A[11][ipsi]= psixx(upperHighpoint[0], upperHighpoint[1], ipsi) \
                         + curvatureHighpointEq*psiy(upperHighpoint[0],
                                                     upperHighpoint[1], ipsi)

    ###########################################################################
    #                                                                         #
    #   Construct the matrix B of the boundary conditions for the particular  #
    #   solutions to the equation                                             #
    #                                                                         #
    ###########################################################################

    if sym:
        B = np.zeros((7,1))
        f = np.array([Acoeff, 1-Acoeff])
        for ipart in range(2):
            B[0] = B[0] + f[ipart]*psipart(outerEqPoint[0],
                                           outerEqPoint[1], ipart)
            B[1] = B[1] + f[ipart]*psipart(innerEqPoint[0],
                                           innerEqPoint[1], ipart)
            B[2] = B[2] + f[ipart]*psipart(upperHighpoint[0],
                                           upperHighpoint[1], ipart)
            B[3] = B[3] + f[ipart]*psipartx(upperHighpoint[0],
                                            upperHighpoint[1], ipart)
            B[4] = B[4] + f[ipart] * (
                psipartyy(outerEqPoint[0], outerEqPoint[1], ipart)
                + curvatureOuterEq*psipartx(outerEqPoint[0],
                                            outerEqPoint[1], ipart) )
            B[5] = B[5] + f[ipart]*(
                psipartyy(innerEqPoint[0], innerEqPoint[1], ipart)
                + curvatureInnerEq*psipartx(innerEqPoint[0],
                                            innerEqPoint[1], ipart) )
            B[6] = B[6] + f[ipart]*(
                psipartxx(upperHighpoint[0], upperHighpoint[1], ipart)
                + curvatureHighpointEq*psiparty(upperHighpoint[0],
                                                upperHighpoint[1], ipart) )
    else:
        B = np.zeros((12,1))
        f = np.array([Acoeff, 1-Acoeff])
        for ipart in range(2):
            B[0] = B[0] + f[ipart]*psipart(outerEqPoint[0],
                                           outerEqPoint[1], ipart)
            B[1] = B[1] + f[ipart]*psipart(innerEqPoint[0],
                                           innerEqPoint[1], ipart)
            B[2] = B[2] + f[ipart]*psipart(upperHighpoint[0],
                                           upperHighpoint[1], ipart)
            B[3] = B[3] + f[ipart]*psipart(lowerXPoint[0], lowerXPoint, ipart)
            B[4] = B[4] + f[ipart]*(
                outerSlope*psipartx(outerEqPoint[0], outerEqPoint[1], ipart)
                + psiparty(outerEqPoint[0], outerEqPoint[1], ipart) )
            B[5] = B[5] + f[ipart]*(
                outerSlope*psipartx(innerEqPoint[0], innerEqPoint[1], ipart)
                + psiparty(innerEqPoint[0], innerEqPoint[1], ipart) )
            B[6] = B[6] + f[ipart]*psipartx(upperHighpoint[0],
                                            upperHighpoint[1], ipart)
            B[7] = B[7] + f[ipart]*psipartx(lowerXPoint[0],
                                            lowerXPoint[1], ipart)
            B[8] = B[8] + f[ipart]*psiparty(lowerXPoint[0],
                                            lowerXPoint[1], ipart)
            B[9] = B[9] + f[ipart]*(
                curvatureOuterEq*psipartx(outerEqPoint[0],
                                          outerEqPoint[1], ipart)
                + psipartyy(outerEqPoint[0], outerEqPoint[1], ipart) )
            B[10] = B[10] + f[ipart]*(
                curvatureInnerEq*psipartx(innerEqPoint[0],
                                          innerEqPoint[1], ipart)
                + psipartyy(innerEqPoint[0], innerEqPoint[1], ipart) )
            B[11] = B[11] + f[ipart]*(
                curvatureHighpointEq*psiparty(upperHighpoint[0],
                                              upperHighpoint[1], ipart)
                + psipartxx(upperHighpoint[0], upperHighpoint[1], ipart) )

    B = np.negative(B)
    ###########################################################################
    #                                                                         #
    #   Solve the linear system for the coefficients C of the general         #
    #   solution to the equation                                              #
    #                                                                         #
    ###########################################################################
    C = np.linalg.lstsq(A, B, rcond=None)[0].flatten()

    if sym: C = np.append(C,[0,0,0,0,0])

    return C

def plot_equilibrium(C, A, R0, epsilon, kappa):
    """Plot the analytical equilibrium.

    Parameters
    ----------
    C : float, array_like
        The c coefficients.
    A : float
        The free parameter in the paper.
    R0 : float
        Major axis R coordinate (not the same as the magnetic axis) [m].
    epsilon : float
        Inverse aspect ratio.
    kappa : float
        Elongation.
    """
    import matplotlib.pyplot  as plt

    R = np.linspace((1-epsilon-0.1), 1+epsilon+0.1, 200)
    Z = np.linspace(-1.2*kappa*epsilon-0.1, 1.2*kappa*epsilon+0.1, 200)
    X,Y = np.meshgrid(R,Z)
    X[X < 0.2] = 0

    Z = psi0(X, Y, C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8],
             C[9], C[10], C[11], A)

    plt.pcolormesh(X*R0, Y*R0, Z)
    plt.contour(X*R0, Y*R0, Z, [0])
    plt.plot([R0, R0],
             np.array([-1.5*kappa*epsilon, 1.5*kappa*epsilon])*R0, '--')
    plt.axis('equal')
    plt.show(block=False)

    return C