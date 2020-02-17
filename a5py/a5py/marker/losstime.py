"""
Evaluate transport coefficients from marker loss-time.

The first-passage time is the time at which a marker crosses a given boundary
for the first time. If the transport is uniform on a infinite interval and the
the process obeys the Fokker-Planck equation, i.e., advection-diffusion
equation, the first-passage time distribution is the so-called inverse Gaussian
distribution.

Launching a delta distribution of markers from a same flux surface, and using
ASCOT5 to calculate the times at which they are lost, allows one to estimate
the transport coefficients via curve fitting.

In most cases the transport is not uniform and the coefficients have a radial
dependency. However, if we assume that the transport is uniform *locally* and
the transport is biased outwards (drift is positive), the radially dependent
coefficients can be solved:

1. Initialize several delta populations at different radial positions.
2. Solve the loss time for each population with ASCOT5.
3. Solve the transport coefficients for each population using the aforementioned
   method. These are not yet the final coefficients.
4. Now consider two first passage time distributions for a radial population i.
   The first one is the first-passage time at the edge, i.e., the loss-time
   which is known (L_i). The second one is the first-passage time (P_i) at the
   radial position where the next (going outwards) delta population is located.
   Because the loss time of the population i+1 is also known, L_{i+1}, we can
   write the following relation: <L_i> = <P_i> + <L_{i+1}> where brackets
   indicate mean value.
5. Solving for the distribution P_i, we can use that to fit the first-passage
   time distribution which gives as the (final) transport coefficients at
   position i. Formally, this is an inverse problem where we try to find
   K_i and D_i (drift and diffusion coefficient) that relation
   L_i = int_{-inf}^{inf} P_i(t-tau; K_i, D_i)*L_{i+1}(tau) dtau
   is fulfilled. If the relation above is confusing, google "convolution of
   probability distributions".
6. For practical purposes, it is easier to use the coefficients calculated at
   step 3 to express L_{i+1} analytically rather than to use the raw data.
"""
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.integrate import cumtrapz

def islost(endcond):
    """
    Small helper function to find markers that are lost: edcond=wall or maxrho.
    """
    return np.logical_or(endcond==32, endcond==128)

def firstpassagetime(t, K, D, inirho, endrho):
    """
    Calculate first passage time distribution for given parameters.

    Args:
        t : array_like <br>
            Time coordinate [s].
        K : array_like <br>
            Drift coefficient [rho/s].
        D : array_like <br>
            Diffusion coefficient [rho^2/s].
        inirho : array_like <br>
            Initial radial coordinate [rho].
        endrho : array_like <br>
            Coordinate for which the first-passage time is evaluated [rho].
    Returns:
        The distribution f(t).
    """
    # Helper variables
    drho = endrho - inirho
    c1   = drho / K
    c2   = 2*drho**2/D

    temp = np.zeros(t.shape)
    x    = t[t>0]
    temp[t>0] = np.sqrt( c2 / ( 2*np.pi*x**3 ) ) \
                * np.exp( -c2*( x - c1 )**2 / ( 2*c1**2*x ) )
    return temp


def firstpassagetimecum(t, K, D, inirho, endrho):
    """
    Calculate cumulative first passage time distribution for given parameters.

    Args:
        t : array_like <br>
            Time coordinate [s].
        K : array_like <br>
            Drift coefficient [rho/s].
        D : array_like <br>
            Diffusion coefficient [rho^2/s].
        inirho : array_like <br>
            Initial radial coordinate [rho].
        endrho : array_like <br>
            Coordinate for which the first-passage time is evaluated [rho].
    Returns:
        The distribution F(t).
    """
    # Helper variables
    drho = endrho - inirho
    c1   = drho/K
    c2   = 2*drho**2/D
    t1   =  np.sqrt(c2 / t) * (t/c1 - 1)
    t2   = -np.sqrt(c2 / t) * (t/c1 + 1)
    return 0.5*( 1 + erf( t1 / np.sqrt(2) ) ) \
        + np.exp( 2*c2/c1 ) * 0.5 * ( 1 + erf( t2 / np.sqrt(2)) )


def convolute(t, K, D, inirho, endrho, L2):
    """
    Convolution of two first-passage time distributions.

    The first distribution is evaluated from the given parameters and the second
    one is given as an input. This is because this function is intended for
    inding parameters K and D.

    Args:
        t : array_like <br>
            Time coordinate [s].
        K : array_like <br>
            Drift coefficient of first first-passage time distribution [rho/s].
        D : array_like <br>
            Diffusion coefficient of second first-passage time distribution
            [rho^2/s].
        inirho : array_like <br>
            Initial radial coordinate [rho].
        endrho : array_like <br>
            Coordinate for which the first-passage time is evaluated [rho].
        L2 : array_like <br>
            The second first-passage time distribution on t grid.
    """
    L1 = np.zeros(t.shape)
    for ti in np.arange(t.size):
        P1 = firstpassagetime(t[ti]-t, K, D, inirho, endrho)
        L1[ti] = np.trapz(P1*L2, t)

    temp = cumtrapz(L1, t, initial=0)
    return temp / temp[-1]


def eval_meanvar_coeffs(run=None, inirho=None, losstime=None, endcond=None,
                        endrho=1.0):
    """
    Evaluate coefficients without any fitting just using E[t] and Var[t].

    Args:
        run : RunNode <br>
            Ascot run node container.
        inirho : float, array_like <br>
            Marker initial rho or R (at OMP) coordinate(s).
        losstime : array_like, optional <br>
            Marker final time (including the confined markers).
        endcond : array_like, optional <br>
            Marker end conditions.
        endrho : float, optional <br>
            Radial position for which the passage time is evaluated. Note, use
            units of [m] here if inirho is also in meters.

    Returns:
        Tuple (rho, drift, diff, lost) where rho is the initial rho (or R)
        coordinate, drift is the drift coefficent [rho/s], diff is the diffusion
        coefficient [rho^2/s], and lost is the fraction of lost markers.
    """
    # Init data.
    dataprovided =     inirho   is not None \
                   and losstime is not None \
                   and endcond  is not None
    if run is not None and not dataprovided:
        losstime = run.endstate["time"]
        inirho   = np.mean(run.inistate["rho"])
        endcond  = run.endstate["endcond"]
        lost     = np.sum(islost(endcond)) / losstime.size
        losstime = losstime[islost(endcond)]
    elif dataprovided:
        inirho   = np.mean(inirho)
        lost     = np.sum(islost(endcond)) / losstime.size
        losstime = losstime[islost(endcond)]
    else:
        raise Exception(
            "Please provide either run node or the data explicitly (not both).")

    # Return NaN if no losses.
    if lost == 0:
        return (inirho, np.nan, np.nan, lost)

    # Evaluate and return coefficients
    c1   = np.average(losstime)
    var  = np.average((losstime-c1)**2)
    c2   = c1**3/var
    drho = endrho - inirho

    coeffs = np.array([0,0])
    coeffs[0] = drho / c1
    coeffs[1] = 2*drho**2 / c2

    return (inirho, coeffs[0], coeffs[1], lost)


def eval_radindep_coeffs(run=None, inirho=None, losstime=None, endcond=None,
                         endrho=1.0):
    """
    Evaluates coefficients assuming no radial dependency.

    Provide either the run node (coefficients are fitted using all markers whose
    end condition is maxtime) or the requested data (inirho, losstime, endcond)
    explicitly.

    Args:
        run : RunNode <br>
            Ascot run node container.
        inirho : float, array_like <br>
            Marker initial rho or R (at OMP) coordinate(s).
        losstime : array_like, optional <br>
            Marker final time (including the confined markers).
        endcond : array_like, optional <br>
            Marker end conditions.
        endrho : float, optional <br>
            Radial position for which the passage time is evaluated. Note, use
            units of [m] here if inirho is also in meters.

    Returns:
        Tuple (rho, drift, diff, lost) where rho is the initial rho (or R)
        coordinate, drift is the drift coefficent [rho/s], diff is the diffusion
        coefficient [rho^2/s], and lost is the fraction of lost markers.
    """

    # Init data.
    dataprovided =     inirho   is not None \
                   and losstime is not None \
                   and endcond  is not None
    if run is not None and not dataprovided:
        inirho   = np.mean(run.inistate["rho"])
        endcond  = run.endstate["endcond"]
        lost     = np.sum(islost(endcond)) / losstime.size
        losstime = losstime[islost(endcond)]
        losstime = np.sort(run.endstate["time"])
    elif dataprovided:
        inirho   = np.mean(inirho)
        lost     = np.sum(islost(endcond)) / losstime.size
        losstime = losstime[islost(endcond)]
        losstime = np.sort(losstime)
    else:
        raise Exception(
            "Please provide either run node or the data explicitly (not both).")

    # Return NaN if no losses.
    if lost == 0:
        return (inirho, np.nan, np.nan, lost)

    # Init fitted function.
    def fitfun(x, a, b):
        """
        Wrapper for firstpassagetimecum().
        """
        return firstpassagetimecum(x, a, b, inirho, endrho)

    # Fit and return the coefficients.
    coeffs, _ = curve_fit(fitfun, losstime,
                          np.linspace(0, 1, losstime.size),
                          p0=[10000, 10000], bounds=(0,np.inf))

    return (inirho, coeffs[0], coeffs[1], lost)


def eval_raddep_coeffs(tgrid, rho2, K2, D2, run=None, inirho=None,
                       losstime=None, endcond=None, endrho=1.0):
    """
    Evaluates coefficients assuming a radial dependency.

    Provide either the run node (coefficients are fitted using all markers whose
    end condition is maxtime) or the requested data (inirho, losstime, endcond)
    explicitly.

    Args:
        tgrid : array_like <br>
            Time-grid used to calculate the convolution. Make sure it spans the
            whole regime of possible loss-times and is dense enough. Logarithmic
            grid is likely to provide better results. For example,
            tgrid=10**np.linspace(-6, -2, 1000) has worked well for RE studies.
        rho2 : float <br>
            Radial coordinate (either [rho] or [R]) for the outer point.
        K2 : float <br>
            Drift coefficient for the outer point [rho or m]/[s].
        D2 : float <br>
            Diffusion coefficient for the outer point [rho^2 or m^2]/[s].
        run : RunNode <br>
            Ascot run node container.
        inirho : float, array_like <br>
            Marker initial rho or R (at OMP) coordinate(s).
        losstime : array_like, optional <br>
            Marker final time (including the confined markers).
        endcond : array_like, optional <br>
            Marker end conditions.
        endrho : float, optional <br>
            Radial position for which the passage time is evaluated. Note, use
            units of [m] here if inirho is also in meters.

    Returns:
        Tuple (rho, drift, diff, lost) where rho is the initial rho (or R)
        coordinate, drift is the drift coefficent [rho/s], diff is the diffusion
        coefficient [rho^2/s], and lost is the fraction of lost markers.
    """
    # Init data.
    dataprovided =     inirho   is not None \
                   and losstime is not None \
                   and endcond  is not None
    if run is not None and not dataprovided:
        inirho   = np.mean(run.inistate["rho"])
        endcond  = run.endstate["endcond"]
        lost     = np.sum(islost(endcond)) / losstime.size
        losstime = losstime[islost(endcond)]
        losstime = np.sort(run.endstate["time"])
    elif dataprovided:
        inirho   = np.mean(inirho)
        lost     = np.sum(islost(endcond)) / losstime.size
        losstime = losstime[islost(endcond)]
        losstime = np.sort(losstime)
    else:
        raise Exception(
            "Please provide either run node or the data explicitly (not both).")

    # Return NaN if no losses.
    if lost == 0:
        return (inirho, np.nan, np.nan, lost)

    # Evaluate loss times on a given grid
    losstimeg = np.interp(tgrid, losstime, np.linspace(0, 1, losstime.size),
                         left=0,right=1)
    L2 = firstpassagetime(tgrid, K2, D2, rho2, endrho)

    # Init fitted function.
    def fitfun(x, a, b):
        return convolute(x, a, b, inirho, rho2, L2)

    # Fit and return the coefficients.
    coeffs, _ = curve_fit(fitfun, tgrid, losstimeg,
                          p0=[10000, 10000], bounds=(0,np.inf))
    return (inirho, coeffs[0], coeffs[1], lost)


def eval_coefficients(runs, method="raddep", lastpoint="raddep", gap=1,
                      endrho=1.0, ids=None, tgrid=None, SIunits=False):
    """
    Evaluate coefficients for given radial positions.

    Runs should be sorted so that the outermost run is given last.

    This function assumes that the data for each separate radial location is
    stored in a different run. However, if all data is in a single run, you can
    dublicate that run and use the ids field to control which markers correspond
    to different radial positions.

    The raddep method cannot be used directly to evaluate the last point, which
    means there is a different parameter to control how the last point is
    evaluated if method="raddep".

    Args:
        runs : array_like <br>
            List of runs for different radial points.
        method : str, optional <br>
            Either raddep, radindep, or meanvar.
        lastpoint : str, optional <br>
            What method to use for the last point (raddep means the last point
            is given the same value as the point before that).
        endrho : float, optional <br>
            Coordinate beyond which markers are lost. Give in meters if
            SIunits=True.
        ids : array_like, optional <br>
            ID of the markers for which the coefficients are evaluated.
            Default is all.
        tgrid : array_like, optional <br>
            Time grid needed if the method is raddep. See eval_raddep_coeffs().
        SIunits : bool, optional <br>
            Flag whether coefficients are in units of rho (default) or R.

    Returns:
        Tuple (rho, drift, diff, lost) where rho is the initial rho (or R)
        coordinate, drift is the drift coefficent [rho/s], diff is the diffusion
        coefficient [rho^2/s], and lost is the fraction of lost markers at each
        radial position.
    """

    MEANVAR  = "meanvar"
    RADDEP   = "raddep"
    RADINDEP = "radindep"

    if method not in [MEANVAR, RADDEP, RADINDEP] \
       or lastpoint not in [MEANVAR, RADDEP, RADINDEP]:
        raise Exception(
            "Please provide a valid method and lastpoint" \
            + " (meanvar,raddep,radindep).")

    if method == RADDEP and tgrid is None:
        raise Exception("Please provide tgrid when method=raddep.")

    rho      = np.zeros((len(runs),))
    drift    = np.zeros((len(runs),))
    diff     = np.zeros((len(runs),))
    lostfrac = np.zeros((len(runs),))

    def extractdata(run, ids):
        """
        Small helper for initializing the run data.
        """
        endcond  = run.endstate["endcond"]
        losstime = run.endstate["time"]
        inirho   = run.inistate["rho"]

        if SIunits:
            inirho = run.inistate["r"]

        if ids is not None:
            inirho   = inirho[ids]
            losstime = losstime[ids]
            endcond  = endcond[ids]

        return (inirho, losstime, endcond)

    if method == MEANVAR:
        for i in range(len(runs)):
            (inirho, losstime, endcond) = extractdata(runs[i], ids[i])

            rho[i], drift[i], diff[i], lostfrac[i] = eval_meanvar_coeffs(
                inirho=inirho, losstime=losstime, endcond=endcond,
                endrho=endrho)

        return (rho, drift, diff, lostfrac)

    # Radindep coefficients are needed for the raddep ones so calculate them.
    for i in range(len(runs)):
        (inirho, losstime, endcond) = extractdata(runs[i], ids[i])
        rho[i], drift[i], diff[i], lostfrac[i] = eval_radindep_coeffs(
            inirho=inirho, losstime=losstime, endcond=endcond,
            endrho=endrho)

    if method == RADINDEP:
        return (rho, drift, diff, lostfrac)

    drift0 = np.copy(drift)
    diff0  = np.copy(diff)

    for i in np.arange(len(runs)-1-gap, -1, -1):
        (inirho, losstime, endcond) = extractdata(runs[i], ids[i])
        rho[i], drift[i], diff[i], lostfrac[i] = eval_raddep_coeffs(
            tgrid=tgrid, rho2=rho[i+gap], K2=drift0[i+gap], D2=diff0[i+gap],
            inirho=inirho, losstime=losstime, endcond=endcond,
            endrho=endrho)

    if lastpoint == RADDEP:
        drift[-1] = drift[-2]
        diff[-1]  = diff[-2]
    elif lastpoint == RADINDEP:
        drift[-1] = drift0[-1]
        diff[-1]  = diff0[-1]
    elif lastpoint == MEANVAR:
        (inirho, losstime, endcond) = extractdata(runs[-1], ids[i])
        rho[-1], drift[-1], diff[-1], lostfrac[-1] = eval_meanvar_coeffs(
            inirho=inirho, losstime=losstime, endcond=endcond,
            endrho=endrho)

    return (rho, drift, diff, lostfrac)
