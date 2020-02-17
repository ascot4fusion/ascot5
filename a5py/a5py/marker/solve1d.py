"""
Solve 1D advection diffusion equation with Monte Carlo method.

The advection-diffusion equation has a reflecting boundary condition at rho0 and
absorbing BC at rho1. The coefficients are stationary but may have spatial
dependency.
"""

import numpy as np
import matplotlib.pyplot as plt

def solve_model(xmin, xmax, simtime, timestep, numberofmarkers,
                driftfun, difffun, distfun, conffun, xgrid=None, tgrid=None):
    """
    Solve 1D advection diffusion equation with boundary conditions.

    Reflective boundary condition is at xmin and absorbing boundary is at xmax,
    where xmin < xmax and x is the radial coordinate.

    Args:
        xmin          : float <br>
            Minimum rho and where the reflecting boundary is located.
        xmax          : float <br>
            Maximum rho and where the absorbing boundary is located.
        simtime         : float <br>
            Simulation time [s].
        timestep        : float <br>
            Simulation time-step [s].
        numberofmarkers : float <br>
            Number of markers. Increase to reduce Monte Carlo noise.
        driftfun        : function <br>
            Function that evaluates drift at the given radial position.
        difffun         : function <br>
            Function that evaluates diffusion at the given radial position.
        distfun         : function <br>
            Function that evaluates initial marker distribution at the given
            radial position.
        conffun         : function <br>
            Function that evaluates which markers are born stationary.
        xgrid           : array_like, optional <br>
            Spatial grid for the density histogram.
        tgrid           : array_like <br>
            Temporal grid for the density histogram.

    Returns:
        Dictionary storing the final x coordinates for the confined and
        un-confined markers alike, losstimes, and the density histogram.
    """
    Nmrk      = numberofmarkers
    xcoords   = np.array([])
    xconfined = np.array([])
    losstime  = simtime*np.ones((Nmrk,))

    if xgrid is None:
        xgrid = np.linspace(xmin, xmax, 10)
    if tgrid is None:
        tgrid = np.linspace(timestep0, simtime, 1000)
    density = np.zeros( (xgrid.size-1, tgrid.size-1) )

    # Evaluate initial profile
    N = 0
    while N < Nmrk:
        # Generate uniform distribution of coordinates and a random number
        # X c [0,1] for each coordinate. Accept coordinates for which
        # X(x) < P(x), where P is the given distribution function.
        x = xmin + (xmax - xmin)*np.random.rand(Nmrk - N,)
        accept  = distfun(x) - np.random.rand(Nmrk - N,)
        x       = x[accept>0]

        # Each coordinate has the probability of C(x) of being confined where
        # C is the given fraction of confined markers.
        conf      = conffun(x) - np.random.rand(x.size,)
        xcoords   = np.append(xcoords, x[conf<=0])
        xconfined = np.append(xconfined, x[conf>0])

        N = xcoords.size + xconfined.size

    # Set confined markers to correct bins
    weights = np.ones(xconfined.shape) * (tgrid[1]-tgrid[0])
    for ti in range(tgrid.size-1):
        density[:,ti] += np.histogram(xconfined, bins=xgrid, weights=weights)[0]

    # The simulation loop
    ismoving = np.ones((Nmrk,), dtype="i8")
    nsteps   = int(simtime/timestep)
    for i in range(nsteps):

        # Init step
        Nmoving  = np.sum(ismoving)
        x0 = xcoords[ismoving]
        K  = driftfun(x0)
        D  = difffun(x0)
        w  = np.random.randn(Nmoving,)

        # Take the step and apply boundary conditions
        x1 = x0 + K*timestep + np.sqrt(2*D*timestep)*w

        reflect = x1 < xmin
        absorb  = x1 >= xmax

        x1[reflect] = 2*xmin - x1[reflect]
        losstime[np.argwhere(ismoving)[absorb]] = (i+1)*timestep

        # Update coordinates and the histogram
        xcoords[ismoving] = x1
        ismoving = xcoords <= xmax

        ti = np.floor( i*timestep / ( tgrid[1] - tgrid[0] ) ).astype("i8")
        weights = np.ones(x1.shape) * timestep
        density[:,ti] += np.histogram(x1, bins=xgrid, weights=weights)[0]

    density /= simtime
    losstime = np.sort(losstime[losstime<simtime])

    out = {"xcoord" : xcoords, "losstime" : losstime, "xconfined" : xconfined,
           "xgrid" : xgrid, "tgrid" : tgrid, "density" : density}
    return out


def solve_model_lincoeff(xmin, xmax, simtime, timestep, numberofmarkers,
                         datagrid, driftdata, diffdata, distdata, confdata,
                         xgrid=None, tgrid=None):
    """
    Solve the 1D model using linear fit on the data.

    If xgrid does not span the whole interval [xmin, xmax], the distribution is
    set to zero outside the grid while the closest value (i.e. first or last) is
    used for the transport coefficients.
    """

    def evaldrift(x):
        """
        Evaluate advection with linear interpolation.
        """
        return np.interp(x, datagrid, driftdata,
                         left=driftdata[0], right=driftdata[-1])

    def evaldiff(x):
        """
        Evaluate diffusion with linear interpolation.
        """
        return np.interp(x, datagrid, diffdata,
                         left=diffdata[0], right=diffdata[-1])

    def evaldist(x):
        """
        Evaluate initial distribution with linear interpolation.
        """
        return np.interp(x, datagrid, distdata, left=0, right=0)

    def evalconf(x):
        """
        Evaluate confined fraction with linear interpolation.
        """
        return np.interp(x, datagrid, confdata, left=0, right=0)

    return solve_model(xmin, xmax, simtime, timestep, numberofmarkers,
                       evaldrift, evaldiff, evaldist, evalconf,
                       xgrid=xgrid, tgrid=tgrid)
