"""
Generate boozer input from a 2D bfield data.

This script takes the psi contours of a 2D equilibrium, and uses
those to construct the corresponding boozer grid. The output
is directly applicable as a boozer input for the ASCOT5.

The script evaluates the q profile in the process.

You need to have B_2DS active before running the script. If you only have
3D data, use write_hdf5_B2DS in ascot5io.B_3DS module to generate 2D
input.

File: generateboozer.py
"""
import sys
import numpy as np

from skimage import measure
from scipy.interpolate import griddata, interp1d
from scipy.integrate import solve_ivp, cumtrapz

import a5py.ascot5io.boozer as boozer
from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy.ascotpy import Ascotpy


def generate(a5, rgrid, zgrid, npsi, nthgeo, nthbzr, zaxis, raxis, psi0, psi1,
             evalq=True, tstep=1e-2, nint=8000):

    # Resolutions for grids used internally
    rres  = 400
    zres  = 800

    # Set up the grids
    R = np.linspace(rgrid[0], rgrid[-1], rres) # For contours
    Z = np.linspace(zgrid[0], zgrid[-1], zres) # For contours

    # Evaluating contours right at the axis or separatrix might cause issues
    psimin = np.amin([psi0, psi1])
    psimax = np.amax([psi0, psi1])
    dpsi   = psimax - psimin
    psipad = 0.01
    psimin += psipad*dpsi
    psimax -= psipad*dpsi

    # For contours and thbzr and nu datas
    psigrid   = np.linspace(psimin, psimax, npsi)
    # For thtbzr data
    thgeogrid = np.linspace(0, 2*np.pi, nthgeo)
    # For nu data
    thbzrgrid = np.linspace(0, 2*np.pi, nthbzr+1)[:-1]
    # Angles at which contour points are evaluated
    thgrid = np.linspace(0, 2*np.pi, nint)

    ## Set up the data tables (psi can be evaluated directly) ##
    psitable = np.squeeze(a5.evaluate(R=rgrid, phi=0, z=zgrid, t=0,
                                      quantity="psi", grid=True))
    thtable  = np.zeros((psigrid.size,thgeogrid.size))
    nutable  = np.zeros((psigrid.size,thbzrgrid.size))

    # Also store q values to an array which can be returned
    qprof = np.zeros(psigrid.shape)
    Iprof = np.zeros(psigrid.shape)
    gprof = np.zeros(psigrid.shape)

    if psimin < psi0:
        psigrid = np.flip(psigrid)

    # Functions needed to trace psi = const. surfaces.
    def event(t, y):
        """Detect the midplane crossing.
        """
        return zaxis - y[1]

    # Terminate when midplane crossing is detected (see docmentation on
    # solve_ivp at scipy) and at right diection. delta is to ensure that
    # we ignore the first crossing by starting the simulation slightly off
    # plane.
    event.terminal  = True
    event.direction = 1 - 2 * (psi1 > psi0)
    delta = -1e-8 * event.direction

    def tracepsi(t, y):
        """Calculate new position when tracing Bpol.
        """
        br = a5.evaluate(R=y[0], phi=0, z=y[1], t=0, quantity="br")
        bz = a5.evaluate(R=y[0], phi=0, z=y[1], t=0, quantity="bz")
        return np.array([br / np.sqrt(br**2 + bz**2),
                         bz / np.sqrt(br**2 + bz**2)]).ravel()

    for i in range(psigrid.size):
        # Find (R,z) location at the outer midplane where psi(R,z) = psi(i)
        rzomp = a5.get_rhotheta_rz(
            np.sqrt((psigrid[i]-psi0) / (psi1 - psi0)), 0, 0, 0)

        # Solve the contour by tracing Bpol. The integration upper limit
        # (2 pi R) is set high as the integration is actually terminated
        # when OMP is crossed (this is taken care by the event. z0 + delta
        # ensures that we don't accidentally terminate it when we start
        # the integration.
        sol = solve_ivp(tracepsi, [0, 2*np.pi*raxis],
                        np.array([rzomp[0], rzomp[1]+delta[0]]).ravel(),
                        max_step=tstep, events=event)
        r = sol.y[0,:]
        z = sol.y[1,:]

        # Interpolate the contour points on the fixed (geometrical) theta
        # grid (at OMP we set thetageom=thetabzr=0)
        theta = np.mod(np.arctan2(z - zaxis, r - raxis), 2*np.pi)
        r = np.interp(thgrid, theta, r, period=2*np.pi)
        z = np.interp(thgrid, theta, z, period=2*np.pi)

        # Magnetic field along the contour (psi can be used to check that
        # the contour was set properly). Drop the last element in r and z
        # as it is the same as first.
        br   = a5.evaluate(R=r[:-1], phi=0, z=z[:-1], t=0, quantity="br")
        bphi = a5.evaluate(R=r[:-1], phi=0, z=z[:-1], t=0, quantity="bphi")
        bz   = a5.evaluate(R=r[:-1], phi=0, z=z[:-1], t=0, quantity="bz")

        bpol  = np.sqrt(br**2 + bz**2)
        bnorm = np.sqrt(br**2 + bphi**2 + bz**2)
        ds    = (np.diff(r) * br + np.diff(z) * bz) / bpol # darc dot e^_pol
        r = r[:-1]; z = z[:-1]

        # The toroidal current term (multiplying this with mu0/2 pi gets
        # enclosed toroidal current)
        Iprof[i] = np.sum( ds * bpol ) / ( 2*np.pi )

        # g = R*Bphi, since Bphi ~ 1/R this is a constant
        gprof[i] = r[0] * bphi[0]

        # The (global) safety factor q(psi)
        qprof[i] = np.sum( ds * gprof[i] / ( r**2 * bpol ) ) / ( 2*np.pi )

        # Boozer coordinate Jacobian is (I - qg) / B^2. Setting it fixes the
        # Boozer poloidal angle which we can now solve.
        jac = (Iprof + qprof*gprof)[i] / bnorm**2
        btheta = np.append(0, np.cumsum( ds / ( jac * bpol ) ))

        # The above Jacobian is for a periodical theta, so theta[-1] should
        # equal to 2 pi already, but normalize it to remove numerical error
        # (note that the new Jacobian would be J / a)
        a = 2*np.pi / btheta[-1]
        btheta *= a
        thtable[i, :] = interp1d(thgrid, btheta, "linear")(thgeogrid)

        # For Boozer toroidal coordinate, we need to integrate the local
        # safety factor along the contour
        nu = gprof[i] * np.append(0, np.cumsum( ds / ( r**2 * bpol ) ) )

        # Interpolate nu used in zeta = phi + nu(psi, theta)
        nutable[i,:] = -interp1d(btheta, nu, 'linear')(thbzrgrid) \
            + qprof[i] * thbzrgrid

    # Flip the data grids to set indices right
    if psimin < psi0:
        thtable = np.flip(thtable,axis=0)
        nutable = np.flip(nutable,axis=0)


    ## Construct the boozer input ##

    # Note that the last contour can be used to define separatrix location
    cr = np.append(r, r[0])
    cz = np.append(z, z[0])
    booz = {
        "psimin":psimin,
        "psimax":psimax,
        "npsi":int(psigrid.size),
        "ntheta":int(thbzrgrid.size),
        "nthetag":int(thgeogrid.size),
        "rmin":rgrid[0],
        "rmax":rgrid[-1],
        "nr":int(rgrid.size),
        "zmin":zgrid[0],
        "zmax":zgrid[-1],
        "nz":int(zgrid.size),
        "r0":raxis,
        "z0":zaxis,
        "psi0":psi0,
        "psi1":psi1,
        "psi_rz":psitable,
        "theta_psithetageom":thtable,
        "nu_psitheta":nutable,
        "nrzs":int(cr.size),
        "rs":cr,
        "zs":cz}

    if evalq:
        rhogrid = np.sqrt( (psigrid - psi0) / (psi1 - psi0) )
        return booz, rhogrid, qprof, gprof, Iprof
    else:
        return booz

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Please give filename of the HDF5 file as an input. Aborting...")
        exit()

    fn = sys.argv[1]

    h5    = Ascot(fn)
    a5    = Ascotpy(fn)
    data  = h5.bfield.active.read()
    rgrid = np.linspace(data["rmin"], data["rmax"], int(data["nr"]))
    zgrid = np.linspace(data["zmin"], data["zmax"], int(data["nz"]))

    a5.init(bfield=True)
    booz, rhoprof, qprof, gprof, Iprof = generate(a5,
                                                  rgrid  = rgrid,
                                                  zgrid  = zgrid,
                                                  npsi   = 100,
                                                  nthgeo = 180,
                                                  nthbzr = 180,
                                                  raxis  = data["axisr"],
                                                  zaxis  = data["axisz"],
                                                  psi0   = data["psi0"],
                                                  psi1   = data["psi1"])
    a5.free(bfield=True)

    print("rho and q(rho):\n")
    print(repr( rhoprof.ravel() ))
    print(repr( qprof.ravel() ))

    boozer.write_hdf5(fn, desc="BOOZER", **booz)
    print("Boozerdata written.")
