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

import a5py.ascot5io.boozer as boozer
from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy.ascotpy import Ascotpy


def generate(a5, rgrid, zgrid, npsi, nthgeo, nthbzr, zaxis, raxis, psi0, psi1,
             evalq=True, psipad=0.005):

    # Resolutions for grids used internally
    rres  = 400
    zres  = 800
    thres = 1000

    # Set up the grids
    R = np.linspace(rgrid[0], rgrid[-1], rres) # For contours
    Z = np.linspace(zgrid[0], zgrid[-1], zres) # For contours

    # Evaluating contours right at the axis or separatrix might cause issues
    psimin = np.amin([psi0, psi1])
    psimax = np.amax([psi0, psi1])
    dpsi   = psimax - psimin

    # For contours and thbzr and nu datas
    psigrid   = np.linspace(psimin+psipad*dpsi, psimax-psipad*dpsi, npsi)
    # For thtbzr data
    thgeogrid = np.linspace(0, 2*np.pi, nthgeo)
    # For nu data
    thbzrgrid = np.linspace(0, 2*np.pi, nthbzr)
    # Angles at which contour points are evaluated
    thgrid = np.linspace(0, 2*np.pi, thres)

    ## Set up the data tables (psi can be evaluated directly) ##
    psitable = np.squeeze(a5.evaluate(R=rgrid, phi=0, z=zgrid, t=0,
                                      quantity="psi", grid=True))
    thtable  = np.zeros((psigrid.size,thgeogrid.size))
    nutable  = np.zeros((psigrid.size,thbzrgrid.size))

    # Also store q values to an array which can be returned
    qprof = np.zeros(psigrid.shape)
    Iprof = np.zeros(psigrid.shape)
    gprof = np.zeros(psigrid.shape)

    # Evaluate psi for finding the contours
    psi  = np.squeeze(a5.evaluate(R=R, phi=0, z=Z, t=0,
                                  quantity="psi", grid=True))

    if psimin < psi0:
        psigrid = np.flip(psigrid)

    for i in range(psigrid.size):
        contours = measure.find_contours(psi, psigrid[i])
        contour = np.array([0])
        for j in range(len(contours)):
            if contour.size < contours[j].size:
                contour = contours[j]

        # Normalize coordinates to R,z values
        cr = contour[:, 0]*(R[-1]-R[0])/R.size + R[0]
        cz = contour[:, 1]*(Z[-1]-Z[0])/Z.size + Z[0]

        # We want to have a contour that "begins" at the OMP
        thdat  = np.mod(np.arctan2(cz[:-1]-zaxis, cr[:-1]-raxis), 2*np.pi)
        r = interp1d(thdat, cr[:-1], kind='cubic', fill_value="extrapolate")
        z = interp1d(thdat, cz[:-1], kind='cubic', fill_value="extrapolate")
        cr = r(thgrid)
        cz = z(thgrid)

        # Contour line lengths and differentials
        dL = np.sqrt( (cr[:-1]-cr[1:])**2 + (cz[:-1]-cz[1:])**2 )
        L = np.append(0, np.cumsum(dL))

        # Magnetic field along the contour
        br    = np.squeeze(a5.evaluate(R=cr, phi=0, z=cz, t=0, quantity="br"))
        bphi  = np.squeeze(a5.evaluate(R=cr, phi=0, z=cz, t=0, quantity="bphi"))
        bz    = np.squeeze(a5.evaluate(R=cr, phi=0, z=cz, t=0, quantity="bz"))
        bpol  = np.sqrt(br*br + bz*bz)
        bnorm = np.sqrt(br*br + bphi*bphi + bz*bz)

        ar = cr[1:] - cr[:-1]
        az = cz[1:] - cz[:-1]
        dL = ( ( br[1:] + br[:-1] )*ar + ( bz[1:] + bz[:-1] )*az ) / (bpol[1:] + bpol[:-1])

        # The toroidal current term
        Iprof[i] = -np.sum(
            0.5*dL*( bpol[1:] + bpol[:-1] )
        ) / (2*np.pi)

        # g = R*Bphi
        gprof[i] = cr[0]*bphi[0]

        # Use trapezoidal rule to evaluate the safety factor q(psi)
        qprof[i] = (0.5 / (2*np.pi)) * np.sum(
              dL*(   gprof[i] / ( cr[1:]**2*bpol[1:]   )
                   + gprof[i] / ( cr[:-1]**2*bpol[:-1] ) ) )

        # Evaluate the boozer theta with the trapezoidal rule (set theta(0) = 0)
        jac = (Iprof + qprof*gprof)[i] / bnorm**2
        theta = np.append(0, np.cumsum(
            0.5*dL*(   1/(jac[1:]*bpol[1:]   )
                     + 1/(jac[:-1]*bpol[:-1] ) )
        ) )

        # Normalize theta to interval [0,2pi]. Note that new jacobian is J/s
        s = theta[-1]
        theta = 2*np.pi*theta / s
        th = interp1d(thgrid, theta, "linear")
        thtable[i, :] = th(thgeogrid)

        nu = np.append(0, 0.5*np.cumsum(
            dL*(   1/(cr[:-1]**2*bpol[:-1])
                 + 1/(cr[1:]**2*bpol[1:]) )
        ) )
        nu = interp1d(theta, nu, 'linear', fill_value="extrapolate")
        nutable[i,:] = -(gprof[i]*nu(thbzrgrid) - qprof[i]*thbzrgrid)

    # Flip the data grids to set indices right
    if psimin < psi0:
        thtable = np.flip(thtable,axis=0)
        nutable = np.flip(nutable,axis=0)


    ## Construct the boozer input ##

    # Note that the last contour can be used to define separatrix location
    cr = np.append(cr, cr[0])
    cz = np.append(cz, cz[0])
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
