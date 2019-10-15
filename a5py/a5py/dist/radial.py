"""
Evaluate radial (1D) distributions from 5D or 6D distributions.

File: radial.py
"""
import copy
import numpy as np
import scipy.constants as const

from . import basic as distmod
from . import conversion as distconv

def eval1d(ascotpy, dist, quantity, rhomin, rhomax, nrho, ma=None, qa=None):

    # Find rho volume
    vol = evalrhovol(
        ascotpy, rhomin, rhomax, nrho,
        dist["r_edges"][0], dist["r_edges"][-1], dist["r_edges"].size,
        dist["phi_edges"][0], dist["phi_edges"][-1], dist["phi_edges"].size,
        dist["z_edges"][0], dist["z_edges"][-1], dist["z_edges"].size )

    # Evaluate the requested quantity in 5D
    if quantity == "density":
        dist = distmod.squeeze(dist, vpar=0, vperp=0, time=0, charge=0)

    if quantity == "epowerdeposition":
        emax = 0.5*ma*np.power(np.maximum( dist["vpar_edges"][-1],
                                           dist["vperp_edges"][-1] ), 2) / const.e
        dist = distconv.convert_vpavpe_to_Exi(dist, ma,
                                              E_edges=np.linspace(0, emax,100),
                                              xi_edges=np.linspace(-1,1,50))
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)
        va = np.sqrt( const.e*2*dist["energy"]/ma )

        distPi = copy.deepcopy(dist)
        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):
                    coefs = ascotpy.eval_collcoefs(ma, qa, dist["r"][ir],
                                                   dist["phi"][ip], dist["z"][iz], 0, va)
                    #print(coefs["F"])
                    #print(coefs["F"].shape)

                    distPi["distribution"][ir, ip, iz, :] = -dist["distribution"][ir, ip, iz, :]\
                                                            *( ma*va*( coefs["F"][0,0,:] + 0*coefs["F"][0,1,:] )\
                                                               +      ma*( coefs["Dpara"][0,0,:] + 0*coefs["Dpara"][0,1,:] ))

        dist = distmod.squeeze(distPi, energy=0)
        #dist = distPi

    # Map  (R,phi,z) to rho
    rho = ascotpy.evaluate(dist["r"], np.deg2rad(dist["phi"]), dist["z"], 0,
                           "rho", grid=True)
    rho = np.squeeze(rho)
    rzpvol = (dist["r_edges"][1] - dist["r_edges"][0]) \
             * (dist["phi_edges"][1] - dist["phi_edges"][0]) \
             * (dist["z_edges"][1] - dist["z_edges"][0])

    dist = dist["distribution"] * rzpvol

    idx   = np.logical_and(rho >= rhomin, rho <= rhomax)

    newdist = {}
    newdist["rho_edges"] = np.linspace(rhomin, rhomax, nrho+1)
    newdist["rho"] = newdist["rho_edges"][:-1] \
                     + (newdist["rho_edges"][1] - newdist["rho_edges"][0]) / 2
    newdist["nrho"] = nrho
    newdist["abscissae"] = ["rho"]

    newdist[quantity]  = np.histogram( rho.ravel(), newdist["rho_edges"],
                                       weights=dist.ravel() )[0]
    newdist[quantity] /= vol.ravel()

    return newdist


def evalrhovol(ascotpy, rhomin, rhomax, nrho, rmin, rmax, nr,
               phimin, phimax, nphi, zmin, zmax, nz):

    # Volume of the test space (rectangular toroid)
    V_tot = 0.5 * (rmax**2 - rmin**2) * np.deg2rad(phimax - phimin) * (zmax - zmin)
    drho = (rhomax - rhomin) / nrho

    points_in_bins = np.zeros((nrho, 1))

    n     = 100000
    n_tot = 0
    nmin  = 2

    while np.any(points_in_bins < nmin):
        r   = np.random.uniform(rmin, rmax, n)
        phi = np.random.uniform(phimin, phimax, n)
        z   = np.random.uniform(zmin, zmax, n)

        rho   = ascotpy.evaluate(r, np.deg2rad(phi), z, 0, "rho")
        ind   = np.logical_and(rho < rhomax, rho >= rhomin)
        i_rho = np.floor((rho[ind]-rhomin)/drho).astype(int)

        np.add.at(points_in_bins, (i_rho, 0), 1)

        n_tot = n_tot + n

    return (points_in_bins/n_tot) * V_tot
