"""
Evaluate radial (1D) distributions from 5D or 6D distributions.

File: radial.py
"""
import numpy as np

from . import basic as distmod

def eval1d(ascotpy, dist, quantity, rhomin, rhomax, nrho):

    # Find rho volume
    vol = evalrhovol(
        ascotpy, rhomin, rhomax, nrho,
        dist["r_edges"][0], dist["r_edges"][-1], dist["r_edges"].size,
        dist["phi_edges"][0], dist["phi_edges"][-1], dist["phi_edges"].size,
        dist["z_edges"][0], dist["z_edges"][-1], dist["z_edges"].size )

    # Evaluate the requested quantity in 5D
    if quantity == "density":
        dist = distmod.squeeze(dist, vpar=0, vperp=0, time=0, charge=0)

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
