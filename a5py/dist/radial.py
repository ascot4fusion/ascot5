"""
Evaluate radial (1D) distributions from 5D or 6D distributions.

File: radial.py
"""
import copy
import numpy as np

from . import basic as distmod
from . import conversion as distconv
import unyt

def eval1d(ascotpy, dist, quantity, rhomin, rhomax, nrho, ma=None, qa=None, vol=None, area=None):
    '''
    Args units:
            ma : float <br>
                Test particle mass [kg].
            qa : float <br>
                Test particle charge [C].
    '''
    
    #if ma is not None:
    #    print('Marker mass {} kg'.format(ma))
    #if qa is not None:
    #    print('Marker charge {} C'.format(qa))
        
    
    # working copy
    dist = copy.deepcopy(dist)

    if vol is None or area is None:
        # Find rho volume
        vol,area = evalrhovol(
            ascotpy, rhomin, rhomax, nrho,
            dist["r_edges"][0], dist["r_edges"][-1], dist["r_edges"].size,
            dist["phi_edges"][0], dist["phi_edges"][-1], dist["phi_edges"].size,
            dist["z_edges"][0], dist["z_edges"][-1], dist["z_edges"].size )

    dist = eval_quantity_5d(ascotpy, dist, quantity, ma, qa)
    
    # Map  (R,phi,z) to rho
    rho = ascotpy.evaluate(dist["r"], np.deg2rad(dist["phi"]), dist["z"], 0,
                           "rho", grid=True)
    
    rho = np.squeeze(rho)
    rzpvol = (dist["r_edges"][1] - dist["r_edges"][0]) \
             * (dist["phi_edges"][1] - dist["phi_edges"][0]) \
             * (dist["z_edges"][1] - dist["z_edges"][0])

    dist = dist["distribution"] * rzpvol

    #idx   = np.logical_and(rho >= rhomin, rho <= rhomax)

    newdist = {}
    newdist["rho_edges"] = np.linspace(rhomin, rhomax, nrho+1)
    newdist["rho"] = newdist["rho_edges"][:-1] \
                     + (newdist["rho_edges"][1] - newdist["rho_edges"][0]) / 2
    newdist["nrho"] = nrho
    newdist["abscissae"] = ["rho"]
    newdist["vol"] = vol.ravel()
    newdist["area"] = area.ravel()

    newdist[quantity]  = np.histogram( rho.ravel(), newdist["rho_edges"],
                                       weights=dist.ravel() )[0]
    newdist[quantity] /= newdist["vol"]
    
    return newdist

def eval1d_fromRhodist(ascotpy, dist, quantity, rmin, rmax, zmin, zmax, ma=None, qa=None, vol=None, area=None):
    '''
    Args units:
            ma : float <br>
                Test particle mass [kg].
            qa : float <br>
                Test particle charge [C].
    '''
    
    #if ma is not None:
    #    print('Marker mass {} kg'.format(ma))
    #if qa is not None:
    #    print('Marker charge {} C'.format(qa))
        


    
    # working copy
    dist = copy.deepcopy(dist)

    if vol is None or area is None:
        # Find rho volume
        vol,area = evalrhovol(
            ascotpy, 
            dist["rho_edges"][0], dist["rho_edges"][-1], dist["rho_edges"].size-1,
            rmin, rmax, None,
            dist["phi_edges"][0], dist["phi_edges"][-1], None,
            zmin, zmax, None )

    dist = eval_quantity_5d_rho(ascotpy, dist, quantity, ma, qa)
    distmod.squeeze(dist, theta=0, phi=0)

    dist["vol"] = vol.ravel()
    dist["area"] = area.ravel()

    # Move from density per rho-slot to density per m3
    dist[quantity] = dist['distribution'] * (dist["rho_edges"][1] - dist["rho_edges"][0]) /  dist["vol"]
    del( dist['distribution'] )

    return dist

def eval_quantity_5d_rho(ascotpy, dist, quantity, ma, qa):
    # Evaluate the requested quantity in 5D
    if quantity == "density":
        dist = distmod.squeeze(dist, ppar=0, pperp=0, time=0, charge=0)

    elif quantity == "energydensity":
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)

        distE = copy.deepcopy(dist)
        for iphi in range(dist["phi"].size):
            for itheta in range(dist["theta"].size):
                for irho in range(dist["rho"].size):
                    distE["distribution"][irho, itheta, iphi, :] = dist["distribution"][irho, itheta, iphi,:] * dist["energy"] * unyt.e

        dist = distmod.squeeze(distE, energy=0)

    elif quantity == "pressure":
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)

        distp = copy.deepcopy(dist)
        for iphi in range(dist["phi"].size):
            for itheta in range(dist["theta"].size):
                for irho in range(dist["rho"].size):
                    velocity2 = 2 * dist["energy"] * unyt.e / ma
                    distp["distribution"][irho, itheta, iphi, :] = dist["distribution"][irho, itheta, iphi,:] * ma * velocity2 / 3.0

        dist = distmod.squeeze(distp, energy=0)

    elif quantity == "toroidalcurrent":
        dist = distmod.squeeze(dist, pperp=0, time=0, charge=0)

        distj = copy.deepcopy(dist)

        for iphi in range(dist["phi"].size):
            for itheta in range(dist["theta"].size):
                                #Calculate R,phi,z for all rho-values for single phi,theta
                (r, z) = ascotpy.get_rhotheta_rz(dist["rho"][:], dist["theta"][itheta], dist["phi"][iphi], time=0.0)
                phi = np.ones_like(r) * dist["phi"][iphi]

                for irho in range(dist["rho"].size):
                    bphi = ascotpy.evaluate( r[irho], phi[irho], z[irho], 0.0, "bphi")
                    bnorm = ascotpy.evaluate(r[irho], phi[irho], z[irho], 0.0, "bnorm")
                    for ippar in range(dist["ppar"].size):
                        distj["distribution"][irho, itheta, iphi, ippar] = dist["distribution"][irho, itheta, iphi, ippar] * dist["ppar"][ippar]/ma * qa * bphi / bnorm
        dist = distmod.squeeze(distj, ppar=0)

    
    elif quantity == "ipowerdeposition" or quantity == "epowerdeposition" or quantity == "powerdeposition":
        
                
        # Convert to pitch-energy distribution.
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        # Integrate over pitch, time and charge dimensions
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)

        # Convert the energy abscissa into velocity
        va = np.sqrt( unyt.e*2*dist["energy"]/ma )

        # Evaluate the ASCOT collision operator at each spatial location for all velocities in va
        distP = copy.deepcopy(dist)
        for iphi in range(dist["phi"].size):
            for itheta in range(dist["theta"].size):

                #Calculate R,phi,z for all rho-values for single phi,theta
                (r, z) = ascotpy.get_rhotheta_rz(dist["rho"][:], dist["theta"][itheta], dist["phi"][iphi], time=0.0)
                phi = np.ones_like(r) * dist["phi"][iphi]

                for irho in range(dist["rho"].size):
                    
                    coefs = ascotpy.eval_collcoefs(ma, qa, r[irho],
                                                   phi[irho], z[irho], 0, va)
                    if   quantity == "ipowerdeposition":
                        distP["distribution"][irho, itheta, iphi, :] = -dist["distribution"][irho, itheta, iphi, :]*sum(coefs["K"][0,1:,:],0) *ma*va
                    elif quantity == "epowerdeposition": 
                        distP["distribution"][irho, itheta, iphi, :] = -dist["distribution"][irho, itheta, iphi, :]*    coefs["K"][0, 0,:]    *ma*va
                    elif quantity == "powerdeposition":
                        distP["distribution"][irho, itheta, iphi, :] = -dist["distribution"][irho, itheta, iphi, :]*sum(coefs["K"][0, :,:],0) *ma*va

        # Integrate over the energy dimension
        dist = distmod.squeeze(distP, energy=0)

    elif quantity == "colltorque":

        # Convert to pitch-energy distribution.
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        # Integrate over pitch, time and charge dimensions
        dist = distmod.squeeze(dist, time=0, charge=0)

        # Convert the energy abscissa into velocity
        va = np.sqrt( unyt.e*2*dist["energy"]/ma.v )

        # Evaluate the ASCOT collision operator at each spatial location for all velocities in va
        distT = copy.deepcopy(dist)
        for iphi in range(dist["phi"].size):
            for itheta in range(dist["theta"].size):

                #Calculate R,phi,z for all rho-values for single phi,theta
                (r, z) = ascotpy.get_rhotheta_rz(dist["rho"][:], dist["theta"][itheta], dist["phi"][iphi], time=0.0)
                phi = np.ones_like(r) * dist["phi"][iphi]

                for irho in range(dist["rho"].size):

                    coefs = ascotpy.eval_collcoefs(ma, qa, r[irho],
                                                   phi[irho], z[irho], 0, va)
                    Bnorm = ascotpy.evaluate(r[irho],phi[irho], z[irho], 0, "bnorm")
                    Bphi = ascotpy.evaluate(r[irho], phi[irho], z[irho], 0, "bphi")

                    v = np.tile(va,(49,1)).T
                    pitch = np.tile(dist["pitch"],(99,1))
                    K = np.tile(sum(coefs["K"][0,:,:],0),(49,1)).T
                    nu = np.tile(sum(coefs["nu"][0,:,:],0),(49,1)).T

                    deflectFreq = nu
                    dpitch = -deflectFreq*pitch
                    dp1 = K * ma
                    dPpara = pitch * dp1 + ma*v*dpitch
                    colltorque = dPpara * r[irho] * Bphi/Bnorm

                    distT["distribution"][irho, itheta, iphi, :] = dist["distribution"][irho, itheta, iphi, :]*-colltorque


        # Integrate over the energy dimension
        dist = distmod.squeeze(distT, energy=0, pitch=0)

    else:
        raise ValueError('Unknown quantity "{}".'.format(quantity))

    return dist


def eval_quantity_5d(ascotpy, dist, quantity, ma, qa):
    # Evaluate the requested quantity in 5D
    if quantity == "density":
        dist = distmod.squeeze(dist, ppar=0, pperp=0, time=0, charge=0)

    elif quantity == "energydensity":
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)

        distE = copy.deepcopy(dist)
        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):
                    distE["distribution"][ir, ip, iz, :] = dist["distribution"][ir,ip,iz,:] * dist["energy"] * unyt.e

        dist = distmod.squeeze(distE, energy=0)

    elif quantity == "pressure":
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)

        distp = copy.deepcopy(dist)
        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):
                    velocity2 = 2 * dist["energy"] * unyt.e / ma
                    distp["distribution"][ir, ip, iz, :] = dist["distribution"][ir,ip,iz,:] * ma * velocity2 / 3.0

        dist = distmod.squeeze(distp, energy=0)

    elif quantity == "poloidalcurrent":
        dist = distmod.squeeze(dist, pperp=0, time=0, charge=0)

        # current distribution
        distj = copy.deepcopy(dist)

        # integration over toroid angles
        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):

                    bphi = ascotpy.evaluate(dist["r"][ir], dist["phi"][ip], dist["z"][iz], 0.0, "bphi")
                    bnorm = ascotpy.evaluate(dist["r"][ir], dist["phi"][ip], dist["z"][iz], 0.0, "bnorm")
                    
                    for ippar in range(dist["ppar"].size):
                        distj["distribution"][ir, ip, iz, ippar] = dist["distribution"][ir,ip,iz,ippar] * dist["ppar"][ippar]/ma * qa * bphi / bnorm
                        
        dist = distmod.squeeze(distj, ppar=0)   
    
    elif quantity == "toroidalcurrent":
        dist = distmod.squeeze(dist, pperp=0, time=0, charge=0)

        distj = copy.deepcopy(dist)

        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):
                    bphi = ascotpy.evaluate(dist["r"][ir], dist["phi"][ip], dist["z"][iz], 0.0, "bphi")
                    bnorm = ascotpy.evaluate(dist["r"][ir], dist["phi"][ip], dist["z"][iz], 0.0, "bnorm")
                    for ippar in range(dist["ppar"].size):
                        distj["distribution"][ir, ip, iz, ippar] = dist["distribution"][ir,ip,iz,ippar] * dist["ppar"][ippar]/ma * qa * bphi / bnorm
        dist = distmod.squeeze(distj, ppar=0)

    elif quantity == "ipowerdeposition" or quantity == "epowerdeposition" or quantity == "powerdeposition":
        
        # Convert to pitch-energy distribution.
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        # Integrate over pitch, time and charge dimensions
        dist = distmod.squeeze(dist, pitch=0, time=0, charge=0)

        # Convert the energy abscissa into velocity
        va = np.sqrt( unyt.e*2*dist["energy"]/ma )

        # Evaluate the ASCOT collision operator at each spatial location for all velocities in va
        distP = copy.deepcopy(dist)
        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):
                    coefs = ascotpy.eval_collcoefs(ma, qa, dist["r"][ir],
                                                   dist["phi"][ip], dist["z"][iz], 0, va)
                    
                    if   quantity == "ipowerdeposition":
                        distP["distribution"][ir, ip, iz, :] = -dist["distribution"][ir, ip, iz, :]*sum(coefs["K"][0,1:,:],0) *ma*va
                    elif quantity == "epowerdeposition": 
                        distP["distribution"][ir, ip, iz, :] = -dist["distribution"][ir, ip, iz, :]*    coefs["K"][0, 0,:]    *ma*va
                    elif quantity == "powerdeposition":
                        distP["distribution"][ir, ip, iz, :] = -dist["distribution"][ir, ip, iz, :]*sum(coefs["K"][0, :,:],0) *ma*va
                    

        # Integrate over the energy dimension
        dist = distmod.squeeze(distP, energy=0)

    elif quantity == "colltorque":

        # Convert to pitch-energy distribution.
        dist = distconv.convert_ppappe_to_Exi(dist, ma, 100, 50)
        # Integrate over pitch, time and charge dimensions
        dist = distmod.squeeze(dist, time=0, charge=0)

        # Convert the energy abscissa into velocity
        va = np.sqrt( unyt.e*2*dist["energy"]/ma )

        # Evaluate the ASCOT collision operator at each spatial location for all velocities in va
        distT = copy.deepcopy(dist)
        for ir in range(dist["r"].size):
            for ip in range(dist["phi"].size):
                for iz in range(dist["z"].size):
                    coefs = ascotpy.eval_collcoefs(ma, qa, dist["r"][ir],
                                                   dist["phi"][ip], dist["z"][iz], 0, va)
                    Bnorm = ascotpy.evaluate(dist["r"][ir],dist["phi"][ip], dist["z"][iz], 0, "bnorm")
                    Bphi = ascotpy.evaluate(dist["r"][ir],dist["phi"][ip], dist["z"][iz], 0, "bphi")

                    v = np.tile(va,(49,1)).T
                    pitch = np.tile(dist["pitch"],(99,1))
                    K = np.tile(sum(coefs["K"][0,:,:],0),(49,1)).T
                    nu = np.tile(sum(coefs["nu"][0,:,:],0),(49,1)).T

                    deflectFreq = nu
                    dpitch = -deflectFreq*pitch
                    dp1 = K * ma
                    dPpara = pitch * dp1 + ma*v*dpitch
                    colltorque = dPpara * dist["r"][ir] * Bphi/Bnorm

                    distT["distribution"][ir, ip, iz, :] = dist["distribution"][ir, ip, iz, :]*colltorque

        # Integrate over the energy dimension
        dist = distmod.squeeze(distT, energy=0, pitch=0)


    else:
        raise ValueError('Unknown quantity "{}".'.format(quantity))

    return dist

def evalrhovol(ascotpy, rhomin, rhomax, nrho, rmin, rmax, nr,
               phimin, phimax, nphi, zmin, zmax, nz):

    # Volume of the test space (rectangular toroid)
    V_tot = 0.5 * (rmax**2 - rmin**2) * np.deg2rad(phimax - phimin) * (zmax - zmin)
    A_tot = (rmax - rmin) * (zmax - zmin)
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

    return (points_in_bins/n_tot) * V_tot, (points_in_bins/n_tot) * A_tot
