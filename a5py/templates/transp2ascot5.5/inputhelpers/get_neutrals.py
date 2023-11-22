# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 17:12:42 2023

@author: Ian Dolby

WIP: COME BACK AND ADD RECOMBINATION AND HALO NEUTRALS... CLEAR UP HOW EXACTLY TO LOAD ALL OF THE NEUTRALS -- i.e sort out after Finland
    ALSO: PROPERLY IMPLEMENT THE FULL CDF FUNCTION, AND ADD A DOCSTRING IN EACH FUNCTION

ASCOT5 neutrals 1D: 
    
    fn : str <br>
            Full path to the HDF5 file.
        rhomin : float <br>
            Minimum value in rho grid [1].
        rhomax : float <br>
            Maximum value in rho grid [1].
        nrho : int <br>
            Number of rho grid points.
        nspecies : int <br>
            Number of neutral species.
        anum : array_like (nspecies,1) <br>
            Neutral species' atomic mass number.
        znum array_like (nspecies,1) <br>
            Neutral species' charge number.
        density array_like (nrho,nspecies) <br>
            Neutral species-wise density [m^-3].
        temperature array_like (nrho,nspecies) <br>
            Neutral species-wise temperature [eV].
        maxwellian array_like (nspecies,1) <br> :
            Whether species distribution is Maxwellian (1) of monoenergetic (0)
        desc : str, optional <br>
            Input description.

"""

import numpy as np
import netCDF4 as nc
import scipy.interpolate
import matplotlib.pyplot as plt
from .extrapolate_profile import extrapolate_profile

from scipy.constants import physical_constants as const



def neutrals_from_plasma_state_dataset(TRANSP_plasma_state_dataset, maxwellian=[1,1], extrapolate=True, extrapolation_type='dummy', rho_pol_R_map=None, dummy_value=5, plotyn=0, **kwargs):
    
    TPSD = TRANSP_plasma_state_dataset
    
    """ First, load the rho_toroidal values and use psipol to convert them to rho_pol. We load PsiRZ and R_axis and Z_axis to find psi0.
    """
    rho_tor_edges = TPSD.variables['rho_gas'][:]
    psipol = TPSD.variables['psipol'][:]  # Get the poloidal flux values which correspond to the rho_toroidal values (the profiles are provided in terms of rho_toroidal)
    psi_interp = scipy.interpolate.RegularGridInterpolator(points=(TPSD.variables['R_grid'][:], TPSD.variables['Z_grid'][:]), values=TPSD.variables['PsiRZ'][:].T, method='cubic')
    psi0  = psi_interp((TPSD.variables['R_axis'][:], TPSD.variables['Z_axis'][:]))
    psi1  = psipol[-1]
    rho_tor_centres = rho_tor_edges[:-1] + np.diff(rho_tor_edges) / 2 # Shift all of the rho-values by half a bin
    # psipol is already defined
    psipol_centres_vs_rho_tor_centres = np.interp(x = rho_tor_centres, xp = rho_tor_edges, fp = psipol)
    rho_pol_centres_vs_rho_tor_centres = np.sqrt((psipol_centres_vs_rho_tor_centres - psi0)/(psi1 - psi0)) # Remember that psi, psi0, and psi1 have been converted to COCOS3, unlike psipol.
    
    
    """
    Each of these next variables has shape [no. of sources, no. of species, no. of rho values].
    
    The temperature is basically the same from source to source.
    We sum over the sources to get the total densities, while averaging over the (mostly identical)
     temperature values for each source to get the temperatures.
    """
    neutral_densities_per_flux = TPSD.variables['n0norm'][:]
    neutral_influxes = TPSD.variables['sc0'][:]
    neutral_temperatures = TPSD.variables['T0sc0'][:] * 1e3 # CONVERT to eV FROM keV
    n_gas = len(TPSD.variables['SGAS_type'][:])
    n_sc0 = len(neutral_influxes)

    # sc0 times n0norm
    neutral_density = np.array( [ np.sum( np.array( [neutral_densities_per_flux[s,g,:] *  neutral_influxes[s] for s in range(n_sc0)] ), axis=0) for g in range(n_gas) ] ) # s for 'source', g for 'gas'
    neutral_temperature = np.array( [ np.mean(neutral_temperatures, axis=0)[g] for g in range(n_gas)] )
    ## ^ These should now have shape (n_gas, nrho).
    
    """
    Now append the boundary values to the bulk profiles, assuming that the proportion of D to H is the same (dn0out is the TOTAL neutral density beyond the LCFS).
    Assume that the temperature is unchanged.
    """
    # H_D_ratio just before LCFS
    neutral_density_norm = np.sum([neutral_density[g,-1] for g in range(n_gas)], axis=0)
    r = neutral_density[:,-1]/neutral_density_norm
    N_tot = TPSD.variables['dn0out'][:]
    
    N_LCFS = r * N_tot
    
    neutral_density = np.append(neutral_density.T, np.atleast_2d(N_LCFS), axis=0)
    neutral_temperature = np.append(neutral_temperature.T, np.atleast_2d(neutral_temperature[:,-1]), axis=0)
    
    """
    Now prepend 0 and append 1 to rho_poloidal, and adjust the profiles accordingly.
    """
    rho_pol_centres_vs_rho_tor_centres = np.append(0, rho_pol_centres_vs_rho_tor_centres) # Now add back the central bin to ensure we still have a value defined at rho=0
    rho_pol_centres_vs_rho_tor_centres = np.append(rho_pol_centres_vs_rho_tor_centres, 1) # Now add the LCFS value so that the boundary values for the profiles can be added.
    
    neutral_density = np.append(np.atleast_2d(neutral_density[0]), neutral_density, axis=0)
    neutral_temperature = np.append(np.atleast_2d(neutral_temperature[0]), neutral_temperature, axis=0)
    
    if extrapolate == True:
        """ Dummy extrapolation means that we append small non-zero values from just beyond the LCFS out to the edge of the grid.
            Constant means that the last value will be copied out to the edge of the grid.
            Linear means that the slope of the last elements will be used, with a warning if the the slope is positive.
            Exponential means that a SOL model will be used to calculate a decay constant for an exponential decay as a function of R_midplane.
        """
        rho_max = kwargs.pop("rho_to_append", 2)
        rho_nsteps = kwargs.pop("rho_nsteps", 10)
        
        extrapolation_args = {'rho_array':rho_pol_centres_vs_rho_tor_centres, 'rho_to_append':rho_max, 'rho_nsteps':rho_nsteps, 'dummy_value':dummy_value, \
                              'extrapolation_type':extrapolation_type, 'plotyn':plotyn, **kwargs}
            
        if extrapolation_type == 'exponential':
            ## set up all the required variables for an exponential decay.
            theta_zero_index = np.min(np.where(TPSD.variables['th_eq'][:] >= 0 ) )
            R_LCFS = TPSD.variables['R_geo'][theta_zero_index, -1]
            Z_axis = TPSD.variables['Z_axis'][:]
            Bpol = np.sqrt( TPSD.variables['BRRZ'][:]**2 + TPSD.variables['BZRZ'][:]**2).T
            Bpol_interp = scipy.interpolate.RegularGridInterpolator(points=(TPSD.variables['R_grid'][:], TPSD.variables['Z_grid'][:]), values=Bpol, method='cubic')
            extrapolation_args['Bpol_LCFS'] = Bpol_interp((R_LCFS, Z_axis))
            extrapolation_args['T_LCFS'] = TPSD.variables['Ti_bdy'][:]
            extrapolation_args['R_LCFS'] = R_LCFS
            extrapolation_args['R_axis'] = TPSD.variables['R_axis'][:]
            
        
        ## NOW actually extrapolate the variables.
        neutral_density_extrap_array = np.zeros((len(neutral_density) + rho_nsteps, n_gas))
        neutral_temperature_extrap_array = np.zeros((len(neutral_temperature) + rho_nsteps, n_gas))
        
        for column in range(n_gas):
            density_profile_to_add = extrapolate_profile(variable_array = neutral_density[:,column].flatten(), **extrapolation_args)['quantity']
            temperature_profile_to_add = extrapolate_profile(variable_array = neutral_temperature[:,column].flatten(), **extrapolation_args)['quantity']
            
            neutral_density_extrap_array[:, column] = density_profile_to_add.reshape(np.shape(neutral_density_extrap_array[:,0]))
            neutral_temperature_extrap_array[:, column] = temperature_profile_to_add.reshape(np.shape(neutral_temperature_extrap_array[:,0]))
            
            rho_pol_centres_vs_rho_tor_centres_extrap = extrapolate_profile(variable_array = neutral_density[:,column].flatten(), **extrapolation_args)['rho'] ## This really only needs to be called once.
        
        
        rho_pol_centres_vs_rho_tor_centres = rho_pol_centres_vs_rho_tor_centres_extrap
        neutral_density = neutral_density_extrap_array
        neutral_temperature = neutral_temperature_extrap_array

    
    """ Now prepare the remaining variables.
    """
    rhomin = np.min(rho_pol_centres_vs_rho_tor_centres)
    rhomax = np.max(rho_pol_centres_vs_rho_tor_centres)
    nrho = len(rho_pol_centres_vs_rho_tor_centres)
    nspecies = n_gas
    mass = TPSD.variables['m_SGAS'][:] / const['atomic mass constant'][0]
    anum = np.round(mass, 0)
    znum = TPSD.variables['q_SGAS'][:] / const['elementary charge'][0]
    maxwellian = np.array(maxwellian)
    
    ## Now we interpolate onto a fine rho_pol grid; the current one is irregular.
    rho_pol_new = np.linspace(rhomin, rhomax, nrho)
    temperature_new = np.array([ np.interp(x=rho_pol_new, xp=rho_pol_centres_vs_rho_tor_centres, fp=neutral_temperature[:,g]) for g in range(nspecies) ]).T
    density_new = np.array([ np.interp(x=rho_pol_new, xp=rho_pol_centres_vs_rho_tor_centres, fp=neutral_density[:,g]) for g in range(nspecies) ]).T
    
    assert np.shape(density_new) == np.shape(neutral_density)
    
    
    ## Now get the species names to put in the description
    species_names = ', '.join(bytes(np.array(TPSD.variables['SGAS_name'][:,0]).flatten()).decode('UTF-8'))
    desc = 'Neutrals from state file, species ' + species_names
    
    if plotyn == 1:
        plt.figure()
        plt.plot(rho_pol_new, density_new[:,0], label='ASCOT5, new rho_pol, species 1', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(rho_pol_centres_vs_rho_tor_centres, neutral_density[:,0], label='ASCOT5, old rho_pol, species 1', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral density [m$^{-3}$]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(rho_pol_new, density_new[:,1], label='ASCOT5, new rho_pol, species 2', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(rho_pol_centres_vs_rho_tor_centres, neutral_density[:,1], label='ASCOT5, old rho_pol, species 2', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral density [m$^{-3}$]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(rho_pol_new, temperature_new[:,0], label='ASCOT5, new rho_pol, species 1', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(rho_pol_centres_vs_rho_tor_centres, neutral_temperature[:,0], label='ASCOT5, old rho_pol, species 1', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral temperature [eV]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(rho_pol_new, temperature_new[:,1], label='ASCOT5, new rho_pol, species 2', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(rho_pol_centres_vs_rho_tor_centres, neutral_temperature[:,1], label='ASCOT5, old rho_pol, species 2', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral temperature [eV]')
        plt.legend()

    
    neutrals_dict = {'rhomin': rhomin, 'rhomax': rhomax, 'nrho': nrho, 'nspecies': nspecies, 'anum': anum, 'znum': znum, 'density': density_new, \
                     'temperature': temperature_new, 'maxwellian': maxwellian, 'desc': desc}
    return neutrals_dict



if __name__ == '__main__':
    plotyn = 1
    
    run_ID = '134020D30'
    TRANSP_plasma_state_directory = '../../134020D30/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    
    neutrals_test_dict_PS = neutrals_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, extrapolate=True, extrapolation_type='constant', rho_pol_R_map=None, dummy_value=1e-0, plotyn=plotyn)
    
    
    if plotyn == 1:
        
        
        
        plt.show()




