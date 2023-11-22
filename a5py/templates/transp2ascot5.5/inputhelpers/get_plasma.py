# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 18:52:15 2023

@author: Ian Dolby


WIP: ADD A DOCSTRING IN EACH FUNCTION

"""

import numpy as np
import netCDF4 as nc
import scipy.interpolate
import matplotlib.pyplot as plt
from .extrapolate_profile import extrapolate_profile

from scipy.constants import physical_constants as const

def plasma_from_plasma_state_dataset(TRANSP_plasma_state_dataset, extrapolate=True, extrapolation_type='dummy', rho_pol_R_map=None, dummy_value=5, plotyn=0, **kwargs):
    
    TPSD = TRANSP_plasma_state_dataset
    
    """ First, load the rho_toroidal values and use psipol to convert them to rho_pol. We load PsiRZ and R_axis and Z_axis to find psi0.
    """
    rho_tor_edges = TPSD.variables['rho'][:]
    psipol = TPSD.variables['psipol'][:]  # Get the poloidal flux values which correspond to the rho_toroidal values (the profiles are provided in terms of rho_toroidal)
    psi_interp = scipy.interpolate.RegularGridInterpolator(points=(TPSD.variables['R_grid'][:], TPSD.variables['Z_grid'][:]), values=TPSD.variables['PsiRZ'][:].T, method='cubic')
    psi0  = psi_interp((TPSD.variables['R_axis'][:], TPSD.variables['Z_axis'][:]))
    psi1  = psipol[-1]
    rho_tor_centres = rho_tor_edges[:-1] + np.diff(rho_tor_edges) / 2 # Shift all of the rho-values by half a bin
    # psipol is already defined
    psipol_centres_vs_rho_tor_centres = np.interp(x = rho_tor_centres, xp = rho_tor_edges, fp = psipol)
    rho_pol_centres_vs_rho_tor_centres = np.sqrt((psipol_centres_vs_rho_tor_centres - psi0)/(psi1 - psi0))
    
    
    """ Secondly, load the bulk plasma values. We assume that all of the bulk plasma is fully ionised, so that znum = |charge|.
         We start by obtaining the index corresponding to the electrons by finding where the charge is equal to -e.
        
        Note that ASCOT5 requires the ion density profiles to have shape (nrho, nspecies), so these variables are transposed. It also makes appending values simpler.
        
        Also note that the plasma rotation is NOT included in this script.
    """
    e = const['elementary charge'][0]
    electron_index = np.where(np.round(TPSD.variables['q_S'][:]/(-e), 3) == 1)[0]
    plasma_indices = np.where(np.round(TPSD.variables['q_S'][:]/(-e), 3) != 1)[0]
    
    ## Now we modify the plasma indices to eliminate any species whose density is 0 to remove species which were not actually used in the simulation.
    plasma_indices = np.array( [idx for idx in plasma_indices if np.round(np.mean(TPSD.variables['ns'][idx]), 0) != 0] )
    
    TE = TPSD.variables['Ts'][electron_index] * 1e3  #CONVERT TO eV from keV
    TI = TPSD.variables['Ti'][:] * 1e3  #CONVERT TO eV from keV
    
    NE = TPSD.variables['ns'][electron_index]   # NO CONVERSION NECESSARY; THESE ARE ALREADY IN m^-3
    NI = TPSD.variables['ns'][plasma_indices].T
    
    # Values at LCFS (or "bdy" for "boundary"):
    TE_bdy = TPSD.variables['Te_bdy'][:] * 1e3  #CONVERT TO eV from keV
    TI_bdy = TPSD.variables['Ti_bdy'][:] * 1e3  #CONVERT TO eV from keV
    
    NE_bdy = TPSD.variables['ns_bdy'][electron_index]
    NI_bdy = TPSD.variables['ns_bdy'][plasma_indices]
    
    
    """ Now append the boundary values to the bulk profiles
    """
    TE = np.append(TE, TE_bdy)
    TI = np.append(TI, TI_bdy)
    
    NE = np.append(NE, NE_bdy)
    NI = np.append(NI, np.atleast_2d(NI_bdy), axis=0)
    
    """ Now prepend 0 and append 1 to rho_poloidal, and adjust the profiles accordingly.
    """
    rho_pol_centres_vs_rho_tor_centres = np.append(0, rho_pol_centres_vs_rho_tor_centres) # Now add back the central bin to ensure we still have a value defined at rho=0
    rho_pol_centres_vs_rho_tor_centres = np.append(rho_pol_centres_vs_rho_tor_centres, 1) # Now add the LCFS value so that the boundary values for the profiles can be added.
    
    TE = np.append(TE[0], TE)
    TI = np.append(TI[0], TI)
    
    NE = np.append(NE[0], NE)
    NI = np.append(np.atleast_2d(NI[0]), NI, axis=0)
    
    
    """ Now get the atomic number, atomic mass, charge, and mass number of each species, assuming znum = |charge|.
    """
    mass = TPSD.variables['m_S'][plasma_indices] / const['atomic mass constant'][0]
    charge = TPSD.variables['q_S'][plasma_indices] / e
    anum = np.round(mass, 0)
    znum = np.abs(charge)
    Nion = len(znum)
    
    
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
            extrapolation_args['T_LCFS'] = TI_bdy
            extrapolation_args['R_LCFS'] = R_LCFS
            extrapolation_args['R_axis'] = TPSD.variables['R_axis'][:]
        
        ## NOW actually extrapolate the variables.
        TE_extrap = extrapolate_profile(variable_array = TE, **extrapolation_args)['quantity']
        TI_extrap = extrapolate_profile(variable_array = TI, **extrapolation_args)['quantity']
        NE_extrap = extrapolate_profile(variable_array = NE, **extrapolation_args)['quantity']
        rho_pol_centres_vs_rho_tor_centres_extrap = extrapolate_profile(variable_array = TE, **extrapolation_args)['rho']
        NI_extrap_array = np.zeros((len(TI) + rho_nsteps, Nion))
        
        for column in range(Nion):
            profile_to_add = extrapolate_profile(variable_array = NI[:,column].flatten(), **extrapolation_args)['quantity']
            NI_extrap_array[:, column] = profile_to_add.reshape(np.shape(NI_extrap_array[:,0]))
        
        rho_pol_centres_vs_rho_tor_centres = rho_pol_centres_vs_rho_tor_centres_extrap
        TE = TE_extrap
        TI = TI_extrap
        NE = NE_extrap
        NI = NI_extrap_array
    
    
    
    
    ## Now calculate the length of rho
    Nrho = len(rho_pol_centres_vs_rho_tor_centres)
    
    if plotyn == 1:
        plt.figure()
        plt.plot(rho_pol_centres_vs_rho_tor_centres, NE, label='ne vs rho_poloidal')
        plt.xlabel('rho')
        plt.ylabel('n_e [m^-3]')
        plt.legend()
    
    
    ## Now get the species names to put in the description
    species_indices_full = np.append(electron_index, plasma_indices)
    species_names = ', '.join(bytes(np.array(TPSD.variables['S_name'][species_indices_full][:,0]).flatten()).decode('UTF-8'))
    desc = 'Plasma from state file, species ' + species_names
    
    plasma_dict = {'nrho': Nrho, 'nion': Nion, 'anum': anum, 'znum': znum, 'mass': mass, 'charge': charge, 'rho': rho_pol_centres_vs_rho_tor_centres, \
                   'edensity': NE, 'etemperature': TE, 'idensity': NI, 'itemperature': TI, 'desc':desc}
    return plasma_dict



def plasma_from_full_CDF(TRANSP_full_CDF_dataset, time=0.65, extrapolate=True, extrapolation_type='dummy', rho_pol_R_map=None, dummy_value=5, plotyn=0, **kwargs):

    TFCD = TRANSP_full_CDF_dataset
    
    time_index = np.min(np.where(TFCD.variables['TIME'][:] > time))
    # Get variables from the netCDF file
    PLFLX = TFCD.variables['PLFLX'][time_index] # A function of XB
    X  = TFCD.variables['X'][time_index,:]  # Zone centres
    XB = TFCD.variables['XB'][time_index] # Zone edges
    XB = np.append(0, XB)  # Prepend the magnetic axis
    XIRSYM = TFCD.variables['XIRSYM'][time_index,:] # X-values corresponding to the R-values contained in RMJSYM -- NOTE that these values are negative on the high field side, but only one side is needed.
    RMJSYM = TFCD.variables['RMJSYM'][time_index,:] * 1e-2 #CONVERT TO M      # R-values for data
    
    
    RHO_pol_vs_rho_tor_edges = np.append(0, np.sqrt((PLFLX - 0)/(PLFLX[-1] - 0))) ## NOTE that we are assuming psi on-axis to be ZERO....
    RHO_pol_vs_rho_tor_centres = np.interp(x=X, xp=XB, fp=RHO_pol_vs_rho_tor_edges)
    
    
    TE = TFCD.variables['TE'][time_index]
    TI = TFCD.variables['TI'][time_index]
    NE = TFCD.variables['NE'][time_index] * 1e6 #CONVERT TO PER METRES CUBED
    Ni_H = TFCD.variables['NH'][time_index] * 1e6
    Ni_D = TFCD.variables['ND'][time_index] * 1e6
    Ni_imp = TFCD.variables['NIMP_SINGL'][time_index] * 1e6
    Te0 = TFCD.variables['TE0'][time_index]
    Ti0 = TFCD.variables['TI0'][time_index]
    
    ## Prepend axis values
    RHO_pol_vs_rho_tor_centres = np.append(0, RHO_pol_vs_rho_tor_centres)
    TE = np.append(Te0, TE)
    TI = np.append(Ti0, TI)
    NE = np.append(NE[0], NE)
    Ni_H = np.append(Ni_H[0], Ni_H)
    Ni_D = np.append(Ni_D[0], Ni_D)
    Ni_imp = np.append(Ni_imp[0], Ni_imp)
    
    
    Te_edge = TFCD.variables['TEEDG'][time_index]
    Ti_edge = TFCD.variables['TIEDG'][time_index]
    
    # COLUMNS: HYDROGEN, DEUTERIUM, IMPURITY
    NI = np.vstack( (Ni_H, Ni_D, Ni_imp) ).T ## Required shape is (nrho, nspecies)
    
    
    """ Now get the atomic number, atomic mass, charge, and mass number of each species, assuming znum = |charge|.
    """
    anum = np.array([ 1, 2, TFCD.variables['AIMP'][time_index]])
    znum = np.array([ 1, 1, TFCD.variables['XZIMP'][time_index]])
    charge = znum
    mass = (anum - znum) * const['neutron relative atomic mass'][0] + znum * const['proton relative atomic mass'][0]
    Nion = len(znum)
    
    
    if extrapolate == True:
        """ Dummy extrapolation means that we append small non-zero values from just beyond the LCFS out to the edge of the grid.
            Constant means that the last value will be copied out to the edge of the grid.
            Linear means that the slope of the last elements will be used, with a warning if the the slope is positive.
            Exponential means that a SOL model will be used to calculate a decay constant for an exponential decay as a function of R_midplane.
        """
        rho_max = kwargs.pop("rho_to_append", 2)
        rho_nsteps = kwargs.pop("rho_nsteps", 10)
        
        extrapolation_args = {'rho_array':RHO_pol_vs_rho_tor_centres, 'rho_to_append':rho_max, 'rho_nsteps':rho_nsteps, 'dummy_value':dummy_value, \
                              'extrapolation_type':extrapolation_type, 'plotyn':plotyn, **kwargs}
        
        if extrapolation_type == 'exponential':
            ## set up all the required variables for an exponential decay.
            LCFS_idx = np.max(np.where(XIRSYM <= 1))
            R_LCFS = RMJSYM[LCFS_idx]
            Z_axis = TFCD.variables['YAXIS'][time_index]
            Bpol = TFCD.variables['BPOL'][time_index][-1]
            extrapolation_args['Bpol_LCFS'] = Bpol
            extrapolation_args['T_LCFS'] = Ti_edge
            extrapolation_args['R_LCFS'] = R_LCFS*1e-2 # CONVERT TO M
            extrapolation_args['R_axis'] = TFCD.variables['RAXIS'][time_index]*1e-2 # CONVERT TO M
        
        ## NOW actually extrapolate the variables.
        TE_extrap = extrapolate_profile(variable_array = TE, **extrapolation_args)['quantity']
        TI_extrap = extrapolate_profile(variable_array = TI, **extrapolation_args)['quantity']
        NE_extrap = extrapolate_profile(variable_array = NE, **extrapolation_args)['quantity']
        rho_pol_centres_vs_rho_tor_centres_extrap = extrapolate_profile(variable_array = TE, **extrapolation_args)['rho']
        NI_extrap_array = np.zeros((len(TI) + rho_nsteps, Nion))
        
        for column in range(Nion):
            profile_to_add = extrapolate_profile(variable_array = NI[:,column].flatten(), **extrapolation_args)['quantity']
            NI_extrap_array[:, column] = profile_to_add.reshape(np.shape(NI_extrap_array[:,0]))
        
        RHO_pol_vs_rho_tor_centres = rho_pol_centres_vs_rho_tor_centres_extrap
        TE = TE_extrap
        TI = TI_extrap
        NE = NE_extrap
        NI = NI_extrap_array
    
    
    Nrho = len(RHO_pol_vs_rho_tor_centres)
    desc = 'Plasma from state file, species H, D, impurity ions with Z = ' + str(znum[2:])  # ( HYDROGEN, DEUTERIUM, IMPURITY)
    
    plasma_dict = {'nrho': Nrho, 'nion': Nion, 'anum': anum, 'znum': znum, 'mass': mass, 'charge': charge, 'rho': RHO_pol_vs_rho_tor_centres, \
                   'edensity': NE, 'etemperature': TE, 'idensity': NI, 'itemperature': TI, 'desc':desc}  
    return plasma_dict





if __name__ == '__main__':
    plotyn = 1
    
    run_ID = '134020D30'
    TRANSP_plasma_state_directory = '../../134020D30/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    TRANSP_full_CDF_filename = TRANSP_plasma_state_directory + run_ID + ".CDF"
    TRANSP_full_CDF_data = nc.Dataset(TRANSP_full_CDF_filename)
    
    plasma_test_dict_PS = plasma_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, extrapolate=True, extrapolation_type='exponential', rho_pol_R_map=None, dummy_value=1e-2, plotyn=0)
    plasma_test_dict_CDF = plasma_from_full_CDF(TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, extrapolate=True, extrapolation_type='exponential', rho_pol_R_map=None, dummy_value=1e-2, plotyn=0)
    
    
    
    if plotyn == 1:
        plt.figure()
        plt.plot(plasma_test_dict_PS['rho'], plasma_test_dict_PS['edensity'], label='NE from PS')
        plt.plot(plasma_test_dict_CDF['rho'], plasma_test_dict_CDF['edensity'], label='NE from CDF')
        plt.legend()
        
        plt.show()




