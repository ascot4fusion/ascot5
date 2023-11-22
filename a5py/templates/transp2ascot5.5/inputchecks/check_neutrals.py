# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:31:12 2023

@author: hvkg84

This script sanity checks the neutrals in the ASCOT5 input file.
Checks:
    1) Are the neutrals the same in the ASCOT5 file as in the TRANSP plasma state file and the full CDF file?

WIP: COME BACK AND IMPLEMENT THE FULL CDF METHOD
    
The comparisons for checking are:
    1) 1D density(rho_pol), A5 vs PS vs CDF
    2) 1D temperature(rho_pol), A5 vs PS vs CDF

NOTE: Maybe in the future, there should be a simple numerical check for each variable.
"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.interpolate
import inputhelpers.get_neutrals as get_neutrals
from a5py import Ascot


def compare_neutrals(ASCOT5_dataset, TRANSP_plasma_state_dataset, TRANSP_full_CDF_dataset, maxwellian=[1,1], time=None, extrapolate=True, extrapolation_type='dummy', dummy_value=5, plotyn=0, **kwargs):
    """
    Start by getting the correct time value for the CDF data. It can also be found from the plasma state file.
    
    neutrals_dict = {'rhomin': rhomin, 'rhomax': rhomax, 'nrho': nrho, 'nspecies': nspecies, 'anum': anum, 'znum': znum, 'density': neutral_density, \
                     'temperature': neutral_temperature, 'maxwellian': maxwellian, 'desc': desc}
    """
    if time == None:
        time = TRANSP_plasma_state_dataset.variables['t1'][:]
    CDF_time_idx = np.max(np.where( TRANSP_full_CDF_dataset.variables['TIME'][:] <= time ))
    
    ### FIRST get the neutral species values from the ASCOT5 file
    A5_temperature = ASCOT5_dataset.data.neutral.active.read_hdf5()['temperature'][:] ## Note that this has shape (nrho, nspecies)
    A5_density = ASCOT5_dataset.data.neutral.active.read_hdf5()['density'][:] ## Note that this has shape (nrho, nspecies)
    A5_rhomin = ASCOT5_dataset.data.neutral.active.read_hdf5()['rhomin'][:].squeeze()
    A5_rhomax = ASCOT5_dataset.data.neutral.active.read_hdf5()['rhomax'][:].squeeze()
    A5_nrho = ASCOT5_dataset.data.neutral.active.read_hdf5()['nrho'][:].squeeze()
    A5_anum = ASCOT5_dataset.data.neutral.active.read_hdf5()['anum'][:]
    A5_znum = ASCOT5_dataset.data.neutral.active.read_hdf5()['znum'][:]
    A5_nspecies = ASCOT5_dataset.data.neutral.active.read_hdf5()['nspecies'][:]
    A5_maxwellian = ASCOT5_dataset.data.neutral.active.read_hdf5()['maxwellian'][:]
    A5_rho_pol = np.linspace(A5_rhomin, A5_rhomax, A5_nrho)
    
    ### THEN get the neutral species values from the plasma state file
    PS_data = get_neutrals.neutrals_from_plasma_state_dataset(TRANSP_plasma_state_dataset, maxwellian=maxwellian, extrapolate=extrapolate, \
                                                          extrapolation_type=extrapolation_type, dummy_value=dummy_value, plotyn=0, **kwargs)
    PS_temperature = PS_data['temperature'][:] ## Note that this has shape (nrho, nspecies)
    PS_density = PS_data['density'][:] ## Note that this has shape (nrho, nspecies)
    PS_rhomin = PS_data['rhomin']
    PS_rhomax = PS_data['rhomax']
    PS_nrho = PS_data['nrho']
    PS_anum = PS_data['anum'][:]
    PS_znum = PS_data['znum'][:]
    PS_nspecies = PS_data['nspecies']
    PS_maxwellian = PS_data['maxwellian'][:]
    PS_rho_pol = np.linspace(PS_rhomin, PS_rhomax, PS_nrho)
    
    ### THEN get the neutral species values from the plasma state file
    # CDF_data = get_neutrals.neutrals_from_full_CDF(TRANSP_full_CDF_dataset, time=time, maxwellian=maxwellian, extrapolate=extrapolate, \
    #                                            extrapolation_type=extrapolation_type, dummy_value=dummy_value, plotyn=0, **kwargs)
    # CDF_temperature = CDF_data['temperature'][:] ## Note that this has shape (nrho, nspecies)
    # CDF_density = CDF_data['density'][:] ## Note that this has shape (nrho, nspecies)
    # CDF_rhomin = CDF_data['rhomin']
    # CDF_rhomax = CDF_data['rhomax']
    # CDF_nrho = CDF_data['nrho']
    # CDF_anum = CDF_data['anum'][:]
    # CDF_znum = CDF_data['znum'][:]
    # CDF_nspecies = CDF_data['nspecies']
    # CDF_maxwellian = CDF_data['maxwellian'][:]
    # CDF_rho_pol = np.linspace(CDF_rhomin, CDF_rhomax, CDF_nrho)
    
    d = 4 ## decimal places for rounding
    tol = 1e-2
    
    other_vars = ['anum', 'znum', 'nspecies', 'maxwellian']
    
    # other_vars_equal = all([(all(np.unique(np.atleast_1d(ASCOT5_dataset.plasma.active.read_hdf5()[key]).squeeze() - np.atleast_1d(PS_data[key]).squeeze() < tol)), \
    #                          all(np.unique( np.atleast_1d(PS_data[key]).squeeze() -  np.atleast_1d(CDF_data[key]).squeeze() < tol))) for key in other_vars])
    
    other_vars_equal = all([(all(np.unique(np.atleast_1d(ASCOT5_dataset.data.neutral.active.read_hdf5()[key]).squeeze() - np.atleast_1d(PS_data[key]).squeeze() < tol))) for key in other_vars])
    print("Other values match within "+str(tol)+":", other_vars_equal)
    
    
    if plotyn == 1:
        plt.figure()
        plt.plot(A5_rho_pol, A5_temperature[:,0], label='ASCOT5, species 1', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_temperature[:,0], label='Plasma state, species 1', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        # plt.plot(CDF_rho_pol, CDF_temperature[:,0], label='CDF, species 1', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral temperature [eV]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_temperature[:,1], label='ASCOT5, species 2', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_temperature[:,1], label='Plasma state, species 2', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        # plt.plot(CDF_rho_pol, CDF_temperature[:,1], label='CDF, species 2', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral temperature [eV]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_density[:,0], label='ASCOT5, species 1', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_density[:,0], label='Plasma state, species 1', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        # plt.plot(CDF_rho_pol, CDF_density[:,0], label='CDF, species 1', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral density [m$^{-3}$]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_density[:,1], label='ASCOT5, species 2', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_density[:,1], label='Plasma state, species 2', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        # plt.plot(CDF_rho_pol, CDF_density[:,1], label='CDF, species 2', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Neutral density [m$^{-3}$]')
        plt.legend()
        
        
        plt.show()
        
        
        
            
    return


if __name__ == '__main__':
    plotyn = 1
    
    run_ID = '134020D30'
    TRANSP_plasma_state_directory = '../../134020D30/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    TRANSP_full_CDF_filename = TRANSP_plasma_state_directory + run_ID + ".CDF"
    TRANSP_full_CDF_data = nc.Dataset(TRANSP_full_CDF_filename)
    
    ASCOT5_filename = "../runs/NSTX_134020_GC_new_tests_8.h5"
    ASCOT5_data = Ascot(ASCOT5_filename)

    comparison = compare_neutrals(ASCOT5_dataset=ASCOT5_data, TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, \
                                TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, time=None, extrapolate=True, \
                                extrapolation_type='constant', dummy_value=1e0, plotyn=plotyn)
    
    if plotyn == 1:
        
        
        plt.show()


