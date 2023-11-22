# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:31:12 2023

@author: Ian Dolby

This script sanity checks the plasma data in the ASCOT5 input file.
Checks:
    A) Is the plasma the same in the ASCOT5 file as in the TRANSP plasma state file and the full CDF file?

The comparisons for checking are:
    1) 1D Te(rho_pol), A5 vs PS vs CDF
    2) 1D Ti(rho_pol), A5 vs PS vs CDF
    3) 1D Ne(rho_pol), A5 vs PS vs CDF
    4) 1D Ni(rho_pol), A5 vs PS vs CDF

"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.interpolate
import inputhelpers.get_plasma as get_plasma
from a5py import Ascot


def compare_plasma(ASCOT5_dataset, TRANSP_plasma_state_dataset, TRANSP_full_CDF_dataset, time=None, extrapolate=True, extrapolation_type='dummy', dummy_value=5, plotyn=0, **kwargs):
    """
    Start by getting the correct time value for the CDF data. It can also be found from the plasma state file.
    """
    if time == None:
        time = TRANSP_plasma_state_dataset.variables['t1'][:]
    CDF_time_idx = np.max(np.where( TRANSP_full_CDF_dataset.variables['TIME'][:] <= time ))
    
    ### FIRST get the plasma values from the ASCOT5 file
    A5_Te = ASCOT5_dataset.data.plasma.active.read()['etemperature'][:]
    A5_Ti = ASCOT5_dataset.data.plasma.active.read()['itemperature'][:]
    A5_Ne = ASCOT5_dataset.data.plasma.active.read()['edensity'][:]
    A5_Ni = ASCOT5_dataset.data.plasma.active.read()['idensity'][:] ## Note that this has shape (nrho, nspecies)
    A5_rho_pol = ASCOT5_dataset.data.plasma.active.read()['rho'][:]
    A5_anum = ASCOT5_dataset.data.plasma.active.read()['anum'][:]
    A5_znum = ASCOT5_dataset.data.plasma.active.read()['znum'][:]
    A5_nion = ASCOT5_dataset.data.plasma.active.read()['nion']#[:]
    A5_mass = ASCOT5_dataset.data.plasma.active.read()['mass'][:]
    A5_charge = ASCOT5_dataset.data.plasma.active.read()['charge'][:]
    
    
    ### THEN get the plasma values from the plasma state file
    PS_data = get_plasma.plasma_from_plasma_state_dataset(TRANSP_plasma_state_dataset, extrapolate=extrapolate, \
                                                          extrapolation_type=extrapolation_type, dummy_value=dummy_value, plotyn=0, **kwargs)
    PS_Te = PS_data['etemperature'][:]
    PS_Ti = PS_data['itemperature'][:]
    PS_Ne = PS_data['edensity'][:]
    PS_Ni = PS_data['idensity'][:] ## Note that this has shape (nrho, nspecies)
    PS_rho_pol = PS_data['rho'][:]
    PS_anum = PS_data['anum'][:]
    PS_znum = PS_data['znum'][:]
    PS_nion = PS_data['nion']
    PS_mass = PS_data['mass'][:]
    PS_charge = PS_data['charge'][:]
    
    
    ### THEN get the plasma values from the plasma state file
    CDF_data = get_plasma.plasma_from_full_CDF(TRANSP_full_CDF_dataset, time=time, extrapolate=extrapolate, \
                                               extrapolation_type=extrapolation_type, dummy_value=dummy_value, plotyn=0, **kwargs)
    CDF_Te = CDF_data['etemperature'][:]
    CDF_Ti = CDF_data['itemperature'][:]
    CDF_Ne = CDF_data['edensity'][:]
    CDF_Ni = CDF_data['idensity'][:] ## Note that this has shape (nrho, nspecies)
    CDF_rho_pol = CDF_data['rho'][:]
    CDF_anum = CDF_data['anum'][:]
    CDF_znum = CDF_data['znum'][:]
    CDF_nion = CDF_data['nion']
    CDF_mass = CDF_data['mass'][:]
    CDF_charge = CDF_data['charge'][:]
    
    d = 4 ## decimal places for rounding
    tol = 1e-2
    
    other_vars = ['anum', 'znum', 'nion', 'mass', 'charge']
    
    other_vars_equal = all([(all(np.unique(np.atleast_1d(ASCOT5_dataset.data.plasma.active.read()[key]).squeeze() - np.atleast_1d(PS_data[key]).squeeze() < tol)), \
                             all(np.unique( np.atleast_1d(PS_data[key]).squeeze() -  np.atleast_1d(CDF_data[key]).squeeze() < tol))) for key in other_vars])

    print("Other values match within "+str(tol)+":", other_vars_equal)
    
    
    if plotyn == 1:
        plt.figure()
        plt.plot(A5_rho_pol, A5_Te, label='ASCOT5', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Te, label='Plasma state', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Te, label='CDF', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Electron temperature [eV]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_Ti, label='ASCOT5', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Ti, label='Plasma state', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Ti, label='CDF', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Ion temperature [eV]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_Ne, label='ASCOT5', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Ne, label='Plasma state', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Ne, label='CDF', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Electron density [m$^{-3}$]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_Ni[:,0], label='ASCOT5, species 1', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Ni[:,0], label='Plasma state, species 1', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Ni[:,0], label='CDF, species 1', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Ion density [m$^{-3}$]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_Ni[:,1], label='ASCOT5, species 2', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Ni[:,1], label='Plasma state, species 2', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Ni[:,1], label='CDF, species 2', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Ion density [m$^{-3}$]')
        plt.legend()
        
        
        plt.figure()
        plt.plot(A5_rho_pol, A5_Ni[:,2], label='ASCOT5, species 3', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Ni[:,2], label='Plasma state, species 3', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Ni[:,2], label='CDF, species 3', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Ion density [m$^{-3}$]')
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
    
    ASCOT5_filename = "../runs/NSTX_134020_GC_new_tests_6.h5"
    ASCOT5_data = Ascot(ASCOT5_filename)

    comparison = compare_plasma(ASCOT5_dataset=ASCOT5_data, TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, \
                                TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, time=None, extrapolate=True, \
                                extrapolation_type='dummy', dummy_value=5, plotyn=plotyn)
    
    if plotyn == 1:
        
        
        plt.show()

