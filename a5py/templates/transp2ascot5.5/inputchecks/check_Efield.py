# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:31:12 2023

@author: hvkg84

This script sanity checks the Efield in the ASCOT5 input file.
Checks:
    1) Is the Efield the same in the ASCOT5 file as in the TRANSP plasma state file and the full CDF file?
    => Comparison: Efield vs rho

"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.interpolate
import inputhelpers.get_Efield as get_Efield
from a5py import Ascot


def compare_Efield(ASCOT5_dataset, TRANSP_plasma_state_dataset, TRANSP_full_CDF_dataset, time=None, plotyn=0):
    """
    Start by getting the correct time value for the CDF data. It can also be found from the plasma state file.
    """
    if time == None:
        time = TRANSP_plasma_state_dataset.variables['t1'][:]
    CDF_time_idx = np.max(np.where( TRANSP_full_CDF_dataset.variables['TIME'][:] <= time ))
    
    ### FIRST get the efield values from the ASCOT5 file
    A5_Efield = -1* np.array([ASCOT5_dataset.data.efield.active.read()['dvdrho']]).squeeze() / ASCOT5_dataset.data.efield.active.read()['reff'].squeeze()
    A5_rho_min = ASCOT5_dataset.data.efield.active.read()['rhomin'].squeeze()
    A5_rho_max = ASCOT5_dataset.data.efield.active.read()['rhomax'].squeeze()
    A5_nrho  = ASCOT5_dataset.data.efield.active.read()['nrho']#.squeeze()
    A5_rho_pol = np.linspace(A5_rho_min, A5_rho_max, A5_nrho)
    
    
    ### THEN get the efield values from the plasma state file
    PS_data = get_Efield.efield_from_plasma_state_dataset(TRANSP_plasma_state_dataset, flux_averaged=False, plotyn=0)
    PS_Efield = -1* np.array(PS_data['dvdrho']) / PS_data['reff']
    PS_rho_min = PS_data['rhomin']
    PS_rho_max = PS_data['rhomax']
    PS_nrho  = PS_data['nrho']
    PS_rho_pol = np.linspace(PS_rho_min, PS_rho_max, PS_nrho)
    
    
    ### FINALLY get the efield values from the full CDF file in two different ways.
    ## (1)
    CDF_data = get_Efield.efield_from_plasma_state_dataset(TRANSP_plasma_state_dataset, flux_averaged=False, plotyn=0)
    CDF_Efield = -1* np.array(CDF_data['dvdrho']) / CDF_data['reff']
    CDF_rho_min = CDF_data['rhomin']
    CDF_rho_max = CDF_data['rhomax']
    CDF_nrho  = CDF_data['nrho']
    CDF_rho_pol = np.linspace(CDF_rho_min, CDF_rho_max, CDF_nrho)
    
    
    ## (2)
    TFCD = TRANSP_full_CDF_dataset
    # Get variables from the netCDF file
    PLFLX = TFCD.variables['PLFLX'][CDF_time_idx] # A function of XB
    rho_pol_edges_vs_rho_tor_edges = np.append(0, np.sqrt((PLFLX - 0)/(PLFLX[-1]  - 0))) ## NOTE that we are assuming psi on-axis to be ZERO....
    XB = TFCD.variables['XB'][CDF_time_idx] # Zone edges
    XB = np.append(0, XB)  # Prepend the magnetic axis
    # ERTOT   NC radial E Field, V/CM  RMAJM
    ERTOT  = TFCD.variables['ERTOT'][CDF_time_idx] * 1e2 # CONVERT to volts per metre
    XILMP = TFCD.variables['XILMP'][CDF_time_idx,:] # THE corresponding toroidal flux values.
    RMAJM = TFCD.variables['RMAJM'][CDF_time_idx,:] * 1e-2 #CONVERT TO m      # R-values for data
    
    ## Now get the indices for the outer midplane (omp)
    RAXIS  = TFCD.variables['RAXIS'][CDF_time_idx] * 1e-2 # CONVERT to metres
    omp_indices = np.where(np.round(RMAJM, 6) >= np.round(RAXIS, 6))
    rho_pol_vs_XILMP = np.interp(x=XILMP[omp_indices], xp=XB, fp=rho_pol_edges_vs_rho_tor_edges) ## Interpolate just in case the flux label values in XB are not all contained in XILMP.
    
    CDF_Efield_2 = ERTOT[omp_indices]
    CDF_rho_pol_2 = rho_pol_vs_XILMP
    
    
    if plotyn == 1:
        plt.figure()
        plt.plot(A5_rho_pol, A5_Efield, label='ASCOT5', marker='o', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_rho_pol, PS_Efield, label='Plasma state', marker='+', markersize=10, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_rho_pol, CDF_Efield, label='CDF', marker='x', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.plot(CDF_rho_pol_2, CDF_Efield_2, label='CDF, directly from file', marker='D', markersize=10, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Radial E field [V/m]')
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

    comparison = compare_Efield(ASCOT5_dataset=ASCOT5_data, TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, \
                                 TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, plotyn=plotyn)
    
    
    if plotyn == 1:
        
        
        plt.show()



