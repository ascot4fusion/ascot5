# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:31:12 2023

@author: Ian Dolby

This script sanity checks the wall in the ASCOT5 input file.
Checks:
    1) Is the wall the same in the ASCOT5 file as in the TRANSP plasma state file and the full CDF file?

This script merely plots the limiters from each dataset.
"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from a5py import Ascot

def compare_wall(ASCOT5_dataset, TRANSP_plasma_state_dataset, TRANSP_full_CDF_dataset, time=None, plotyn=0):
    """
    Start by getting the correct time value for the CDF data (although, the limiter should not be time-dependent). 
        It can also be found from the plasma state file.
    """
    if time == None:
        time = TRANSP_plasma_state_dataset.variables['t1'][:]
    CDF_time_idx = np.max(np.where( TRANSP_full_CDF_dataset.variables['TIME'][:] <= time ))
    
    A5_wall_R = np.array([ASCOT5_dataset.data.wall.active.read()['r']]).squeeze()
    A5_wall_Z = np.array([ASCOT5_dataset.data.wall.active.read()['z']]).squeeze()
    PS_wall_R = np.array(TRANSP_plasma_state_dataset.variables['rlim'][:])
    PS_wall_Z = np.array(TRANSP_plasma_state_dataset.variables['zlim'][:])
    CDF_wall_R = np.array(TRANSP_full_CDF_dataset.variables['RLIM'][CDF_time_idx]) * 1e-2 #CONVERT to m
    CDF_wall_Z = np.array(TRANSP_full_CDF_dataset.variables['YLIM'][CDF_time_idx]) * 1e-2 #CONVERT to m
    
    if plotyn == 1:
        plt.figure()
        plt.plot(A5_wall_R, A5_wall_Z, label='ASCOT5 wall', marker='o', markersize=13, markeredgewidth=2, fillstyle='none')
        plt.plot(PS_wall_R, PS_wall_Z, label='Plasma state wall', marker='+', markersize=13, markeredgewidth=2, fillstyle='none')
        plt.plot(CDF_wall_R, CDF_wall_Z, label='CDF wall', marker='x', markersize=13, markeredgewidth=1, fillstyle='none')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
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

    comparison = compare_wall(ASCOT5_dataset=ASCOT5_data, TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, \
                                 TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, plotyn=plotyn)
    
    
    if plotyn == 1:
        
        
        plt.show()
