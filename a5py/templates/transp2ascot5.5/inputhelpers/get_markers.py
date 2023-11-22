# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 18:15:26 2023

@author: hvkg84

The purpose of this script is to load a NUBEAM FI birth marker set and copy the data
to an ASCOT5 hdf5 file.


ASCOT5 stores marker data in a dictionary with the following categories/names:
    ['anum', 'charge', 'energy', 'mass', 'n', 'phi', 'pitch', 'r', 'time', 'weight', 'z', 'zeta', 'znum', 'ids']
    
    'anum' is the atomic mass number, etc. 'phi' is the TOROIDAL ANGLE in DEGREES, and 'zeta' is the GYRO ANGLE in RADIANS.

NUBEAM stores marker birth data in a similar dictionary, with the following categories/names:
    ['mclabel', 'bs_r_D_MCBEAM', 'bs_z_D_MCBEAM', 'bs_rgc_D_MCBEAM', 'bs_zgc_D_MCBEAM', 'bs_xksid_D_MCBEAM', 
     'bs_einj_D_MCBEAM', 'bs_wght_D_MCBEAM', 'bs_zeta_D_MCBEAM', 'bs_time_D_MCBEAM', 'bs_ib_D_MCBEAM']

    'mclabel': the name.
    'bs_r_D_MCBEAM', 'bs_z_D_MCBEAM': the DEPOSITION r and z values IN CM.
    'bs_rgc_D_MCBEAM', 'bs_zgc_D_MCBEAM': the GC r and z values IN CM.
    'bs_xksid_D_MCBEAM': the PITCH, v_par / v, DEFINED RELATIVE TO THE PLASMA CURRENT, NOT THE BFIELD.
    'bs_einj_D_MCBEAM': the ENERGY in eV.
    'bs_wght_D_MCBEAM': the marker weights.
    'bs_zeta_D_MCBEAM': the TOROIDAL ANGLE at DEPOSITION, in DEGREES.
    'bs_time_D_MCBEAM': the TIME at deposition, in seconds.
    'bs_ib_D_MCBEAM': the ID of the NBI injector.
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import copy
import os

from os import listdir
from os.path import isfile, join
from scipy.constants import physical_constants as const

def GC_markers_from_birth_dataset(TRANSP_FI_birth_dataset, sigma_I_times_sigma_B0 = -1, AVGTIM=None, plotyn=0):
    
    TFBD = TRANSP_FI_birth_dataset
    # THESE must be returned as normal arrays rather than masked arrays
    #CHANGE THEM FROM CENTIMETRES TO METRES
    r_vals = np.array(TFBD.variables['bs_r_D_MCBEAM'][:])/100 #Now in metres
    z_vals = np.array(TFBD.variables['bs_z_D_MCBEAM'][:])/100
    r_GC_vals = np.array(TFBD.variables['bs_rgc_D_MCBEAM'][:])/100 #Now in metres
    z_GC_vals = np.array(TFBD.variables['bs_zgc_D_MCBEAM'][:])/100
    weight_vals = np.array(TFBD.variables['bs_wght_D_MCBEAM'][:])
    energy_vals =  np.array(TFBD.variables['bs_einj_D_MCBEAM'][:]) #In eV
    pitch_vals  =  np.array(TFBD.variables['bs_xksid_D_MCBEAM'][:])
    phi_vals    =  np.array(TFBD.variables['bs_zeta_D_MCBEAM'][:])
    time_vals  =  np.array(TFBD.variables['bs_time_D_MCBEAM'][:]) #Time should be in seconds for both
    time_interval = np.max(time_vals) - np.min(time_vals) if AVGTIM == None else AVGTIM # If AVGTIM has not been provided, roughly get the time interval from the markers.
    
    # Re-express as variables ASCOT5 uses
    anum   = 2 * np.ones_like(r_vals)
    charge = 1 * np.ones_like(r_vals)
    energy = energy_vals
    mass   = const["deuteron relative atomic mass"][0] * np.ones_like(r_vals)
    n      = np.size(r_vals)
    phi    = phi_vals
    pitch  = pitch_vals * sigma_I_times_sigma_B0 # Adjust for NUBEAM'S pitch CONVENTION where pitch = sigma_I_times_sigma_B0 * (v dot B)/(|v|*|B|)
    r      = r_vals
    rgc    = r_GC_vals
    time   = np.zeros_like(time_vals)
    weight = weight_vals
    z      = z_vals
    zgc    = z_GC_vals
    znum   = 1 * np.ones_like(r_vals) #znum array; Deuterium
    ids    = np.linspace(start=1, stop= n, num= n)
    
    zeta   = 2*np.pi*np.random.random(r_vals.shape) # create a random gyrophase for each particle
    
    #Scale weights so that the full power is represented for ASCOT5
    weight = np.array(weight)/time_interval
    
    #Now create the dictionary for ASCOT5
    ASCOT5_full_dictionary = {'anum':anum, 'charge':charge, 'energy':energy, 'mass':mass, 'n':n, 'phi':phi, \
                         'pitch':pitch, 'r':r_GC_vals, 'time':time, 'weight':weight, 'z':z_GC_vals, \
                         'zeta':zeta, 'znum':znum, 'ids':ids}
    return ASCOT5_full_dictionary


def GC_markers_from_multiple_birth_datasets(directory, sigma_I_times_sigma_B0 = -1, outtim_indices=None, marker_range=None, AVGTIM=None, plotyn=0):
    """
    This function automatically takes all of the birth files in the specified directory, runs the function to retrieve an
     ASCOT5 marker dictionary from a given file, and returns an ASCOT5 marker dictionary with all of the markers from each 
     of the runs.
     
    The variable 'outtim_indices' can be set to an array or tuple to specify only some outtim indices for which to collect the marker info.
    
    The filenames take the form: 134020D08_birth.cdf3
    
    ASCOT5 stores marker data in a dictionary with the following categories/names:
        ['anum', 'charge', 'energy', 'mass', 'n', 'phi', 'pitch', 'r', 'time', 'weight', 'z', 'zeta', 'znum', 'ids']
    
    The weight data needs to be modified to account for the fact that n times as many particles are now being used; 
     we must divide the weights by the number of time windows.
    """
    TRANSP_run_directory = directory
    birth_files = [f for f in listdir(TRANSP_run_directory) if isfile(join(TRANSP_run_directory, f)) and f[:-1].endswith(".cdf") and '_birth' in f]
    
    if outtim_indices != None:
        birth_files = [file for file in birth_files if file[-1] in str(np.atleast_1d(outtim_indices))]
    
    sorted_birth_files = sorted(birth_files)
    ASCOT5_full_dictionary = {}
    
    
    # Now add all of the marker data dictionaries into a single dictionary
    for file in sorted_birth_files:
        TRANSP_data = nc.Dataset(join(TRANSP_run_directory, file))
        ASCOT5_marker_dictionary = GC_markers_from_birth_dataset(TRANSP_FI_birth_dataset=TRANSP_data, sigma_I_times_sigma_B0 = sigma_I_times_sigma_B0, AVGTIM=AVGTIM)
        
        for key in ASCOT5_marker_dictionary:
            # If the full dictionary already contains this key, append the current dictionary's data to it, otherwise add the key and the data.
            if key in ASCOT5_full_dictionary:
                ASCOT5_full_dictionary[key] = np.append(ASCOT5_full_dictionary[key], ASCOT5_marker_dictionary[key])
            else:
                ASCOT5_full_dictionary[key] = ASCOT5_marker_dictionary[key]
    
        
    
    # Now correct the 'n', 'weight' and 'ids' values 
    ASCOT5_full_dictionary['n'] = np.sum(ASCOT5_full_dictionary['n'])
    ASCOT5_full_dictionary['ids'] = np.linspace(start=1, stop=ASCOT5_full_dictionary['n'], num=ASCOT5_full_dictionary['n']).astype(int)    
    ASCOT5_full_dictionary['weight'] = ASCOT5_full_dictionary['weight']/len(sorted_birth_files)
    
    if marker_range != None:
        # If a range of markers to include has been specified, e.g only the first 1k particles, then return only that range.
        for key in ASCOT5_full_dictionary:
            if key != 'n':
                ASCOT5_full_dictionary[key] =  ASCOT5_full_dictionary[key][marker_range[0]:marker_range[1]]
        # Now adjust the total number of particles.
        ASCOT5_full_dictionary['n'] = len(ASCOT5_full_dictionary['ids'])
        
    return ASCOT5_full_dictionary


if __name__ == '__main__':
    
    directory_with_TRANSP_birth_files = '../../134020D22'
    
    ASCOT5_markers_from_NUBEAM = GC_markers_from_multiple_birth_datasets(directory=directory_with_TRANSP_birth_files, sigma_I_times_sigma_B0 = -1, \
                                                                                   outtim_indices=1, marker_range=[None, None], plotyn=0)
    
    #Use ASCOT5's output-writing function to write the new set of markers
    mrk_gc.write_hdf5(ASCOT5_filename, **ASCOT5_markers_from_NUBEAM, desc='GC_from_1st_NUBEAM_file') 
    print('**SUCCESS**: markers loaded successfully into file ' + ASCOT5_filename + '.')
