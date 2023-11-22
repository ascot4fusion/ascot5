# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:31:12 2023

@author: Ian Dolby

This script sanity checks the markers in the ASCOT5 input file.
Checks:
    1) Are the marker values the same in the ASCOT5 file as in the TRANSP birth file? (adjusted for power)
    2) Do the implied power values match which the information in the TRANSP plasma state file and the full CDF file?
    3) Are the marker Larmor radii reasonable? Plot the cumulative histogram of the Larmor radii.

The comparisons for checking are:
    1) Are the values all identical to 3 significant figures? (after adjustements for power)
    2) Sum of energy * weight / time_window from each injector
    3) Plot the cumulative histogram of the Larmor radii, using the bfield in the A5 file.

NOTE: The user should know what power values are expected.
"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.interpolate
import inputhelpers.get_markers as get_markers
from a5py import Ascot


def compare_markers(ASCOT5_dataset, directory_with_TRANSP_birth_files, sigma_I_times_sigma_B0 = -1, outtim_indices=1, marker_range=[None, None], plotyn=0):
    """
    ASCOT5_full_dictionary = {'anum':anum, 'charge':charge, 'energy':energy, 'mass':mass, 'n':n, 'phi':phi, \
                         'pitch':pitch, 'r':r_GC_vals, 'time':time, 'weight':weight, 'z':z_GC_vals, \
                         'zeta':zeta, 'znum':znum, 'ids':ids}
    """
    # A5_data = ASCOT5_dataset.marker.active.read()
    A5_data = ASCOT5_dataset.data.marker.active.read()
    birth_data = get_markers.GC_markers_from_multiple_birth_datasets(directory=directory_with_TRANSP_birth_files, \
                                                                     sigma_I_times_sigma_B0 = sigma_I_times_sigma_B0, \
                                                                     outtim_indices=outtim_indices, marker_range=marker_range, plotyn=plotyn)
    A5_data_copy = {**A5_data}
    del(A5_data_copy['zeta'])
    vars_to_check = A5_data_copy.keys()
    tol = 1e-2
    vars_equal = all([(all(np.unique(np.atleast_1d(A5_data[key]).squeeze() - np.atleast_1d(birth_data[key]).squeeze() < tol))) for key in vars_to_check])
    print("All values except randomised gyrophase match within "+str(tol)+" (absolute difference)?", vars_equal)

    return vars_equal



if __name__ == '__main__':
    
    directory_with_TRANSP_birth_files = '../../134020D30'
    
    
    ASCOT5_filename = "../runs/NSTX_134020_GC_new_tests_8.h5"
    ASCOT5_data = Ascot(ASCOT5_filename)
    
    
    compare_markers(ASCOT5_dataset=ASCOT5_data, directory_with_TRANSP_birth_files=directory_with_TRANSP_birth_files, \
                    sigma_I_times_sigma_B0 = -1, outtim_indices=1, marker_range=[None, None], plotyn=0)

