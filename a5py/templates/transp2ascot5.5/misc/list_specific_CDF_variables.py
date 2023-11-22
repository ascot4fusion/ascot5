# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:44:26 2023

@author: hvkg84

This script generates a list of variables, whose names or descriptions contain certain key phrases, from the TRANSP CDF files,
 and outputs the list to a csv file.
This is mainly useful for verifying that a given feature is on or off.

FEATURES:
    
    CX: Setting nlbeamcx to false should eliminate any CX component of the FIs slowing down. This will not eliminate halo neutrals, or
         CX in the bulk plasma (switching these off would be ridiculous).
    Rotation:
    FLR/GFLR:


"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

########## INPUT FILES ##########
plotyn = 1
saveyn = 0

run_ID = '134020D30'
TRANSP_plasma_state_directory = '....//134020D30/'

TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)

TRANSP_full_CDF_filename = TRANSP_plasma_state_directory + run_ID + ".CDF"
TRANSP_full_CDF_data = nc.Dataset(TRANSP_full_CDF_filename)

TPSDV = TRANSP_plasma_state_data.variables
TCFDV = TRANSP_full_CDF_data.variables

time = 0.65
time_idx = np.min(np.where(TCFDV['TIME'][:] > time))

########## TERMS TO LOOK FOR ##########
# search_terms = ['TIME']
search_terms = ['WALL', 'LIM', 'BOUN']
# search_terms = ['CX', 'HALO']
# search_terms = ['FLR', 'Larmor']

########## GENERATE THE DICTIONARY OF VARIABLES WHICH HAVE ACTUAL DATA ##########

full_CDF_variables = [ key for key in TCFDV.keys() if (hasattr(TCFDV[key], 'units') and hasattr(TCFDV[key], 'long_name'))]
plasma_state_variables = [ key for key in TPSDV.keys() if (hasattr(TPSDV[key], 'units') and hasattr(TPSDV[key], 'long_name'))]


########## ADD THE VARIABLES WHOSE NAME OR DESCRIPTION CONTAINS THE SEARCH TERMS, TOGETHER WITH VALUES, UNITS, ETC ##########

full_CDF_search_result = np.array([ [key, np.max(TCFDV[key][time_idx]), np.min(TCFDV[key][time_idx]), TCFDV[key].units.strip(), \
                            TCFDV[key].long_name.strip(), np.max(np.abs(TCFDV[key][time_idx])) == 0] for key in full_CDF_variables \
                            if any((term in str(key).upper() or term in str(TCFDV[key].long_name).upper()) for term in search_terms) ])

plasma_state_search_result = np.array([ [key, np.max(TPSDV[key]), np.min(TPSDV[key]), TPSDV[key].units.strip(), \
                            TPSDV[key].long_name.strip(), np.max(np.abs(TPSDV[key])) == 0] for key in plasma_state_variables \
                            if any((term in str(key).upper() or term in str(TPSDV[key].long_name).upper()) for term in search_terms) ])

########## SORT THE ARRAYS BY VARIABLE NAME ##########
full_CDF_search_result_argsort = np.argsort(np.char.upper(full_CDF_search_result[:,0]))
plasma_state_search_result_argsort = np.argsort(np.char.upper(plasma_state_search_result[:,0]))
full_CDF_search_result_sorted = full_CDF_search_result[full_CDF_search_result_argsort]
plasma_state_search_result_sorted = plasma_state_search_result[plasma_state_search_result_argsort]

########## ADD THE COLUMN DESCRIPTIONS ##########
variable_descriptions = np.array([ 'Variable', 'Max value', 'Min value', 'Units', 'Description', 'Equals zero' ])

full_CDF_search_result_full = np.vstack((variable_descriptions, full_CDF_search_result_sorted))
plasma_state_search_result_full = np.vstack((variable_descriptions, plasma_state_search_result_sorted))

########## SAVE RESULTS TO CSV ##########
if saveyn == 1:
    np.savetxt( run_ID + '_'.join(search_terms) + '_vars_in_full_CDF.csv', full_CDF_search_result_full, delimiter=',', fmt='%s' )
    np.savetxt( run_ID + '_'.join(search_terms) + '_vars_in_plasma_state.csv', plasma_state_search_result_full, delimiter=',', fmt='%s' )





