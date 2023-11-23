#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ian.j.dolby@durham.ac.uk

This script sets up an ASCOT5 run from TRANSP files.
It requires:
    A) A TRANSP plasma state file
    B) A NUBEAM marker birth file
    C) If running input checks, a TRANSP full CDF file

"""
import sys

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

import inputhelpers.get_Bfield as get_Bfield
import inputhelpers.get_Efield as get_Efield
import inputhelpers.get_plasma as get_plasma
import inputhelpers.get_neutrals as get_neutrals
import inputhelpers.get_wall as get_wall
import inputhelpers.get_markers as get_markers

import inputchecks.check_Bfield as check_Bfield
import inputchecks.check_Efield as check_Efield
import inputchecks.check_plasma as check_plasma
import inputchecks.check_neutrals as check_neutrals
import inputchecks.check_wall as check_wall
import inputchecks.check_markers as check_markers


from a5py import Ascot
from a5py.ascot5io.options import Opt
from scipy.constants import physical_constants as const



def make_ascot5_slowingdownrun(a5, TRANSP_marker_directory=None, TRANSP_plasma_state_dataset=None, TRANSP_full_CDF_dataset=None, marker_outtim=1, check_input=True, plotyn=0, checking_plotyn=1):
    """==================================== GET INPUT FROM TRANSP FILES AND WRITE TO ASCOT5 FILE ============================================
    """
    
    TPSD = TRANSP_plasma_state_dataset
    TFCD = TRANSP_full_CDF_dataset
    
    ### Bfield
    Bfield_dict = get_Bfield.bfield_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TPSD, COCOS_in=5, COCOS_out=3)
    a5.data.create_input("B_2DS", **Bfield_dict, desc="from TRANSP")
    ### Efield
    Efield_dict = get_Efield.efield_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TPSD)
    a5.data.create_input("E_1DS", **Efield_dict)
    ### Plasma
    plasma_dict = get_plasma.plasma_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TPSD, extrapolate=True, extrapolation_type='dummy', dummy_value=1)
    a5.data.create_input("plasma_1D", **plasma_dict)
    ### Neutrals
    neutrals_dict = get_neutrals.neutrals_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TPSD, maxwellian=[1,1], extrapolate=True, \
                                                                    extrapolation_type='constant')
    a5.data.create_input("N0_1D", **neutrals_dict)
    ### Wall 
    wall_dict = get_wall.wall_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TPSD)
    a5.data.create_input("wall_2D", **wall_dict)
    ### Markers, given the directory containing the birth files and the possibility of selecting a specific outtim.
    marker_dict = get_markers.GC_markers_from_multiple_birth_datasets(directory=TRANSP_marker_directory, sigma_I_times_sigma_B0 = -1, \
                                                                      outtim_indices=marker_outtim, plotyn=plotyn)
    a5.data.create_input("gc", **marker_dict, desc='GC from NUBEAM birth file ' + str(marker_outtim))
    
    ### Dummy input
    a5.data.create_input("Boozer")
    a5.data.create_input("MHD_STAT")
    a5.data.create_input("asigma_loc")
    
    """======================================== CHECK INPUT MATCHES CONTENT OF TRANSP FILES ==========================================
    """
    if check_input == True:
        checking_plotyn = checking_plotyn
        ## Bfield_check
        print("Checking Bfield.")
        check_Bfield.compare_bfield(ASCOT5_dataset=a5, TRANSP_plasma_state_dataset=TPSD, TRANSP_full_CDF_dataset=TFCD, plotyn=checking_plotyn)
        ## Efield_check
        print("Checking Efield.")
        check_Efield.compare_Efield(ASCOT5_dataset=a5, TRANSP_plasma_state_dataset=TPSD, TRANSP_full_CDF_dataset=TFCD, plotyn=checking_plotyn)
        ## plasma check
        print("Checking plasma.")
        check_plasma.compare_plasma(ASCOT5_dataset=a5, TRANSP_plasma_state_dataset=TPSD, TRANSP_full_CDF_dataset=TFCD, extrapolate=True, \
                                    extrapolation_type='dummy', dummy_value=1, plotyn=checking_plotyn)
        ## neutrals check
        print("Checking neutrals.")
        check_neutrals.compare_neutrals(ASCOT5_dataset=a5, TRANSP_plasma_state_dataset=TPSD, TRANSP_full_CDF_dataset=TFCD, extrapolate=True, \
                                    extrapolation_type='constant', plotyn=checking_plotyn)
        ## wall check
        print("Checking wall.")
        check_wall.compare_wall(ASCOT5_dataset=a5, TRANSP_plasma_state_dataset=TPSD, TRANSP_full_CDF_dataset=TFCD, plotyn=checking_plotyn)
        ## marker check
        print("Checking markers:")
        check_markers.compare_markers(ASCOT5_dataset=a5, directory_with_TRANSP_birth_files=TRANSP_marker_directory, \
                        sigma_I_times_sigma_B0 = -1, outtim_indices=marker_outtim, marker_range=[None, None], plotyn=checking_plotyn)
        print("End of checks.")
        
        
    
    
    #################### #################### ####################   OPTIONS   #################### ################### ####################    
    m_d = const['deuteron mass'][0]
    
    settings={}
    # Do you want 2D or 3D wall model
    settings["wall_make_3D"] =False
    # Simulate guiding centers (GC) (if false, simulate gyro-orbits (GO))
    settings["sim_gc_mode"] = False
    # If this is a guiding center simulation, do we use adaptive step?
    settings["sim_use_adaptivestep"] = True
    # Try out the hybrid model for wall collision checks
    settings["sim_use_hybrid"] = True
    # Record GC position even if the simulation is GO simulation
    settings["sim_recordGCasGO"] = False
    # Magnetic field input, choose whether you want
    # - Analytical or spline interpolated
    # - Axisymmetric or 3D
    settings["bfield_make_3D"] = False
    # Set options to null state
    o = {}
    o["SIM_MODE"]               = 2
    o["FIXEDSTEP_USERDEFINED"]  = 1e-9
    
    if settings["sim_gc_mode"]:
        o["SIM_MODE"]              = 2
        o["FIXEDSTEP_USERDEFINED"] = 1e-7
        

    if settings["sim_use_hybrid"]:
        o["SIM_MODE"]               = 3
        o["FIXEDSTEP_USERDEFINED"]  = 5.0e-8

    if settings["sim_use_adaptivestep"]:
        o["ENABLE_ADAPTIVE"] = 1

    if settings["sim_recordGCasGO"]:
        o["RECORD_MODE"] = 0

    o["FIXEDSTEP_USE_USERDEFINED"] = 0
    o["ADAPTIVE_TOL_ORBIT"]        = 1.0e-9
    o["ADAPTIVE_TOL_CCOL"]         = 1.0e-4
    o["ADAPTIVE_MAX_DRHO"]         = 1.0
    o["ADAPTIVE_MAX_DPHI"]         = 1.0

    o["ENDCOND_SIMTIMELIM"] = 1
    o["ENDCOND_CPUTIMELIM"] = 1
    o["ENDCOND_RHOLIM"]     = 0
    o["ENDCOND_ENERGYLIM"]  = 1
    o["ENDCOND_WALLHIT"]    = 1
    
    #MAXSIMTIME divided by FIXEDSTEP_USERDEFINED GIVES ME the recorded orbit points

    o["ENDCOND_LIM_SIMTIME"]  = 0.15
    o["ENDCOND_MAX_MILEAGE"]  = 1
    o["ENDCOND_MAX_CPUTIME"]  = 1e6
    o["ENDCOND_MAX_RHO"]      = 1
    o["ENDCOND_MIN_RHO"]      = 0
    o["ENDCOND_MIN_ENERGY"]   = 10
    o["ENDCOND_MIN_THERMAL"]  = 1.5
    
   
    o["ENABLE_ORBIT_FOLLOWING"]    = 1
    o["ENABLE_COULOMB_COLLISIONS"] = 1


    # All distributions on
    o["ENABLE_DIST_5D"]    = 1
    o["ENABLE_DIST_6D"]    = 0
    o["ENABLE_DIST_RHO5D"] = 0
    o["ENABLE_DIST_RHO6D"] = 0

    o["DIST_MIN_R"]    = 0.18
    o["DIST_MAX_R"]    = 1.6
    o["DIST_NBIN_R"]   = 50
    
    o["DIST_MIN_PHI"]  = 0
    o["DIST_MAX_PHI"]  = 360
    o["DIST_NBIN_PHI"] = 1

    o["DIST_MIN_Z"]    = -1.7
    o["DIST_MAX_Z"]    = 1.7
    o["DIST_NBIN_Z"]   = 50

    o["DIST_MIN_RHO"]  = 0
    o["DIST_MAX_RHO"]  = 1.5
    o["DIST_NBIN_RHO"] = 50

    o["DIST_MIN_THETA"]  = 0
    o["DIST_MAX_THETA"]  = 360
    o["DIST_NBIN_THETA"] = 10

    o["DIST_MIN_PPA"]  = -3e6 * m_d
    o["DIST_MAX_PPA"]  = 3e6 * m_d
    o["DIST_NBIN_PPA"] = 50

    o["DIST_MIN_PPE"]  = 0
    o["DIST_MAX_PPE"]  = 3e6 * m_d
    o["DIST_NBIN_PPE"] = 50

    o["DIST_MIN_PR"]    = -3e6 * m_d
    o["DIST_MAX_PR"]    = 3e6 * m_d
    o["DIST_NBIN_PR"]   = 50

    o["DIST_MIN_PPHI"]  = -3e6 * m_d
    o["DIST_MAX_PPHI"]  = 3e6 * m_d
    o["DIST_NBIN_PPHI"] = 50

    o["DIST_MIN_PZ"]    = -3e6 * m_d
    o["DIST_MAX_PZ"]    = 3e6 * m_d
    o["DIST_NBIN_PZ"]   = 50

    o["DIST_MIN_TIME"]    = 0
    o["DIST_MAX_TIME"]    = 1
    o["DIST_NBIN_TIME"]   = 1
    
    o["DIST_MIN_CHARGE"]    = -100
    o["DIST_MAX_CHARGE"]    = 100
    o["DIST_NBIN_CHARGE"]   = 1

    o["ENABLE_ORBITWRITE"]    = 0
    o["ORBITWRITE_MODE"]      = 1
    o["ORBITWRITE_NPOINT"]    = 1e3
    o["ORBITWRITE_INTERVAL"]  = 0
    
    opt = Opt.get_default()
    opt.update(**o)
    a5.data.create_input("opt", **opt)
    



if __name__ == '__main__':
    run_ID = '134020D31'
    TRANSP_CDF_directory = '../../../../ASCOT_files/134020D31/'
    
    TRANSP_plasma_state_filename = TRANSP_CDF_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    TRANSP_full_CDF_filename = TRANSP_CDF_directory + run_ID + ".CDF"
    TRANSP_full_CDF_data = nc.Dataset(TRANSP_full_CDF_filename)
    
    
    filename_to_create = "./runs/NSTX_134020_D31_GC_test.h5"
    
    ASCOT5_object = Ascot("./" + filename_to_create, create=True)
    
    make_ascot5_slowingdownrun(a5 = ASCOT5_object, TRANSP_marker_directory=TRANSP_CDF_directory, TRANSP_plasma_state_dataset=TRANSP_data, \
                               TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, check_input=True, checking_plotyn=1)
    print("** SUCCESS**: File " + filename_to_create + " generated successfully.")
