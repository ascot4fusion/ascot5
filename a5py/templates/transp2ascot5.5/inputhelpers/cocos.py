# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:36:29 2023

@author: Ian Dolby

Below are functions which take the TRANSP plasma state output file (FI_OUTTIM must be set in the namelist)
 and deduce the COCOS convention used, and then transform it to whatever COCOS convention is required by the user.

TRANSP appears to use COCOS5, but it may be that someday the COCOS used in its output files can be changed by the user,
 so these functions are based entirely on variables contained in the plasma state output file.

Much of this code is based on the eqdsk2input.py script written by Matteo Vallar, which can be found 
 in the ASCOT5 GitLab repository in python/a5py/a5py/preprocessing/eqdsk2input.py .
 
The COCOS convention is laid out in the paper:
    
    https://www.sciencedirect.com/science/article/pii/S0010465512002962

    Notes: a) Table 4 is not updated. 
           b) EFIT has usually cocos=7, but if people try to keep q>0, then cocos changes.
"""

import numpy as np
import scipy as sc
import netCDF4 as nc
import matplotlib.pyplot as plt


def cocos_check(TRANSP_plasma_state_dataset, COCOS, quiet=False):
    """cocos check
    This function does not trust the user and calculates the COCOS of the input from the variable values.
    
    The COCOS variable is the COCOS integer for the input.
    """
    cocos_dict = fill_cocosdict(COCOS) # Get the signs corresponding to the input COCOS.
    TPSD = TRANSP_plasma_state_dataset
    
    ### Get sigma_RphiZ from the sign of Bphi and the direction relative to anti-clockwise.
    Bphi_ccw = TPSD.variables['kccw_Bphi'][:]
    BphiRZ   = TPSD.variables['BphiRZ'][:]
    sign_b0 = np.sign(BphiRZ[0,0])
    sigma_RphiZ = Bphi_ccw * sign_b0 ## +ve value means phi is anti-clockwise
    
    ### Get sigma_Ip from sigma_RphiZ and the direction of the current relative to anti-clockwise.
    Jphi_ccw = TPSD.variables['kccw_Jphi'][:]
    sigma_Ip = Jphi_ccw * sigma_RphiZ
    
    ### Now get sigma_Bp from sigma_Ip and the sign of psi_edge - psi_axis
    psipol = TPSD.variables['psipol'][:] # As a function of rho_eq
    psi_near_axis  = psipol[0]
    psi1  = psipol[-1]
    psi_increase = np.sign(psi1-psi_near_axis) # 1 if increasing, -1 if decreasing
    sigma_Bp = psi_increase/sigma_Ip
    
    ### Now get sigma_rhothetaphi from the flux surface dZ/dtheta value provided, at the outer midplane. (theta should always be measured from the horizontal axis)
    theta_zero_index = np.min(np.where(TPSD.variables['th_eq'][:] >= 0 ) )
    dZ_dtheta = TPSD.variables['dZ_geo_dTH'][theta_zero_index,-1]
    # If dZ_dtheta is positive, theta is anti-clockwise. If theta is anti-clockwise and phi is anticlockwise, sigma_rhothetaphi = -1,
    #  so positive dZ_dtheta and positive sigma_RphiZ implies negative sigma_rhothetaphi.
    sigma_rhothetaphi = - np.sign(dZ_dtheta) * sigma_RphiZ

    
    ### Now check for consistency with the safety factor.
    q        = TPSD.variables['q_eq'][:]
    sign_q = np.sign(q[0])
    
    if sign_q*cocos_dict['sigma_rhothetaphi']*sigma_Ip*sign_b0<0:
        print(f"The sign of q is not consistent. sign(q)={sign_q}, but it should be sigma_rhothetaphi*sigma_Ip*sign_b0={cocos_dict['sigma_rhothetaphi']}*{sigma_Ip}*{sign_b0}={cocos_dict['sigma_rhothetaphi']*sigma_Ip*sign_b0}")
    

    ### Now we find if psi has 2pi value by calculating B_R from the provided PsiRZ grid and 
    ###  checking the absolute value (which is independent of sigma_Bp) against the provided B_R values.
    BRRZ = TPSD.variables['BRRZ'][:].T ## TRANSPOSE so that the indices correspond to R,Z, and not the other way around
    PsiRZ = TPSD.variables['PsiRZ'][:].T
    R = TPSD.variables['R_grid'][:]
    Z = TPSD.variables['Z_grid'][:]
    R_diff = np.diff(R)
    Z_diff = np.diff(Z)
    R_mg, Z_mg = np.meshgrid(R, Z, indexing='ij')
    psi_prime = np.gradient(PsiRZ, R_diff[0], Z_diff[0], edge_order=2)
    BR_abs_from_grid = np.abs(1/R_mg * psi_prime[1])
    BRRZ_abs = np.abs(BRRZ)
    
    # If PsiRZ is the STREAM FUNCTION, then it has already been divided by 2pi and the diff without 2pi will be closer.
    # If it is the POLOIDAL FLUX, then we need to divide by 2pi when calculating B_R, and the diff WITH 2pi will be closer.
    mean_abs_B_R_diff_no2pi = np.mean(BRRZ_abs - BR_abs_from_grid)
    mean_abs_B_R_diff_2pi = np.mean(BRRZ_abs - BR_abs_from_grid/(2*np.pi))
    
    if np.abs(mean_abs_B_R_diff_2pi) < np.abs(mean_abs_B_R_diff_no2pi):
        cocos_gt_10 = True
    
    
    # print('Bphi_ccw', Bphi_ccw)
    # print('Jphi_ccw', Jphi_ccw)
    # print('sign_q', sign_q)
    # print('sign_b0', sign_b0)
    # print('sigma_Ip', sigma_Ip)
    # print('psi_increase', psi_increase)
    # print('sigma_Bp', sigma_Bp)
    # print('sigma_RphiZ', sigma_RphiZ)
    # print('sigma_rhothetaphi', sigma_rhothetaphi)
    
    
    cocos_read = _sign_to_cocos(sigma_Bp, sigma_RphiZ, sigma_rhothetaphi, cocos_gt_10=False)

    if COCOS not in cocos_read:
        error_msg = "=== \n"
        error_msg += f"You said cocos {COCOS}, I read cocos {cocos_read[0:2]} (depending on direction of phi)"
        error_msg += " \nWe strongly suggest to check the cocos you have.\nThis is fundamental for correct B field creation.\n"
        error_msg += "=== \n"
        print(error_msg)
        error_msg = f" COCOS in file {cocos_read} is different from the one in input [{COCOS}]"
        raise Exception(error_msg)

    if quiet==False:
        print("Good! Your cocos matches!")
    
    
    cocos_sigma_dict = {'sigma_Ip': sigma_Ip, 'sign_b0': sign_b0}
    
    return cocos_sigma_dict


def _sign_to_cocos(sigma_Bp, sigma_RphiZ, sigma_rhothetaphi, cocos_gt_10):
    """
    Associating the sign with the correct cocos
    
    Parameters:
        sigma_Bp (int): sigma_Bp as in cocos definition
        sigma_rhothetaphi (int): sigma_rhothetaphi as in cocos definition
        sigma_RphiZ (int): sigma_RphiZ as in cocos definition
        cocos_gt_10 (bool): true if the cocos should be >10
    Attributes:
        cocos_read (array): cocos inferred from input
    
    """
    cocos_read=[0,0]
    if sigma_Bp>0:
        if sigma_rhothetaphi>0:
            if sigma_RphiZ > 0:
                cocos_read = [1, 11]
            else:
                cocos_read = [2, 12]
        
        else:
            if sigma_RphiZ > 0:
                cocos_read = [5, 15]
            else:
                cocos_read = [6, 16]
    
    else:
        if sigma_rhothetaphi>0:
            if sigma_RphiZ > 0:
                cocos_read = [7, 17]
            else:
                cocos_read = [8, 18]

        else:
            if sigma_RphiZ > 0:
                cocos_read = [3, 13]
            else:
                cocos_read = [4, 14]

    cocos_read = np.array(cocos_read)
    if cocos_gt_10:
        cocos_read = cocos_read[cocos_read>10]
    else:
        cocos_read = cocos_read[cocos_read<10]
    return cocos_read



def cocos_transform(TRANSP_plasma_state_dataset, COCOS_in, COCOS_out=3, quiet=False):
    """
    This function converts the magnetic input from their starting cocos to a given cocos (needed by ascot5)


    Parameters:
        COCOS_in (int): input cocos
        COCOS_out (int): output cocos
        
    """
    TPSD = TRANSP_plasma_state_dataset
    
    if quiet==False:
        print("COCOS tranformation from "+str(COCOS_in)+" to " +str(COCOS_out))
        print('Checking specified input COCOS against the provided data')
        
    cocos_check(TPSD, COCOS=COCOS_in, quiet=quiet)
    
    cocosin = fill_cocosdict(COCOS_in)
    cocosout = fill_cocosdict(COCOS_out)


    q        = TPSD.variables['q_eq'][:]
    Ip       = TPSD.variables['curt'][-1] # In Amps
    BphiRZ   = TPSD.variables['BphiRZ'][:] # In Amps
    psi      = TPSD.variables['PsiRZ'][:]
    psipol   = TPSD.variables['psipol'][:] # As a function of rho_eq
    psi0     = psipol[0]
    psi1     = psipol[-1]
    B0_abs   = TPSD.variables['B_axis'][:] # In Amps

    sign_q = np.sign(q[0])
    sign_b0 = np.sign(BphiRZ[0,0])


    # Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
    sigma_Bp_eff = cocosin['sigma_Bp'] * cocosout['sigma_Bp']
    exp_Bp_eff   = cocosout['exp_Bp']  - cocosin['exp_Bp']
    sigma_rhothetaphi_eff = cocosin['sigma_rhothetaphi'] * cocosout['sigma_rhothetaphi']
    
    sigma_Ip_eff = cocosin['sigma_RphiZ']*cocosout['sigma_RphiZ']
    sigma_B0_eff = cocosin['sigma_RphiZ']*cocosout['sigma_RphiZ']

    ### Define input
    psirz_in   = psi
    psiaxis_in = psi0
    psiedge_in = psi1

    q_in  = q
    b0_in = B0_abs * sign_b0
    ip_in = Ip
    
    
    ### Transform
    _fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2*np.pi)**exp_Bp_eff
    psirz     = psirz_in   * _fact_psi 
    psiaxis   = psiaxis_in * _fact_psi
    psiedge   = psiedge_in * _fact_psi
    
    q  = q_in  * sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff
    b0 = b0_in * sigma_B0_eff
    ip = ip_in * sigma_Ip_eff
    
    bphi_out = BphiRZ * sigma_B0_eff

    ### Define output
    plasma_state_out = {}
    plasma_state_out['psi']   = psirz
    plasma_state_out['psi0']  = psiaxis
    plasma_state_out['psi1']  = psiedge
    
    plasma_state_out['q']     = q
    plasma_state_out['B0']    = b0
    plasma_state_out['Ip']    = ip
    plasma_state_out['bphi']  = bphi_out
    plasma_state_out['COCOS_in']  = COCOS_in
    plasma_state_out['COCOS_out']  = COCOS_out
    plasma_state_out['sigma_B0_eff']  = sigma_B0_eff ## Output these variables to allow conversion of other datasets if desired.
    plasma_state_out['sigma_Ip_eff']  = sigma_Ip_eff
    plasma_state_out['sigma_Bp_eff']  = sigma_Bp_eff
    plasma_state_out['COCOS_in_signs']  = cocosin
    plasma_state_out['COCOS_out_signs']  = cocosout
    return plasma_state_out


def fill_cocosdict(COCOS):
    """
    Function to fill the dictionary with the COCOS variables

    Parameters:
        COCOS (int): input cocos
    Args:
        cocosdict (dict): dictionary with cocos variables
    
    """
    cocos_keys = ['sigma_Bp', 'sigma_RphiZ', 'sigma_rhothetaphi',\
          'sign_q_pos', 'sign_pprime_pos', 'exp_Bp']
    cocosdict = dict.fromkeys(cocos_keys)

    cocosdict['exp_Bp'] = 1 if COCOS > 10 else 0

    if COCOS==1 or COCOS==11:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==2 or COCOS==12:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==3 or COCOS==13:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==4 or COCOS==14:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==5 or COCOS==15:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==6 or COCOS==16:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==7 or COCOS==17:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==8 or COCOS==18:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = +1
    else:
        raise ValueError(f"COCOS {COCOS} does not exist \n")

    return cocosdict

if __name__ == '__main__':
    
    plotyn=0
    
    run_ID = '134020D22'
    TRANSP_plasma_state_directory = '../../134020D22/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    cocos_check(TRANSP_plasma_state_data, COCOS=5)
    
    output_equilibrium_dict = cocos_transform(TRANSP_plasma_state_data, COCOS_in=5, COCOS_out=15)
    
    