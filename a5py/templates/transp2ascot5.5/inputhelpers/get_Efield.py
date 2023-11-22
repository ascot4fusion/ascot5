# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 18:15:34 2023

@author: Ian Dolby

Args:
        fn : str <br>
            Full path to the HDF5 file.
        nrho : int <br>
            Number of rho slots in data.
        rhomin : float <br>
            Minimum rho value.
        rhomax : float <br>
            Maximum rho value.
        dvdrho : array_like (nrho,1) <br>
            Derivative of electric potential WRT rho_poloidal [V].
            If reff = 1, this is essentially equal to dv/drho
        reff : float <br>
            Effective minor radius of the plasma [m].
            This is the reciprocal of drho_pol_dr.
        desc : str, optional <br>
            Input description.

NOTE: As of 04/11/2023, there is a small difference between PLFLX from the full CDF file and psipol from the plasma state file, which seems to lead to
                        the Efield graphs being very slightly different ~ 1% at most in the significant regions.

"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import copy
import scipy.interpolate

from scipy.constants import physical_constants as const

def efield_from_plasma_state_dataset(TRANSP_plasma_state_dataset, flux_averaged=False, plotyn=0):
    """ Since Efield = -grad(Phi_electric), Efield = - partial(Phi_electric)/partial(rho) * grad(rho), by the chain rule, since Phi_electric is a flux quantity.
    Phi_electric = 1/charge * Epot, but as Epot is in keV, we also multiply by 1e3 * 1.6e-19 to get Phi_electric in units of VOLTS.
        Therefore Phi_electric [V] = 1e3 * Epot [J/C (= V)].
    
    The method which is consistent with the data in the full CDF file is to obtain these values only for the midplane. However, 
      the "grho1" variable in the plasma state file is the FLUX SURAFCE AVERAGE grad(rho), and using the EFIELD calculated using this is NOT
       consistent with using the outer midplane R-values.
       On the midplane, the rate of change of rho with respect to Z must be around 0; therefore, on the midplane, |grad(rho)| = |partial(rho)/partial(R)|
       The difference is roughly a scale factor of around 1.4.
    """   
    TPSD = TRANSP_plasma_state_dataset
    
    ## START by getting rho_pol at the rho_tor values provided. This is only for output to ASCOT5; the calculations use rho_tor.
    rho_tor_edges = TPSD.variables['rho_eq'][:]
    psipol = TPSD.variables['psipol'][:]  # Get the poloidal flux values which correspond to the rho_toroidal values (the profiles are provided in terms of rho_toroidal)
    psi_interp = scipy.interpolate.RegularGridInterpolator(points=(TPSD.variables['R_grid'][:], TPSD.variables['Z_grid'][:]), values=TPSD.variables['PsiRZ'][:].T, method='cubic')
    psi0  = psi_interp((TPSD.variables['R_axis'][:], TPSD.variables['Z_axis'][:]))
    psi1  = psipol[-1]
    rho_pol_edges_vs_rho_tor_edges = np.sqrt((psipol - psi0)/(psi1 - psi0))
    
    
    Epot_radial = TPSD.variables['Epot'][:]  # The electric potential ENERGY in keV at rho_tor_centres
    dEpot_drho_pol = np.gradient(Epot_radial, rho_pol_edges_vs_rho_tor_edges, edge_order=2)
    dPhi_electric_drho_pol = dEpot_drho_pol * 1e3 # CONVERT to VOLTS
    
    R = TPSD.variables['R_geo'][:,:] ## SHAPE (nTHETA, nRHO) I.E Theta is the 0th index, rho is the 1st.
    Z = TPSD.variables['Z_geo'][:,:]
    R_axis = TPSD.variables['R_axis'][:]
    Z_axis = TPSD.variables['Z_axis'][:]
    r = np.sqrt( (R-R_axis)**2 + (Z-Z_axis)**2 )
    
    if flux_averaged == True:
        ## IF we want the r_eff averaged over the flux surface (used in drho_pol/dr = 1/r_eff),
        #   then we can take r_avg, and find the gradient of the line of best fit which passes through (r=0, rho=0),
        #   which is equal to < r*rho >/<r^2> (see Simple Linear Regression)
        r_avg = np.mean(r[:,:], axis=0)
        Exp_r_avg_times_rho = np.mean(r_avg * rho_pol_edges_vs_rho_tor_edges)
        Exp_r_avg_squared = np.mean(r_avg**2)
        drho_pol_dr = Exp_r_avg_times_rho / Exp_r_avg_squared
        r_eff = 1/drho_pol_dr
        
    
    else:
        ## ELSE we will get the values for the OUTER MIDPLANE
        theta_zero_index = np.min(np.where(TPSD.variables['th_eq'][:] >= 0 ) )
        r_omp = r[theta_zero_index,:]
        ## NOTE THAT using the omp R-values with an EFIELD calculated using the FLUX SURAFCE AVERAGE grad(rho) is NOT consistent.
        # On the midplane, the rate of change of rho with respect to Z must be around 0; therefore, on the midplane, |grad(rho)| = |partial(rho)/partial(R)|
        # Therefore the same calculation with dPhi_electric_drho_tor with this new grad(rho) is needed.
        Exp_r_omp_times_rho = np.mean(r_omp * rho_pol_edges_vs_rho_tor_edges)
        Exp_r_omp_squared = np.mean(r_omp**2)
        drho_pol_dr = Exp_r_omp_times_rho / Exp_r_omp_squared
        r_eff = 1/drho_pol_dr
    
    nrho = len(rho_pol_edges_vs_rho_tor_edges)
    rhomin, rhomax = np.min(rho_pol_edges_vs_rho_tor_edges), np.max(rho_pol_edges_vs_rho_tor_edges)
    
    ### Finally, we must interpolate at the regular rho_pol grid we have defined, since our rho values are not evenly spaced.
    dPhi_electric_drho_pol_new = np.interp(x = np.linspace(rhomin, rhomax, nrho), xp = rho_pol_edges_vs_rho_tor_edges, fp = dPhi_electric_drho_pol)
    
    efield_dict = {'nrho': nrho, 'rhomin': rhomin, 'rhomax': rhomax, 'dvdrho': dPhi_electric_drho_pol_new, 'reff': r_eff}
    
    return efield_dict


def efield_from_full_CDF(TRANSP_full_CDF_dataset, time=0.65, plotyn=0):
    """ 
    No flux averaging in this case.
    
    """ 
    TFCD = TRANSP_full_CDF_dataset
    
    time_index = np.min(np.where(TFCD.variables['TIME'][:] > time))
    
    # Get variables from the netCDF file
    PLFLX = TFCD.variables['PLFLX'][time_index] # A function of XB
    rho_pol_edges_vs_rho_tor_edges = np.append(0, np.sqrt((PLFLX - 0)/(PLFLX[-1]  - 0))) ## NOTE that we are assuming psi on-axis to be ZERO....
    
    XB = TFCD.variables['XB'][time_index] # Zone edges
    XB = np.append(0, XB)  # Prepend the magnetic axis
    XILMP = TFCD.variables['XILMP'][time_index,:] # X-values corresponding to the R-values contained in RMAJM -- NOTE that these values are negative on the high field side, but only one side is needed.
    RMAJM = TFCD.variables['RMAJM'][time_index,:] * 1e-2 #CONVERT TO M      # R-values for data
    ## Now get the indices for the outer midplane (omp)
    RAXIS  = TFCD.variables['RAXIS'][time_index] * 1e-2 # CONVERT to metres
    omp_indices = np.where(np.round(RMAJM, 6) >= np.round(RAXIS, 6))
    
    rho_pol_vs_XILMP = np.interp(x=XILMP[omp_indices], xp=XB, fp=rho_pol_edges_vs_rho_tor_edges) ## Interpolate just in case the flux label values in XB are not all contained in XILMP.

    # ERTOT   NC radial E Field                 V/CM  RMAJM
    ERTOT  = TFCD.variables['ERTOT'][time_index][omp_indices] * 1e2 # CONVERT to volts per metre
    
    r = RMAJM[omp_indices] - RAXIS
    
    Exp_r_times_rho = np.mean(r * rho_pol_vs_XILMP)
    Exp_r_squared = np.mean(r**2)
    drho_pol_dr = Exp_r_times_rho / Exp_r_squared
    r_eff = 1/drho_pol_dr
    
    ## NOW we must multiply ERTOT ( equal to -dV/dr at the omp ) by -1* dr/drho_pol, because ASCOT5 requires dV/drho_pol as input
    dr_drho_pol = np.gradient(r, rho_pol_vs_XILMP, edge_order =2)
    dV_dRHO_POL = ERTOT * -1 * dr_drho_pol
    
    nrho = len(rho_pol_vs_XILMP)
    rhomin, rhomax = np.min(rho_pol_vs_XILMP), np.max(rho_pol_vs_XILMP)
    
    if plotyn == 1:
        plt.figure()
        plt.plot(rho_pol_vs_XILMP, ERTOT, label='CDF, ERTOT')
        plt.plot(rho_pol_vs_XILMP, -dV_dRHO_POL/r_eff, label='CDF, from dict vars')
        plt.xlabel('rho_pol')
        plt.ylabel('Radial E field [V/m]')
        plt.legend()
    
    ### Finally, we must interpolate at the regular rho_pol grid we have defined, since our rho values are not evenly spaced.
    dV_dRHO_POL_new = np.interp(x = np.linspace(rhomin, rhomax, nrho), xp = rho_pol_vs_XILMP, fp = dV_dRHO_POL)

    efield_dict = {'nrho': nrho, 'rhomin': rhomin, 'rhomax': rhomax, 'dvdrho': dV_dRHO_POL_new, 'reff': r_eff}
    
    return efield_dict




if __name__ == '__main__':
    plotyn = 1
    
    run_ID = '134020D30'
    TRANSP_plasma_state_directory = '../../134020D30/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    TRANSP_full_CDF_filename = TRANSP_plasma_state_directory + run_ID + ".CDF"
    TRANSP_full_CDF_data = nc.Dataset(TRANSP_full_CDF_filename)
    
    efield_test_dict_2 = efield_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, flux_averaged=False, plotyn=plotyn)
    efield_test_dict_CDF = efield_from_full_CDF(TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, time=0.65, plotyn=plotyn)
    
    rho_PS = np.linspace(efield_test_dict_2['rhomin'], efield_test_dict_2['rhomax'], efield_test_dict_2['nrho'])
    rho_CDF = np.linspace(efield_test_dict_CDF['rhomin'], efield_test_dict_CDF['rhomax'], efield_test_dict_CDF['nrho'])
    
    
    if plotyn == 1:
        # plt.figure()
        # plt.plot(efield_test_dict_2['dvdrho'])
        # plt.plot(efield_test_dict_CDF['dvdrho'])
        
        plt.figure()
        plt.plot(rho_PS, -1*efield_test_dict_2['dvdrho']/efield_test_dict_2['reff'], label='PS',  marker='x', markersize=13, markeredgewidth=2, fillstyle='none')
        plt.plot(rho_CDF, -1*efield_test_dict_CDF['dvdrho']/efield_test_dict_CDF['reff'], label='CDF', marker='+', markersize=13, markeredgewidth=1, fillstyle='none')
        plt.xlabel('rho_pol')
        plt.ylabel('Radial E field [V/m]')
        plt.legend()
        
        
        plt.show()

