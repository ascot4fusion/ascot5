max# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 18:16:12 2023

@author: Ian Dolby

WIP: THE exponential extrapolation method currently assumes that rho_pol(R) is roughly linear at the outer midplane. 
        This is so that less information has to be passed from the plasma and neutrals scripts, but it is of course possible to map them properly.
        IF I try to require the rho_pol <-> R mapping, this will produce problems for the FULL_CDF functions, which do not provide this information.

EXTRAPOLATION methods:
'Dummy' extrapolation means that we append small non-zero values from just beyond the LCFS out to the edge of the grid.
'Constant' means that the last value will be copied out to the edge of the grid.
'Linear' means that the slope of the last elements will be used, with a warning if the the slope is positive.
'Exponential' means that a SOL model will be used to calculate a decay constant for an exponential decay as a function of R_midplane;
       the SOL width is found in R.J. Goldston 2012 Nucl. Fusion 52 013009 .
    

This means that the extrapolation function needs:   1) The variable array to be extended, 
                                                    2) The rho array which must also be extended,
                                                    3) The rho value to be appended,
                                                    4) The R-values at the (outer) midplane IF the method is 'exponential', since the exponential
                                                        is in terms of R (although rho_poloidal is basically a linear function of R at the omp...),
                                                    5) The poloidal magnetic field at the LCFS for the decay constant IF the method is 'exponential',
                                                    6) The bulk ion temperature at the LCFS for the decay constant IF the method is 'exponential'.
"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import copy
import scipy.interpolate

from scipy.constants import physical_constants as const


def extrapolate_profile(variable_array, rho_array, rho_to_append=2, rho_nsteps=10, dummy_value=5, extrapolation_type='dummy', R_LCFS=None, Bpol_LCFS=None, T_LCFS=None, plotyn=0, **kwargs):
    """ SHOULD this function automatically detect 2D variable arrays and extrapolate intelligently, or is it clearer to just extrapolate e.g the ion density
         individually? (The ion density has shape (nrho, nion)).
    
    Below we ensure that the extrapolation_type string is not case-sensitive.
    
    The keyword args are intended for the variables needed by the SOL model in the exponential extrapolation section.
    """
    f = 1e-4
    
    if extrapolation_type.lower() == 'dummy':
        rho_extension_min = rho_array[-1]*(1 + f) ## We want to set the variable to the dummy value as soon as possible, unless we also want some linear extrapolation.
        rho_extension = np.linspace(rho_extension_min, rho_to_append, rho_nsteps)
        variable_extension = dummy_value * np.ones_like(rho_extension)
        
        rho_array_new = np.append(rho_array, rho_extension)
        variable_array_new = np.append(variable_array, variable_extension)
    
    elif extrapolation_type.lower() == 'constant':
        rho_extension_min = rho_array[-1]*(1 + f) ## We want to set the variable to the new value as soon as possible, unless we also want some linear extrapolation.
        rho_extension = np.linspace(rho_extension_min, rho_to_append, rho_nsteps)
        variable_extension = variable_array[-1] * np.ones_like(rho_extension)
        
        rho_array_new = np.append(rho_array, rho_extension)
        variable_array_new = np.append(variable_array, variable_extension)
        
    elif extrapolation_type.lower() == 'linear':
        rho_extension_min = rho_array[-1]*(1 + f) ## We want to set the variable to the new value as soon as possible, unless we also want some linear extrapolation.
        rho_extension = np.linspace(rho_extension_min, rho_to_append, rho_nsteps)
        
        variable_extension_slope = np.gradient(variable_array, rho_array)[-1]
        
        
        if variable_extension_slope >= 0:
            ## IF for some reason the slope is positive, get the slope for the last three elements instead. IF this is still positive, return an error.
            variable_extension_slope = (variable_array[-1] - variable_array[-3])/(rho_array[-1] - rho_array[-3])
            if variable_extension_slope > 0:
                print("** ERROR** Gradient is positive.")
                raise ValueError
        
        if variable_extension_slope == 0:
            print("NOTE: the gradient is zero. Extrapolating anyway.")
        
        variable_c =  variable_array[-1] - rho_array[-1] * variable_extension_slope
        variable_extension = rho_extension * variable_extension_slope  +  variable_c
        
        ## THEN we need to find where this becomes negative and replace the subsequent values with the dummy value.
        last_positive_idx = np.max(np.where(variable_extension >= 0))
        last_positive_val = variable_extension[last_positive_idx]
        variable_extension = np.where(variable_extension < 0, 0, variable_extension)
        
        
        rho_array_new = np.append(rho_array, rho_extension)
        variable_array_new = np.append(variable_array, variable_extension)
        
    
    elif extrapolation_type.lower() == 'exponential':
        ## FIRST check that the R, B and T values are present.
        if R_LCFS==None or Bpol_LCFS==None or T_LCFS==None:
            print("Please provide the R-values at the outer midplane, the poloidal B-field, and the bulk ion temperature at the LCFS")
            raise ValueError
        
        rho_extension_min = rho_array[-1]*(1 + f) ## We want to set the variable to the new value as soon as possible, unless we also want some linear extrapolation.
        rho_extension = np.linspace(rho_extension_min, rho_to_append, rho_nsteps)
        
        # kwargs = ion_mass=2, charge=1, major_radius=0.85, minor_radius=0.65
        ion_mass_amu = kwargs.pop("ion_mass", 2) 
        charge = kwargs.pop("charge", 1)
        R0 = kwargs.pop("major_radius", 0.85)
        a = kwargs.pop("minor_radius", 0.65)
        R_axis = kwargs.pop("R_axis", 1)
        
        decay_constant = SOL_width(T_separatrix=T_LCFS, B_P_separatrix=Bpol_LCFS, major_radius=R0, minor_radius=a, ion_mass_amu=ion_mass_amu, charge=charge)
        
        delta_R_delta_rho_pol = (R_LCFS - R_axis)/(rho_array[-1] - rho_array[0]) ## Get the rough gradient, treating rho_pol(R) as linear.
        R_extension = rho_extension * delta_R_delta_rho_pol + R_LCFS
        
        ## Now we fit an exponential for the provided quantity (Q): Q = A*exp[ -lambda_SOL * (R - R_LCFS) ] + C,
        ## for which A = T_LCFS - T_infinity, and C = T_infinity IF we want the decay to asymptotically approach zero.
        ## IF, however, we want it to approach a specific value at R_max, we have A = (T(R_max) - T_LCFS)/(exp_max - 1) and C = (T_LCFS * exp_max - T(R_max))/(exp_max - 1)
        T_Rmax = dummy_value
        exp_max = np.exp( -decay_constant*(R_extension[-1] - R_LCFS) )
        A = (T_Rmax - T_LCFS)/(exp_max - 1)
        C = (T_LCFS*exp_max - T_Rmax)/(exp_max - 1)
        
        
        variable_extension = A * np.exp( -decay_constant * (R_extension - R_LCFS) ) + C
        
        rho_array_new = np.append(rho_array, rho_extension)
        variable_array_new = np.append(variable_array, variable_extension)
    
    else:
        print("** ERROR**: Unrecognised extrapolation type.")
    
    extended_vars = {'rho': rho_array_new, 'quantity': variable_array_new}
    return extended_vars


def SOL_width(T_separatrix, B_P_separatrix, major_radius=0.85, minor_radius=0.65, ion_mass_amu=2, charge=1):
    ### The SOL width can be roughly determined from the bulk plasma. In the heuristic drift model [1], it is given by
    # lambda approx= 2a/R * r_L_POL, where a is the minor radius, R is the major radius,
    #  and r_L_POL is the poloidal ion gyroradius, for which v_perp and omega_c are evaluated
    #  using the POLOIDAL field B_P_separatrix [T].
    # v_perp is calculated from the ion temperature at the separatrix as v_perp = sqrt[ pi * T_i_sep / mass ],
    # which can be obtained by averaging v_perp over a Maxwellian velocity distribution.
    # This assumes an isotropic velocity distribution, at which point this integral is the
    # same for any chosen direction, so it is valid for any field direction.
    # [1]: R.J. Goldston 2012 Nucl. Fusion 52 013009
    aspect_ratio = major_radius/minor_radius
    omega_c = charge * const['elementary charge'][0] * B_P_separatrix / (ion_mass_amu * const['atomic mass constant'][0])
    v_perp = np.sqrt(np.pi * T_separatrix * const['elementary charge'][0] / (ion_mass_amu * const['atomic mass constant'][0]))
    
    lambda_SOL = (2/aspect_ratio) * v_perp / omega_c
    
    return lambda_SOL



if __name__ == '__main__':
    plotyn = 1
    
    run_ID = '134020D30'
    TRANSP_plasma_state_directory = '../../134020D30/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    
    TPSD = TRANSP_plasma_state_data
    
    rho_tor_edges = TPSD.variables['rho'][:]
    psipol = TPSD.variables['psipol'][:]  # Get the poloidal flux values which correspond to the rho_toroidal values (the profiles are provided in terms of rho_toroidal)
    psi_interp = scipy.interpolate.RegularGridInterpolator(points=(TPSD.variables['R_grid'][:], TPSD.variables['Z_grid'][:]), values=TPSD.variables['PsiRZ'][:].T, method='cubic')
    psi0  = psi_interp((TPSD.variables['R_axis'][:], TPSD.variables['Z_axis'][:]))
    psi1  = psipol[-1]
    rho_tor_centres = rho_tor_edges[:-1] + np.diff(rho_tor_edges) / 2 # Shift all of the rho-values by half a bin
    # psipol is already defined
    psipol_centres_vs_rho_tor_centres = np.interp(x = rho_tor_centres, xp = rho_tor_edges, fp = psipol)
    rho_pol_centres_vs_rho_tor_centres = np.sqrt((psipol_centres_vs_rho_tor_centres - psi0)/(psi1 - psi0))
    
    
    TE = TPSD.variables['Ts'][0] * 1e3  #CONVERT TO eV from keV
    TE_bdy = TPSD.variables['Te_bdy'][:] * 1e3  #CONVERT TO eV from keV
    rho_pol_centres_vs_rho_tor_centres = np.append(0, rho_pol_centres_vs_rho_tor_centres) # Now add back the central bin to ensure we still have a value defined at rho=0
    rho_pol_centres_vs_rho_tor_centres = np.append(rho_pol_centres_vs_rho_tor_centres, 1) # Now add the LCFS value so that the boundary values for the profiles can be added.
    TE = np.append(TE[0], TE)
    TE = np.append(TE, TE_bdy)
    RHO = rho_pol_centres_vs_rho_tor_centres
    
    theta_zero_index = np.min(np.where(TPSD.variables['th_eq'][:] >= 0 ) )
    
    R_LCFS = TPSD.variables['R_geo'][theta_zero_index, -1]
    R_axis = TPSD.variables['R_axis'][:]
    Z_axis = TPSD.variables['Z_axis'][:]
    Bpol = np.sqrt( TPSD.variables['BRRZ'][:]**2 + TPSD.variables['BZRZ'][:]**2).T
    
    Bpol_interp = scipy.interpolate.RegularGridInterpolator(points=(TPSD.variables['R_grid'][:], TPSD.variables['Z_grid'][:]), values=Bpol, method='cubic')
    Bpol_LCFS = Bpol_interp((R_LCFS, Z_axis))
    
    T_LCFS = TPSD.variables['Ti_bdy'][:] * 1e3  #CONVERT TO eV from keV
    
    # print(np.shape(TE))
    
    extended_TE_vars = extrapolate_profile(variable_array=TE, rho_array=RHO, rho_to_append=2, rho_nsteps=10, dummy_value=1, extrapolation_type='exponential', R_LCFS=R_LCFS, Bpol_LCFS=Bpol_LCFS, T_LCFS=T_LCFS, R_axis=R_axis, plotyn=plotyn)
    
    if plotyn == 1:
        plt.figure()
        plt.plot(extended_TE_vars['rho'], extended_TE_vars['quantity'])
        plt.ylabel('Electron temperature [eV]')
        plt.xlabel('rho_pol')
        
        
        plt.show()
    
