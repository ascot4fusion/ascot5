# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:20:57 2023

@author: Ian Dolby

This script checks that the bfield data in the ASCOT5 input is identical to that in the input files, 
 and provides sanity checks that A) the rho_pol <--> R mapping matches that in the full CDF, and B) Other Bfield quantities match.

A5 = ASCOT5 hdf5 file, PS = plasma state file, CDF = full CDF file.

The comparisons for checking are:
    1) 1D rho_pol(R), A5 vs PS vs CDF
    2) 2D bfield & psi profiles, A5 vs PS
    3) 1D <Bpol> profiles, A vs, PS vs CDF
    4) 1D Psi profile, A vs, PS vs CDF


NOTE: Since this script was designed, the full Psi(R,Z) grid is available in the TRANSP output, but as older TRANSP files will not contain it, I have not used it.
"""
import numpy as np
import scipy as sc
import scipy.interpolate
import netCDF4 as nc
import matplotlib.pyplot as plt
import inputhelpers.get_Bfield as get_Bfield
import skimage

from a5py import Ascot

from scipy.constants import physical_constants as const


def get_ASCOT5_bfield_vars(ASCOT5_dataset, plotyn=0):
    """ This function prepares the variables which will be used in the comparisons.
    
    The required variables are:
        A) rho_pol(R_midplane) [1D]
        B) psi(R,Z) [2D]
        C) bphi(R,Z) [2D]
        D) br(R,Z) [2D]
        E) bz(R,Z) [2D]
        F) rho(R,Z) [2D]
        G) Bpol(R) flux-surface-averaged [1D] NOTE that thet flux-surface averaging must be done after 
            getting the flux surface boundaries from the plasma state file, so we do not calculate it here.
        H) psi(R_midplane) [1D]
    """
    # bfield = ASCOT5_dataset.bfield.active.read()
    bfield = ASCOT5_dataset.data.bfield.active.read()
    
    #Then we unload the magnetic field; the unit of length is METRES, and psi is the stream function = the poloidal flux divided by 2pi.
    axisr = bfield['axisr']
    axisz = bfield['axisz']
    bphi  = bfield['bphi']
    br    = bfield['br']
    bz    = bfield['bz']
    psi  = bfield['psi']
    psi0 = bfield['psi0'][0]
    psi1 = bfield['psi1'][0]
    Rmin = bfield['rmin'][0]
    Rmax = bfield['rmax'][0]
    nR   = bfield['nr'][0]
    Zmin = bfield['zmin'][0]
    Zmax = bfield['zmax'][0]
    nZ   = bfield['nz'][0]
    
    ## Now generate the B-field R & Z components from psi
    bfield_r = np.linspace(start=Rmin, stop=Rmax, num=nR)
    bfield_z = np.linspace(start=Zmin, stop=Zmax, num=nZ)
    r_diff = np.diff(bfield_r)[0]
    z_diff = np.diff(bfield_z)[0]
    psi_prime = np.gradient(psi, r_diff, z_diff, edge_order=2)

    # Calculate B_R, B_Z and B_phi from R and psi.
    r_mg, z_mg = np.meshgrid(bfield_r, bfield_z, indexing='ij')
    B_R = ((-1/r_mg)*psi_prime[1])
    B_Z = ((1/r_mg)*psi_prime[0])
    Bpol = np.sqrt(B_R**2 + B_Z**2)
    
    ## Now generate the midplane rho<->R map by interpolating the rho at the R-Z midplane values.
    R_array_2 = r_mg.flatten() # This gives us the R and Z values corresponding to psi in the form of a flattened meshgrid.
    Z_array_2 = z_mg.flatten()
    grid_R_Z_pairs = (R_array_2, Z_array_2)
    
    midplane_Z_array = axisz * np.ones_like(bfield_r)
    midplane_R_array = bfield_r
    midplane_R_Z_pairs = (midplane_R_array, midplane_Z_array) # These are the NEW positions, so midplane_Z_array should be the same size as bfield_r
    
    rho = np.sqrt( (psi - psi0) / (psi1 - psi0) )
    rho_flat = rho.flatten()
    psi_flat = psi.flatten()
    rho_midplane = scipy.interpolate.griddata(grid_R_Z_pairs, rho_flat, midplane_R_Z_pairs, method='cubic')
    psi_midplane = scipy.interpolate.griddata(grid_R_Z_pairs, psi_flat, midplane_R_Z_pairs, method='cubic')
    
    if plotyn == 1:
        plt.figure()
        plt.plot(midplane_R_array, rho_midplane, label='rho vs R_midplane from ASCOT5 file')
        plt.xlabel('R [m]')
        plt.ylabel('rho')
        plt.legend()


    bfield_dict = {'R_array': bfield_r, 'Z_array': bfield_z, 'R_midplane': midplane_R_array, 'rho_pol_midplane': rho_midplane, 'psi': psi, \
                   'bphi': bphi, 'br': B_R, 'bz': B_Z, 'rho_pol': rho, 'psi_midplane': psi_midplane, 'Bpol': Bpol, 'br_in_file': br, \
                   'bz_in_file': bz}
    return bfield_dict


def get_plasma_state_bfield_vars(TRANSP_plasma_state_dataset, COCOS_in=5, COCOS_out=3, plotyn=0):
    """ This function prepares the variables which will be used in the comparisons.
    
    The required variables are:
        A) rho_pol(R_midplane) [1D]
        B) psi(R,Z) [2D]
        C) bphi(R,Z) [2D]
        D) br(R,Z) [2D]
        E) bz(R,Z) [2D]
        F) rho(R,Z) [2D]
        G) Bpol(R) flux-surface-averaged [1D]
        H) psi(R_midplane) [1D]
    
    The flux-surface boundaries will also be returned as output, as rho_pol values evaluated at the rho_tor bin edges.
    """
    TPSD = TRANSP_plasma_state_dataset
    bfield = get_Bfield.bfield_from_plasma_state_dataset(TRANSP_plasma_state_dataset=TPSD, COCOS_in=COCOS_in, COCOS_out=COCOS_out, quiet=True)
    
    #Then we unload the magnetic field; the unit of length is METRES, and psi is the stream function = the poloidal flux divided by 2pi.
    axisr = bfield['axisr']
    axisz = bfield['axisz']
    bphi  = bfield['bphi']
    br    = bfield['br']
    bz    = bfield['bz']
    psi  = bfield['psi']
    psi0 = bfield['psi0']
    psi1 = bfield['psi1']
    Rmin = bfield['rmin']
    Rmax = bfield['rmax']
    nR   = bfield['nr']
    Zmin = bfield['zmin']
    Zmax = bfield['zmax']
    nZ   = bfield['nz']
    
    ## Now generate the B-field R & Z components from psi
    bfield_r = np.linspace(start=Rmin, stop=Rmax, num=nR)
    bfield_z = np.linspace(start=Zmin, stop=Zmax, num=nZ)
    r_diff = np.diff(bfield_r)[0]
    z_diff = np.diff(bfield_z)[0]
    psi_prime = np.gradient(psi, r_diff, z_diff, edge_order=2)

    # Calculate B_R, B_Z and B_phi from R and psi.
    r_mg, z_mg = np.meshgrid(bfield_r, bfield_z, indexing='ij')
    B_R = ((-1/r_mg)*psi_prime[1])
    B_Z = ((1/r_mg)*psi_prime[0])
    
    ## Now generate the midplane rho<->R map by interpolating the rho at the R-Z midplane values.
    R_array_2 = r_mg.flatten() # This gives us the R and Z values corresponding to psi in the form of a flattened meshgrid.
    Z_array_2 = z_mg.flatten()
    grid_R_Z_pairs = (R_array_2, Z_array_2)
    
    midplane_Z_array = axisz * np.ones_like(bfield_r)
    midplane_R_array = bfield_r
    midplane_R_Z_pairs = (midplane_R_array, midplane_Z_array) # These are the NEW positions, so midplane_Z_array should be the same size as bfield_r
    
    rho = np.sqrt( (psi - psi0) / (psi1 - psi0) )
    rho_flat = rho.flatten()
    psi_flat = psi.flatten()
    rho_midplane = scipy.interpolate.griddata(grid_R_Z_pairs, rho_flat, midplane_R_Z_pairs, method='cubic')
    psi_midplane = scipy.interpolate.griddata(grid_R_Z_pairs, psi_flat, midplane_R_Z_pairs, method='cubic')
    
    ### NOW we get variables from the plasma state file for Bpol comparisons.
    psipol = TPSD.variables['psipol'][:]  ## REMEMBER that this has not been adjusted for COCOS
    sigma_Bp_eff = np.sign(psipol[-1]) * np.sign(psi1)
    rho_pol_at_rho_tor_edges = np.sqrt((psipol*sigma_Bp_eff - psi0)/(psi1 - psi0))
    
    ## FOR our average over flux surfaces, we use the R and Z values which have already been mapped to rho_tor_edges and theta_edges.
    ## Although we have B_R and B_Z calculated above, we will get them directly from the plasma state file again for the interpolator, to minimise error.
    R_flux_surfaces = TPSD.variables['R_geo'][:]
    Z_flux_surfaces = TPSD.variables['Z_geo'][:]
    ## Remember that the values in the plasma state file have shape (nZ, nR)
    flux_surface_R_Z_pairs = (R_flux_surfaces.flatten(), Z_flux_surfaces.flatten())
    BpolRZ = np.sqrt(TPSD.variables['BRRZ'][:]**2 + TPSD.variables['BZRZ'][:]**2).T
    Bpol_geo_flat = scipy.interpolate.griddata(grid_R_Z_pairs, BpolRZ.flatten(), flux_surface_R_Z_pairs, method='cubic')
    Bpol_geo = np.reshape(Bpol_geo_flat, np.shape(R_flux_surfaces)) ## Reshape into the form (nTheta, nrho_tor)
    
    Bpol_avg = np.mean(Bpol_geo, axis=0) ## Average over all theta values. This should have shape nrho_tor.
    
    
    if plotyn == 1:
        plt.figure()
        plt.plot(midplane_R_array, rho_midplane, label='rho vs R_midplane from plasma state file')
        plt.xlabel('R [m]')
        plt.ylabel('rho')
        plt.legend()
        
    
    bfield_dict = {'R_array': bfield_r, 'Z_array': bfield_z, 'R_midplane': midplane_R_array, 'rho_pol_midplane': rho_midplane, 'psi': psi, \
                   'bphi': bphi, 'br': B_R, 'bz': B_Z, 'rho_pol': rho, 'psi_midplane': psi_midplane, 'Bpol': BpolRZ, 'Bpol_avg': Bpol_avg, \
                   'rho_pol_at_rho_tor_edges': rho_pol_at_rho_tor_edges}
    return bfield_dict


def get_full_CDF_bfield_components(TRANSP_full_CDF_dataset, time=0.65, plotyn=0):
    """ This function prepares the variables which will be used in the comparisons.
    
    The required variables are:
        A) rho_pol(R_midplane) [1D]
        B) Bpol(R) flux-surface-averaged [1D]
        C) psi(R_midplane) [1D]
    
    The flux-surface boundaries will also be returned as output, as rho_pol values evaluated at the rho_tor bin edges.
    """
    time_index = np.min(np.where(TRANSP_full_CDF_dataset.variables['TIME'][:] > time))
    
    # Get variables from the netCDF file
    PLFLX = TRANSP_full_CDF_dataset.variables['PLFLX'][time_index] # A function of XB
    BPOL = TRANSP_full_CDF_dataset.variables['BPOL'][time_index] # A function of XB
    X  = TRANSP_full_CDF_dataset.variables['X'][time_index,:]  # Zone centres
    XB = TRANSP_full_CDF_dataset.variables['XB'][time_index]
    XB = np.append(0, XB)  # Prepend the magnetic axis
    XIRSYM = TRANSP_full_CDF_dataset.variables['XIRSYM'][time_index,:] # X-values corresponding to the R-values contained in RMJSYM -- NOTE that these values are negative on the high field side, but only one side is needed.
    RMJSYM = TRANSP_full_CDF_dataset.variables['RMJSYM'][time_index,:] * 1e-2 #CONVERT TO M      # R-values for data
    
    X_in_XIRSYM_mask = [XIRSYMval in X for XIRSYMval in np.abs(XIRSYM) ] # XIRSYM but with the value True when it is also in X, and False if not.
    RMJSYM_at_X_vals = RMJSYM[X_in_XIRSYM_mask]
    
    XB_in_XIRSYM_mask = [XIRSYMval in XB for XIRSYMval in np.abs(XIRSYM) ] # XIRSYM but with the value True when it is also in XB, and False if not.
    RMJSYM_at_XB_vals = RMJSYM[XB_in_XIRSYM_mask]
    
    R_array_edges = RMJSYM_at_XB_vals
    R_array_centres = RMJSYM_at_X_vals
    
    
    RHO_pol_vs_rho_tor_edges = np.append(0, np.sqrt((PLFLX - 0)/(np.max(PLFLX)  - 0)))
    RHO_pol_vs_rho_tor_edges_full = np.append(np.flip(RHO_pol_vs_rho_tor_edges[1:]), RHO_pol_vs_rho_tor_edges)
    
    RHO_pol_vs_rho_tor_centres      = np.interp(x=X, xp=XB, fp=RHO_pol_vs_rho_tor_edges)
    RHO_pol_vs_rho_tor_centres_full = np.interp(x=RMJSYM_at_X_vals, xp=RMJSYM_at_XB_vals, fp=RHO_pol_vs_rho_tor_edges_full)
    
    XB_for_plot = np.abs(XIRSYM[XB_in_XIRSYM_mask])
    X_for_plot = np.abs(XIRSYM[X_in_XIRSYM_mask])
    
    RAXIS = TRANSP_full_CDF_dataset.variables['RAXIS'][time_index]  * 1e-2 #CONVERT TO M
    ZAXIS = TRANSP_full_CDF_dataset.variables['YAXIS'][time_index]  * 1e-2 #CONVERT TO M
    
    RAXIS_index = np.max(np.where(RMJSYM_at_XB_vals <= RAXIS))
    
    """
    Now we load the X values so that RHO_poloidal_full_CDF can be adjusted to allow plotting of the plasma quantities, which are provided in terms of X.
    """
    
    # Prepend the axis value:
    BPOL = np.append(0, BPOL)
    PLFLX = np.append(0, PLFLX)
    R_array_edges = RMJSYM_at_XB_vals
    R_array_edges_omp = R_array_edges[RAXIS_index:]
    
    bfield_dict = {'time_index': time_index,   'XB': XB,                                 'XIRSYM': XIRSYM,     'X': X, 'RAXIS': RAXIS,     'ZAXIS': ZAXIS, \
                           'RMJSYM': RMJSYM,   'XB_in_XIRSYM_mask': XB_in_XIRSYM_mask,   'RMJSYM_at_XB_vals': RMJSYM_at_XB_vals,  \
                           'BPOL': BPOL,  'PLFLX': PLFLX,  'RHO_pol_vs_rho_tor_edges': RHO_pol_vs_rho_tor_edges,  'RHO_pol_vs_rho_tor_edges_full': RHO_pol_vs_rho_tor_edges_full, \
                           'X_in_XIRSYM_mask': X_in_XIRSYM_mask, 'RMJSYM_at_X_vals': RMJSYM_at_X_vals, 'XB_for_plot': XB_for_plot, 'X_for_plot': X_for_plot, \
                           'RHO_pol_vs_rho_tor_centres': RHO_pol_vs_rho_tor_centres,  'RHO_pol_vs_rho_tor_centres_full': RHO_pol_vs_rho_tor_centres_full, \
                           'R_array_edges': R_array_edges, 'R_array_centres': R_array_centres, 'R_array_edges_omp': R_array_edges_omp}
    return bfield_dict


def compare_bfield(ASCOT5_dataset, TRANSP_plasma_state_dataset, TRANSP_full_CDF_dataset, COCOS_in=5, COCOS_out=3, time=None, plotyn=0):
    """ In this function, we make the following comparisons:
            1) 1D rho_pol(R), A5 vs PS vs CDF
            2) 2D bfield & psi profiles, A5 vs PS
            3) 1D <Bpol> profiles, A vs, PS vs CDF
            4) 1D Psi profile, A vs, PS vs CDF
        
        Comparing <Bpol> requires a little bit of work because we need to average over rho_pol contours, but 
         the values of rho_pol to use are only provided by the other two datasets.
    A5: bfield_dict = {'R_array': bfield_r, 'Z_array': bfield_z, 'R_midplane': midplane_R_array, 'rho_pol_midplane': rho_midplane, 'psi': psi, \
                       'bphi': bphi, 'br': B_R, 'bz': B_Z, 'rho_pol': rho, 'psi_midplane': psi_midplane, 'Bpol': Bpol, , 'br_in_file': br, \
                       'bz_in_file': bz}
    
    PS: bfield_dict = {'R_array': bfield_r, 'Z_array': bfield_z, 'R_midplane': midplane_R_array, 'rho_pol_midplane': rho_midplane, 'psi': psi, \
                   'bphi': bphi, 'br': B_R, 'bz': B_Z, 'rho_pol': rho, 'psi_midplane': psi_midplane, 'Bpol': BpolRZ, 'Bpol_avg': Bpol_avg, \
                   'rho_pol_at_rho_tor_edges': rho_pol_at_rho_tor_edges, 'Bpol_2': Bpol_2, 'contour_points': contour_points, 'Bpol_avg_2': Bpol_avg_2}
    
    CDF: bfield_dict = {'time_index': time_index,   'XB': XB,                                 'XIRSYM': XIRSYM,     'X': X, 'RAXIS': RAXIS,     'ZAXIS': ZAXIS, \
                           'RMJSYM': RMJSYM,   'XB_in_XIRSYM_mask': XB_in_XIRSYM_mask,   'RMJSYM_at_XB_vals': RMJSYM_at_XB_vals,  \
                           'BPOL': BPOL,  'PLFLX': PLFLX,  'RHO_pol_vs_rho_tor_edges': RHO_pol_vs_rho_tor_edges,  'RHO_pol_vs_rho_tor_edges_full': RHO_pol_vs_rho_tor_edges_full, \
                           'X_in_XIRSYM_mask': X_in_XIRSYM_mask, 'RMJSYM_at_X_vals': RMJSYM_at_X_vals, 'XB_for_plot': XB_for_plot, 'X_for_plot': X_for_plot, \
                           'RHO_pol_vs_rho_tor_centres': RHO_pol_vs_rho_tor_centres,  'RHO_pol_vs_rho_tor_centres_full': RHO_pol_vs_rho_tor_centres_full, \
                           'R_array_edges': R_array_edges, 'R_array_centres': R_array_centres}
    
    """
    if time == None:
        time = TRANSP_plasma_state_dataset.variables['t1'][:]
    
    ### CALL the functions for preparing the Bfield data
    A5 = get_ASCOT5_bfield_vars(ASCOT5_dataset=ASCOT5_dataset, plotyn=0)
    PS = get_plasma_state_bfield_vars(TRANSP_plasma_state_dataset=TRANSP_plasma_state_dataset, COCOS_in=COCOS_in, COCOS_out=COCOS_out, plotyn=0)
    CDF = get_full_CDF_bfield_components(TRANSP_full_CDF_dataset=TRANSP_full_CDF_dataset, time=time, plotyn=0)
    
    ### FIRST, we interpolate the 2D profiles onto the ASCOT5 file's grid
    B_R_PS_interp = scipy.interpolate.RegularGridInterpolator(points=(PS['R_array'], PS['Z_array']), values=PS['br'], bounds_error=False, fill_value=np.nan)
    B_phi_PS_interp = scipy.interpolate.RegularGridInterpolator(points=(PS['R_array'], PS['Z_array']), values=PS['bphi'], bounds_error=False, fill_value=np.nan)
    B_Z_PS_interp = scipy.interpolate.RegularGridInterpolator(points=(PS['R_array'], PS['Z_array']), values=PS['bz'], bounds_error=False, fill_value=np.nan)
    psi_PS_interp = scipy.interpolate.RegularGridInterpolator(points=(PS['R_array'], PS['Z_array']), values=PS['psi'], bounds_error=False, fill_value=np.nan)
    rho_pol_PS_interp = scipy.interpolate.RegularGridInterpolator(points=(PS['R_array'], PS['Z_array']), values=PS['rho_pol'], bounds_error=False, fill_value=np.nan)
    
    A5_R_mg, A5_Z_mg = np.meshgrid(A5['R_array'], A5['Z_array'], indexing='ij')
    A5_R_Z_pairs = (A5_R_mg.flatten(), A5_Z_mg.flatten())
    
    B_R_PS_at_A5 = np.reshape(B_R_PS_interp(A5_R_Z_pairs), np.shape(A5_R_mg))
    B_phi_PS_at_A5 = np.reshape(B_phi_PS_interp(A5_R_Z_pairs), np.shape(A5_R_mg))
    B_Z_PS_at_A5 = np.reshape(B_Z_PS_interp(A5_R_Z_pairs), np.shape(A5_R_mg))
    psi_PS_at_A5 = np.reshape(psi_PS_interp(A5_R_Z_pairs), np.shape(A5_R_mg))
    rho_pol_PS_at_A5 = np.reshape(rho_pol_PS_interp(A5_R_Z_pairs), np.shape(A5_R_mg))
    

    ### THEN WE PREPARE Bpol_avg, starting by getting the rho contours along with to average
    PS_rho_pol = PS['rho_pol_at_rho_tor_edges']
    CDF_rho_pol = CDF['RHO_pol_vs_rho_tor_edges']
    
    contour_points_PS = []
    contour_points_CDF = []
    A5_Bpol_avg_PS = np.zeros_like(PS_rho_pol)
    A5_Bpol_avg_CDF = np.zeros_like(CDF_rho_pol)
    A5_Bpol_interp = scipy.interpolate.RegularGridInterpolator(points=(A5['R_array'], A5['Z_array']), values=A5['Bpol'], bounds_error=False, fill_value=np.nan)
    A5_rho = A5['rho_pol']
    A5_R = A5['R_array']
    A5_Z = A5['Z_array']
    nR = len(A5_R); nZ = len(A5_Z)
    
    ### START with PS rho_pol levels
    for i in range(1, len(A5_Bpol_avg_PS)):
        ## We ignore the rho=0 value, since there should not be any contour there.
        indices = np.array(skimage.measure.find_contours(A5_rho, PS_rho_pol[i])).squeeze() ## Z values first, R values second
        R_indices = indices.T[0]
        Z_indices = indices.T[1]
        
        R_points = np.interp(x=R_indices, xp = range(nR), fp=A5_R)
        Z_points = np.interp(x=Z_indices, xp = range(nZ), fp=A5_Z)
        points = np.array([R_points, Z_points])
        contour_points_PS.append(points) ## contour_points is a LIST
        Bpol_at_points = A5_Bpol_interp(tuple(points))
        average_bpol = np.mean(Bpol_at_points)
        A5_Bpol_avg_PS[i] = average_bpol
    
    ### THEN use the CDF rho_pol levels
    for i in range(1, len(A5_Bpol_avg_CDF)):
        ## We ignore the rho=0 value, since there should not be any contour there.
        indices = np.array(skimage.measure.find_contours(A5_rho, CDF_rho_pol[i])).squeeze() ## Z values first, R values second
        R_indices = indices.T[0]
        Z_indices = indices.T[1]
        
        R_points = np.interp(x=R_indices, xp = range(nR), fp=A5_R)
        Z_points = np.interp(x=Z_indices, xp = range(nZ), fp=A5_Z)
        points = np.array([R_points, Z_points])
        contour_points_CDF.append(points) ## contour_points is a LIST
        Bpol_at_points = A5_Bpol_interp(tuple(points))
        average_bpol = np.mean(Bpol_at_points)
        A5_Bpol_avg_CDF[i] = average_bpol
    
    
    ### NOW we PLOT
    if plotyn == 1:
        #(1) 1D rho_pol(R), A5 vs PS vs CDF
        plt.figure()
        plt.plot(A5['R_midplane'], A5['rho_pol_midplane'], label='ASCOT5 file')
        plt.plot(PS['R_midplane'], PS['rho_pol_midplane'], label='Plasma state file')
        plt.plot(CDF['R_array_edges'], CDF['RHO_pol_vs_rho_tor_edges_full'], label='CDF file')
        plt.xlabel('R [m]')
        plt.ylabel('rho_pol')
        plt.legend()
        
        #(2) 2D bfield & psi profiles, A5 vs PS
        # First, we get the extreme values for the colourbar
        
        plt.figure()
        plt.pcolor(A5['R_array'], A5['Z_array'], ((A5['br'] - B_R_PS_at_A5)/A5['br']).T, vmin=-1, vmax=1, cmap='bwr_r')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('BR, (ASCOT5 - Plasma state)/ASCOT5')
        plt.colorbar()
        
        plt.figure()
        plt.pcolor(A5['R_array'], A5['Z_array'], ((A5['bphi'] - B_phi_PS_at_A5)/A5['bphi']).T, vmin=-1, vmax=1, cmap='bwr_r')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('Bphi, (ASCOT5 - Plasma state)/ASCOT5')
        plt.colorbar()
        
        plt.figure()
        plt.pcolor(A5['R_array'], A5['Z_array'], ((A5['bz'] - B_Z_PS_at_A5)/A5['bz']).T, vmin=-1, vmax=1, cmap='bwr_r')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('BZ, (ASCOT5 - Plasma state)/ASCOT5')
        plt.colorbar()
        
        plt.figure()
        plt.pcolor(A5['R_array'], A5['Z_array'], ((A5['psi'] - psi_PS_at_A5)/A5['psi']).T, vmin=-1, vmax=1, cmap='bwr_r')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('psi, (ASCOT5 - Plasma state)/ASCOT5')
        plt.colorbar()
        
        plt.figure()
        plt.pcolor(A5['R_array'], A5['Z_array'], ((A5['rho_pol'] - rho_pol_PS_at_A5)/A5['rho_pol']).T, vmin=-1, vmax=1, cmap='bwr_r')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('rho_pol, (ASCOT5 - Plasma state)/ASCOT5')
        plt.colorbar()
        
        #(3) 1D <Bpol>(R), A5 vs PS vs CDF
        plt.figure()
        plt.plot(PS_rho_pol, A5_Bpol_avg_PS, label='ASCOT5 file, using PS rho_pol')
        plt.plot(CDF_rho_pol, A5_Bpol_avg_CDF, label='ASCOT5 file, using CDF rho_pol')
        plt.plot(PS_rho_pol, PS['Bpol_avg'], label='Plasma state file')
        plt.plot(CDF_rho_pol, CDF['BPOL'], label='CDF file')
        plt.xlabel('rho_pol')
        plt.ylabel('Bpol_avg [T]')
        plt.legend()
        
        #(4) 1D psi(R), A5 vs PS vs CDF
        plt.figure()
        plt.plot(A5['R_midplane'], np.abs(A5['psi_midplane']), label='ASCOT5 file')
        plt.plot(PS['R_midplane'], np.abs(PS['psi_midplane']), label='Plasma state file')
        plt.plot(CDF['R_array_edges_omp'], np.abs(CDF['PLFLX']), label='CDF file')
        plt.xlabel('R [m]')
        plt.ylabel('|psi| [Wb/rad]')
        plt.legend()
        
        
        plt.show()
    
    ## Return all of the variables we compare in case debugging is needed.
    comparison_variables_dict = {'R_midplane_A5': A5['R_midplane'],   'R_midplane_PS': PS['R_midplane'],     'R_midplane_CDF': CDF['R_array_edges_omp'], \
                                 'rho_pol_midplane_A5': A5['rho_pol_midplane'], 'rho_pol_midplane_PS': PS['rho_pol_midplane'], \
                                 'rho_pol_midplane_CDF': CDF['RHO_pol_vs_rho_tor_edges_full'], \
                                 'br_A5': A5['br'], 'br_PS': B_R_PS_at_A5, 'bphi_A5': A5['bphi'],  'bphi_PS': B_phi_PS_at_A5, \
                                 'bz_A5': A5['bz'], 'bz_PS': B_Z_PS_at_A5, 'psi_A5': A5['psi'],  'psi_PS': psi_PS_at_A5, \
                                 'rho_A5': A5['rho_pol'], 'rho_PS': rho_pol_PS_at_A5, 'rho_midplane_PS': PS_rho_pol,  'rho_midplane_CDF': CDF_rho_pol, \
                                 'Bpol_avg_A5_from_PS': A5_Bpol_avg_PS,   'Bpol_avg_A5_from_CDF': A5_Bpol_avg_CDF, 'Bpol_avg_PS': PS['Bpol_avg'], \
                                 'Bpol_avg_CDF': CDF['BPOL'], \
                                 'psi_midplane_A5': A5['R_midplane'], 'psi_midplane_PS': PS['psi_midplane'], 'psi_midplane_CDF': CDF['PLFLX']}
    return comparison_variables_dict


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
    
    # bfield_test_dict_A5 = get_ASCOT5_bfield_vars(ASCOT5_dataset=ASCOT5_data, plotyn=plotyn)
    # bfield_test_dict_PS = get_plasma_state_bfield_vars(TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, COCOS_in=5, COCOS_out=3, plotyn=plotyn)
    # bfield_test_dict_CDF = get_full_CDF_bfield_components(TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, time=0.65, plotyn=plotyn)
    
    comparison = compare_bfield(ASCOT5_dataset=ASCOT5_data, TRANSP_plasma_state_dataset=TRANSP_plasma_state_data, \
                                 TRANSP_full_CDF_dataset=TRANSP_full_CDF_data, plotyn=plotyn)
    
    
    if plotyn == 1:
        # plt.figure()
        # plt.pcolor(bfield_test_dict_PS['R_array'], bfield_test_dict_PS['Z_array'], bfield_test_dict_PS['Bpol'].T)
        # plt.xlabel('R [m]')
        # plt.ylabel('Z [m]')
        # plt.colorbar()
                
        # plt.figure()
        # plt.plot(bfield_test_dict_PS['rho_pol_at_rho_tor_edges'], bfield_test_dict_PS['Bpol_avg'], label='Bpol_avg vs rho_pol_edges from plasma state file')
        # plt.plot(bfield_test_dict_CDF['RHO_pol_vs_rho_tor_edges'], bfield_test_dict_CDF['BPOL'], label='Bpol_avg vs rho_pol_edges from full CDF file')
        # plt.xlabel('rho_pol')
        # plt.ylabel('Bpol_avg [T]')
        # plt.legend()
        
        
        plt.show()


