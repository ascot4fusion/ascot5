# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 18:15:45 2023

@author: Ian Dolby

"""
import netCDF4 as nc
import matplotlib.pyplot as plt

def wall_from_plasma_state_dataset(TRANSP_plasma_state_dataset, plotyn=0):
    """
    Args:
        fn : str <br>
            Full path to the HDF5 file.
        nelements : int <br>
            Number of wall segments.
        r : array_like (n,1) <br>
            R coordinates of wall segment vertices [m].
        z : array_like (n,1) <br>
            z coordinates of wall segment vertices [m].
        desc : str, optional <br>
            Input description.

    """
    
    TPSD = TRANSP_plasma_state_dataset
    
    wall_R = TPSD.variables['rlim'][:] # Already in m
    wall_Z = TPSD.variables['zlim'][:]
    nelements = len(wall_R)
    
    if plotyn == 1:
        plt.figure()
        plt.plot(wall_R, wall_Z, label='Limiter from plasma state file')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.legend()
    
    desc = 'Wall from plasma state file'
    
    plasma_dict = {'nelements': nelements, 'r': wall_R, 'z': wall_Z, 'desc':desc}  
    return plasma_dict

if __name__ == '__main__':
    
    plotyn=1
    
    run_ID = '134020D22'
    TRANSP_plasma_state_directory = '../../134020D22/'
    
    TRANSP_plasma_state_filename = TRANSP_plasma_state_directory + run_ID + "_ps_ts1_state.cdf"
    TRANSP_plasma_state_data = nc.Dataset(TRANSP_plasma_state_filename)
    
    ps_wall = wall_from_plasma_state_dataset(TRANSP_plasma_state_dataset = TRANSP_plasma_state_data, plotyn=plotyn)