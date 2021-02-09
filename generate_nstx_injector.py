"""
Generate the NSTX neutral beam injector geometry for BBNBI5:
based on the AUG geometry file made by Juan Francisco Rivero


A. Sperduti <andrea.sperduti@physics.uu.se>
"""

import numpy as np
from scipy.constants import physical_constants as const
import a5py.marker.interpret as a5interpret
import a5py.ascot5io.nbi as nbi

def generate(fn='nstx_nbi.h5', inj=[1,2,3,4,5,6], specie = "D", desc="NSTX_NBI"):
    """
    Generate NSTX injector.

    Args:
        fn : str
            Full path to the HDF5 file.
        inj : list of int, default = all
            List of NSTX injectors from 1 to 6 to be generated.
        desc : str, optional
            Input description.
        specie : {"D","H"}
            Specie injected by the NBI. Deuterium "D" (default) 
         Rcenter, phicenter, zcenter: float
           Position of the center point of the injector grid [m], [deg], [m]
        thetasteer : float
            Angle around the horizontal axis by which the beam center line is oriented)
        phisteer : float
            Angle around the vertical axis by which the beam center line is oriented)
        tilt : float
            Angle around the horizontal axis by which half of the beamline grid is tilted 
            and looks like figure 1 in https://doi.org/10.1016/j.cpc.2014.10.024
        focal_length : float
            Distance to the focus point from the center of the grid
        div_h,div_v,div_halo_frac,div_halo_h and div_halo_v 
            Horizontal and vertical divergences, including the halo contribution
        tilt : float, optional
            Vertical tilt angle of the beam centerline [deg]
            
    """
    
    inj = np.array(inj)-1  
    Rcenter = np.array([11.3189, #1A tanrad = 0.69 m  
                        11.3573, #1B tanrad = 0.59 m
                        11.3887, #1C tanrad = 0.49 m
                        11.3164, #2A tanrad = 1.30 m
                        11.3932, #2B tanrad = 1.20 m
                        11.4595])#2C tanrad = 1.10 m
    
    phicenter = np.array([56.91, #1A
                          60.38, #1B
                          63.85, #1C
                          103.09, #2A
                          106.54, #2B
                          109.98])*np.pi/180 #2C
    
    zcenter = np.array([0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0])
    
    thetasteer = np.array([0.0,
                           0.0,
                           0.0,
                           0.0,
                           0.0,
                           0.0])
    
    phisteer = np.array([-3.4887,
                         -2.9740,
                         -2.4130,
                         -6.5538,
                         -6.0131,
                         -5.4835])*np.pi/180
    
    focal_length_v = np.array([0.988E3,
                               0.988E3,
                               0.988E3,
                               0.988E3,
                               0.988E3,
                               0.988E3])
    focal_length_h = np.array([0.988E3,
                               0.988E3,
                               0.988E3,
                               0.988E3,
                               0.988E3,
                               0.988E3])
    tilt = np.pi/180.
    div_h = np.array([4.94E-3,
                      4.94E-3,
                      4.94E-3,
                      4.94E-3,
                      4.94E-3,
                      4.94E-3])
    
    div_v = np.array([1.36E-2,
                      1.36E-2,
                      1.36E-2,
                      1.36E-2,
                      1.36E-2,
                      1.36E-2,
                      1.36E-2,
                      1.36E-2])
    div_halo_frac = np.array([0, #unknown data
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0])
    div_halo_h= np.array([1e-10, #unknown data
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10])
    div_halo_v= np.array([1e-10, #unknown data
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10,
                          1e-10])
    
    anum = 2
    znum = 1
    mass = a5interpret.mass_kg("deuterium")
    energy = np.array([1.0e5,
                       1.0e5,
                       1.0e5,
                       1.0e5,
                       1.0e5,
                       1.0e5])*const["elementary charge"][0]
    power = np.array([2.5,
                      2.5,
                      2.5,
                      2.5,
                      2.5,
                      2.5])*1e6
    efrac = np.array([[0.65,0.25,0.1],
                      [0.65,0.25,0.1],
                      [0.65,0.25,0.1],
                      [0.65,0.25,0.1],
                      [0.65,0.25,0.1],
                      [0.65,0.25,0.1]])
    
     # Make nbi input of the selected injectors
    out = []
    for i in inj:
     ## Define beamlets grid
            grid = define_grid(Rcenter[i],phicenter[i],zcenter[i],
                               thetasteer[i],phisteer[i],
                               tilt,focal_length_h[i],focal_length_v[i])

            out.append({})
            out[-1]["id"] = i+1
            out[-1]["nbeamlet"] = grid["xyz"][:,0].size
            out[-1]["beamletx"] = grid["xyz"][:,0]
            out[-1]["beamlety"] = grid["xyz"][:,1]
            out[-1]["beamletz"] = grid["xyz"][:,2]
            out[-1]["beamletdx"] = grid["dxdydz"][:,0]
            out[-1]["beamletdy"] = grid["dxdydz"][:,1]
            out[-1]["beamletdz"] = grid["dxdydz"][:,2]
            out[-1]["div_h"] = div_h[i]
            out[-1]["div_v"] = div_v[i]
            out[-1]["div_halo_frac"] = div_halo_frac[i]
            out[-1]["div_halo_h"] = div_halo_h[i]
            out[-1]["div_halo_v"] = div_halo_v[i]
            out[-1]["anum"] = anum
            out[-1]["znum"] = znum
            out[-1]["mass"] = mass
            out[-1]["energy"] = energy[i]
            out[-1]["efrac"] = efrac[i,:]
            out[-1]["power"] = power[i] 
    return nbi.write_hdf5(fn, out, desc)

def define_grid(Rcenter,phicenter,zcenter,
                thetasteer,phisteer,
                tilt,fl_h,fl_v):

    """
    Define beamlets grid in local coordinates  (centre of 
    grid at x=0, y=0, z=0) and then allocate and orientate 
    the grid according to the given location and orientation
    of the beam.

    Args:
        Rcenter, phicenter, zcenter : floats [m,rad,m].
            Location of the grid center in cilindrical coordinates.
        thetasteer: float [rad]
            Angle around the horizontal axis by which the beam center
            line is oriented.
        phisteer: float [rad]
            Angle around the vertical axis by which the beam center
            line is oriented.
        tilt: float [rad]
            Angle around the horizontal axis by which half of the beamline
            grid is tilted. This angle makes grid to be tilted in halves,
            producing a shape similar to the one shown below:
                \
                 \
                  \
                  /
                 /
                /
        fl_h, fl_v: floats [m,m]
            Horizontal and vertical focal lengths.
    """

    grid = {}

    ## Define beamlets position in local coordinates

    # Divide each gridhalf to 6 pieces: 2 of which create the main
    # diamond pattern, and 4 that add the four top rows.
    # Generate a rows-by-columns matrix of beamlets where
    # the extreme points are given as (front view):
    #
    # TopLeft ---------- LastCol
    #  \                       \
    #   \                       \
    #    \                       \
    #   LastRow ------------------   

    Rows = np.array([10,10,10])
    Cols = np.array([10,10,10])

    TopLeft = np.array([[-0.11, 0.3025],
                        [-0.11, 0.2825],
                        [-0.09, 0.3025]])

    LastRow = np.array([[-0.11, 0],
                        [-0.11, 0.02],
                        [-0.09, 0]])


    LastCol = np.array([[0.11, 0.3025],
                        [0.11, 0.2825],
                        [0.09, 0.3025]])

    RowStep = ((LastRow - TopLeft).T/np.where(Rows>1,Rows-1,1)).T
    ColumnStep = ((LastCol - TopLeft).T/np.where(Cols>1,Cols-1,1)).T

    ## Create grid in 2D

    n = 0
    nbeamlets = np.sum(Rows*Cols)
    grid["xyz"] = np.zeros((nbeamlets,2))
    for k in range(Rows.size):
        for i in range(Rows[k]):
            for j in range(Cols[k]):
                grid["xyz"][i+j*Rows[k]+n,:] = TopLeft[k] + i*RowStep[k] + j*ColumnStep[k]
        n = n + Rows[k]*Cols[k]

    ## Add a third dimension to the grid

    grid["xyz"] = np.array([np.zeros(nbeamlets),grid["xyz"][:,0],grid["xyz"][:,1]]).T
    
    ## Tilt the gridhalf

    A_tilt = np.array([[np.cos(tilt), 0, -np.sin(tilt)],
                       [0,1,0],
                       [np.sin(tilt), 0, np.cos(tilt)]])

    grid["xyz"] = np.matmul(A_tilt,grid["xyz"].T).T

    ## Second gridhalf as a mirror of the first

    grid["xyz"] = np.append(grid["xyz"],
                           np.array([grid["xyz"][:,0],grid["xyz"][:,1],-grid["xyz"][:,2]]).T,
                           axis=0)
    nbeamlets = 2*nbeamlets

    ## Define the beamlets direction in local coordinates

    grid["dxdydz"] = np.array([-np.ones(nbeamlets),
                               -grid["xyz"][:,1]/(fl_h - grid["xyz"][:,0]),
                               -grid["xyz"][:,2]/(fl_v - grid["xyz"][:,0])])

    grid["dxdydz"] = (grid["dxdydz"]/np.linalg.norm(grid["dxdydz"],axis=0)).T

    ## Reorientate and allocate the grid to its specific orientation and location

    A_thetasteer = np.array([[np.cos(thetasteer), 0, -np.sin(thetasteer)],
                       [0,1,0],
                       [np.sin(thetasteer), 0, np.cos(thetasteer)]])

    A_phisteer = np.array([[np.cos(phisteer),-np.sin(phisteer),0],
                           [np.sin(phisteer),np.cos(phisteer),0],
                           [0,0,1]])

    A_phicenter = np.array([[np.cos(phicenter),-np.sin(phicenter),0],
                      [np.sin(phicenter),np.cos(phicenter),0],
                      [0,0,1]])

    grid["xyz"] = np.matmul(A_thetasteer,grid["xyz"].T).T
    grid["xyz"] = np.matmul(A_phisteer,grid["xyz"].T).T

    grid["xyz"][:,0] = grid["xyz"][:,0] + Rcenter
    grid["xyz"][:,2] = grid["xyz"][:,2] + zcenter

    grid["xyz"] = np.matmul(A_phicenter,grid["xyz"].T).T

    grid["dxdydz"] = np.matmul(A_thetasteer,grid["dxdydz"].T).T
    grid["dxdydz"] = np.matmul(A_phisteer,grid["dxdydz"].T).T
    grid["dxdydz"] = np.matmul(A_phicenter,grid["dxdydz"].T).T
    
    return grid


