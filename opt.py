"""
Human readable options file for generating ASCOT5 options.

For users: Modify the values here (in class opt) and use
script a5setoptions (in a5py) to assing these options to 
your HDF5 file.

For developers: Each time you add new options, make sure
you define its type in generateopt function below.
"""
import h5py
import numpy as np


class opt():

    ##################################
    # SIMULATION MODE AND TIME STEPS #
    ##################################

    # SIM_MODE == 1, 2, or 3
    # - Simulation mode: 
    #   1 - Gyro-orbit 
    #   2 - Guiding center
    #   3 - Hybrid
    #   4 - Magnetic field lines
    #
    # ENABLE_ADAPTIVE == 0, 1, or 2
    # - Use adaptive time-step in guiding center and hybrid 
    #   (only in the guiding center part of the simulation) modes
    #   0 - Use fixed time-step
    #   1 - Use adaptive time-step
    #   2 - Use acceleration (valid only in axisymmetric fields) !!TODO Not implemented!!
    #

    SIM_MODE        = 3
    ENABLE_ADAPTIVE = 1
    RECORD_GO_AS_GC = 0


    ## Fixed time-step specific options ##
    #
    # FIXEDSTEP_USE_USERDEFINED == 0 or 1
    # - Use user-defined time step instead of one based on 
    #   gyro-frequency (see below)
    # 
    # FIXEDSTEP_USERDEFINED == 1e-8
    # - User-defined fixed time-step in seconds
    #
    # FIXEDSTEP_NSTEPS_PER_GYROTIME == 10
    # - Time-step determined from dt = (1/gyrofreq)/N where 
    #   N is this parameter 
    #

    FIXEDSTEP_USE_USERDEFINED     = 1
    FIXEDSTEP_USERDEFINED         = 1e-10
    FIXEDSTEP_NSTEPS_PER_GYROTIME = 50


    ## Adaptive time-step specific options ##
    # ADAPTIVE_TOL_ORBIT == 1e-8
    # - Relative error tolerance for orbit following
    #
    # ADAPTIVE_TOL_CCOL == 1e-2
    # - Relative error tolerance for Coulomb collisions
    # 
    # ADAPTIVE_MAX_DRHO == 1.0
    # - Maximum allowed change in rho during one time-step
    #
    # ADAPTIVE_MAX_DPHI == 1.0
    # - Maximum allowed change in toroidal angle during
    #   one time-step
    #
    # ADAPTIVE_MAX_ACC == 10.0 !!TODO Not implemented!!
    # - Maximum acceleration factor
    #

    ADAPTIVE_TOL_ORBIT = 1e-8
    ADAPTIVE_TOL_CCOL  = 1e-2
    ADAPTIVE_MAX_DRHO  = 1.0
    ADAPTIVE_MAX_DPHI  = 1.0
    ADAPTIVE_MAX_ACC   = 10.0

    ##################
    # END CONDITIONS #
    ##################

    ## Enable/Disable end conditions
    # ENDCOND_SIMTIMELIM == 1
    # - Terminate when maximum simulation (marker specific) time
    #   is reached
    #
    # ENDCOND_CPUTIMELIM == 0 
    # - Terminate when marker specific cpu time is reached
    #
    # ENDCOND_RHOLIM == 0
    # - Terminate when marker goes below or over user-specified 
    #   rho (or stellarator equivalent) values
    #
    # ENDCOND_ENERGYLIM  == 1
    # - Terminate when marker energy is below a user-specified value
    #
    # ENDCOND_WALLHIT == 1
    # - Terminate when marker impacts wall
    #
    # ENDCOND_MAXORBS == 0
    # - Terminate when marker has completed user-specified number of
    #   toroidal or poloidal orbits
    #

    ENDCOND_SIMTIMELIM = 1
    ENDCOND_CPUTIMELIM = 1
    ENDCOND_RHOLIM     = 0
    ENDCOND_ENERGYLIM  = 1
    ENDCOND_WALLHIT    = 1
    ENDCOND_MAXORBS    = 0

    ## Defining values for different end conditions
    # ENDCOND_MAX_SIM_TIME == 1.0
    # - Maximum simulation time in seconds
    #
    # ENDCOND_MAX_CPU_TIME == 1.0e5
    # - Maximum cpu time allowed for each marker
    #
    # ENDCOND_MAX_RHO == 1.0
    # - Maximum rho value
    #
    # ENDCOND_MIN_RHO == 0.0
    # - Minimum rho value
    #
    # ENDCOND_MIN_ENERGY == 10.0
    # - Minimum energy in electron volts
    #
    # ENDCOND_MIN_ENERGY_TIMES_THERMAL == 2.0
    # - Minimum energy determined from N x local 
    #   (electron) thermal energy where N is this parameter
    #
    # ENDCOND_MAX_TOROIDALORBS == 1000
    # - Maximum number of toroidal orbits
    #
    # ENDCOND_MAX_POLOIDALORBS == 500
    # - Maximum number of poloidal orbits
    #

    ENDCOND_MAX_SIM_TIME             = 1.0e0
    ENDCOND_MAX_CPU_TIME             = 1.0e4
    ENDCOND_MAX_RHO                  = 1.0
    ENDCOND_MIN_RHO                  = 0.0
    ENDCOND_MIN_ENERGY               = 10.0e3
    ENDCOND_MIN_ENERGY_TIMES_THERMAL = 2.0
    ENDCOND_MAX_TOROIDALORBS         = 1000
    ENDCOND_MAX_POLOIDALORBS         = 500


    ###########
    # PHYSICS #
    ###########

    ## Define what physics are included in simulations ##
    # ENABLE_ORBIT_FOLLOWING == 1
    # - Trace markers in an electromagnetic field
    #
    # ENABLE_COULOMB_COLLISIONS == 0
    # - Markers experience Coulomb collisions with background 
    #   plasma
    #

    ENABLE_ORBIT_FOLLOWING    = 1
    ENABLE_COULOMB_COLLISIONS = 1


    ###############
    # DIAGNOSTICS #
    ###############

    ## Distribution specific options ##
    # ENABLE_R_phi_z_vpa_vpe_t_q_DIST == 0 or 1
    # - Enables collection of marker distribution histogram in coordinates
    #     R     - major radius in meters [m]
    #     phi   - toroidal angle in degrees [deg]
    #     z     - z-coordinate in meters [m]
    #     vpara - velocity component parallel to magnetic field in meters per second [m/s]
    #     vperp - velocity component perpendicular to magnetic field in meters per second [m/s]
    #     t     - time in seconds [s]
    #     q     - charge in multiples of elementary charge [e] (must be an integer value)
    #
    # ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST == 0 or 1 !!TODO Not implemented!!
    # - Enables collection of marker distribution histogram in coordinates
    #     R     - major radius in meters [m]
    #     phi   - toroidal angle in degrees [deg]
    #     z     - z-coordinate in meters [m]
    #     vR    - velocity R-component in meters per second [m/s]
    #     vphi  - velocity phi-component in meters per second [m/s]
    #     vz    - velocity z-component in meters per second [m/s]
    #     t     - time in seconds [s]
    #     q     - charge in multiples of elementary charge [e] (must be an integer value)
    #
    # Distribution histograms are defined with parameters
    # MIN_X  - Smallest bin edge for coordinate X
    # MAX_X  - Largest bin edge for coordinate X
    # NBIN_X - Number of bins interval [MIN_X MAX_X] is divided into
    #

    ENABLE_R_phi_z_vpa_vpe_t_q_DIST    = 1
    ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST = 0

    DIST_MIN_R    = 3.0
    DIST_MAX_R    = 8.5
    DIST_NBIN_R   = 20

    DIST_MIN_phi  = 0
    DIST_MAX_phi  = 360
    DIST_NBIN_phi = 1

    DIST_MIN_z    = -4.25
    DIST_MAX_z    = 3.6
    DIST_NBIN_z   = 40

    DIST_MIN_vpa  = -1.5e7
    DIST_MAX_vpa  = 1.5e7
    DIST_NBIN_vpa = 40

    DIST_MIN_vpe  = 0
    DIST_MAX_vpe  = 1.5e7
    DIST_NBIN_vpe = 20

    DIST_MIN_vR    = -1.5e7
    DIST_MAX_vR    = 1.5e7
    DIST_NBIN_vR   = 40

    DIST_MIN_vphi  = -1.5e7
    DIST_MAX_vphi  = 1.5e7
    DIST_NBIN_vphi = 40

    DIST_MIN_vz    = -1.5e7
    DIST_MAX_vz    = 1.5e7
    DIST_NBIN_vz   = 40

    DIST_MIN_t    = 0
    DIST_MAX_t    = 100
    DIST_NBIN_t   = 1

    DIST_MIN_q    = -100
    DIST_MAX_q    = 100
    DIST_NBIN_q   = 1


    ## Orbit writing specific options ##
    # ENABLE_ORBITWRITE == 0 or 1
    # - Write exact marker state
    #
    # ORBITWRITE_MODE == 0, 1, or 2
    # - When marker state is written
    #   0 - When marker crosses a plane (Poincare-plot)
    #   1 - Between given time intervals
    #   2 - Keep last positions stored and write them
    #       at the end of simulation
    #
    # ORBITWRITE_NTOROIDALPLOTS == 1
    # - Number of toroidal plots collected
    #
    # ORBITWRITE_TOROIDALANGLES == [0, 180]
    # - Poloidal angles defining the toroidal planes where
    #   toroidal plots are collected
    #
    # ORBITWRITE_NPOLOIDALPLOTS == 1
    # - Number of poloidal plots collected
    #
    # ORBITWRITE_POLOIDALANGLES == [0, 180]
    # - Toroidal angles defining the poloidal planes where
    #   poloidal plots are collected
    #
    # ORBITWRITE_INTERVAL == 1e-6
    # - Time interval [s] for writing marker state
    #
    # ORBITWRITE_LASTNPOINTS == 100
    # - Keep last N positions stored and write them when
    #   marker is terminated
    #

    ENABLE_ORBITWRITE         = 0
    ORBITWRITE_MODE           = 2
    ORBITWRITE_NTOROIDALPLOTS = 2
    ORBITWRITE_TOROIDALANGLES = [0, 180]
    ORBITWRITE_NPOLOIDALPLOTS = 2
    ORBITWRITE_POLOIDALANGLES = [0, 180]
    ORBITWRITE_INTERVAL       = 1e-9
    ORBITWRITE_LASTNPOINTS    = 100

    ## Debug options ##
    # ENABLE_DEBUGDIST == 0 or 1 !!TODO Not implemented!!
    # - Collect debug histograms: 
    #   - time-step distribution
    #   - reason for time-step rejection distribution
    #   - acceleration factor distribution
    #

    ENABLE_DEBUGDIST = 0

## END OF OPTIONS ##


#############################################################
## A python functions for writing the options in HDF5 file ##
##     DO NOT MODIFY THIS IF YOU ARE NOT A DEVELOPER       ##

def generateopt():
    """
    Convert human readable options to python dictionary.

    The options to be converted is the one specified above.
    All fields should be converted to numpy arrays of specific 
    type.
    """

    # This ugly thing converts class variables into dictionary
    f = {key:value for key, value in opt.__dict__.items() if not key.startswith('__') and not callable(key)}
    
    ## Simulation mode ##
    f["SIM_MODE"]        = settype(f["SIM_MODE"],'i4')
    f["ENABLE_ADAPTIVE"] = settype(f["ENABLE_ADAPTIVE"],'i4')
    f["RECORD_GO_AS_GC"] = settype(f["RECORD_GO_AS_GC"],'i4')

    ## Fixed time-step specific options ##
    f["FIXEDSTEP_USE_USERDEFINED"]     = settype(f["FIXEDSTEP_USE_USERDEFINED"],'i4')
    f["FIXEDSTEP_USERDEFINED"]         = settype(f["FIXEDSTEP_USERDEFINED"],'f8')
    f["FIXEDSTEP_NSTEPS_PER_GYROTIME"] = settype(f["FIXEDSTEP_NSTEPS_PER_GYROTIME"],'i4')

    ## Adaptive time-step specific options ##
    f["ADAPTIVE_TOL_ORBIT"] = settype(f["ADAPTIVE_TOL_ORBIT"],'f8')
    f["ADAPTIVE_TOL_CCOL"]  = settype(f["ADAPTIVE_TOL_CCOL"],'f8')
    f["ADAPTIVE_MAX_DRHO"]  = settype(f["ADAPTIVE_MAX_DRHO"],'f8')
    f["ADAPTIVE_MAX_DPHI"]  = settype(f["ADAPTIVE_MAX_DPHI"],'f8')
    f["ADAPTIVE_MAX_ACC"]   = settype(f["ADAPTIVE_MAX_ACC"],'f8')

    ## End conditions ##
    f["ENDCOND_SIMTIMELIM"] = settype(f["ENDCOND_SIMTIMELIM"],'i4')
    f["ENDCOND_CPUTIMELIM"] = settype(f["ENDCOND_CPUTIMELIM"],'i4')
    f["ENDCOND_RHOLIM"]     = settype(f["ENDCOND_RHOLIM"],'i4')
    f["ENDCOND_ENERGYLIM"]  = settype(f["ENDCOND_ENERGYLIM"],'i4')
    f["ENDCOND_WALLHIT"]    = settype(f["ENDCOND_WALLHIT"],'i4')
    f["ENDCOND_MAXORBS"]    = settype(f["ENDCOND_MAXORBS"],'i4')

    f["ENDCOND_MAX_SIM_TIME"]             = settype(f["ENDCOND_MAX_SIM_TIME"],'f8')
    f["ENDCOND_MAX_CPU_TIME"]             = settype(f["ENDCOND_MAX_CPU_TIME"],'f8')
    f["ENDCOND_MAX_RHO"]                  = settype(f["ENDCOND_MAX_RHO"],'f8')
    f["ENDCOND_MIN_RHO"]                  = settype(f["ENDCOND_MIN_RHO"],'f8')
    f["ENDCOND_MIN_ENERGY"]               = settype(f["ENDCOND_MIN_ENERGY"],'f8')
    f["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = settype(f["ENDCOND_MIN_ENERGY_TIMES_THERMAL"],'f8')
    f["ENDCOND_MAX_TOROIDALORBS"]         = settype(f["ENDCOND_MAX_TOROIDALORBS"],'i4')
    f["ENDCOND_MAX_POLOIDALORBS"]         = settype(f["ENDCOND_MAX_POLOIDALORBS"],'i4')

    ## Physics ##
    f["ENABLE_ORBIT_FOLLOWING"]    = settype(f["ENABLE_ORBIT_FOLLOWING"],'i4')
    f["ENABLE_COULOMB_COLLISIONS"] = settype(f["ENABLE_COULOMB_COLLISIONS"],'i4')

    ## Distributions ##
    f["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"]   = settype(f["ENABLE_R_phi_z_vpa_vpe_t_q_DIST"],'i4')

    f["ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST"]    = settype(f["ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST"],'i4')
    
    f["DIST_MIN_R"]    = settype(f["DIST_MIN_R"],'f8')
    f["DIST_MAX_R"]    = settype(f["DIST_MAX_R"],'f8')
    f["DIST_NBIN_R"]   = settype(f["DIST_NBIN_R"],'i4')
    
    f["DIST_MIN_phi"]  = settype(f["DIST_MIN_phi"],'f8')
    f["DIST_MAX_phi"]  = settype(f["DIST_MAX_phi"],'f8')
    f["DIST_NBIN_phi"] = settype(f["DIST_NBIN_phi"],'i4')

    f["DIST_MIN_z"]    = settype(f["DIST_MIN_z"],'f8')
    f["DIST_MAX_z"]    = settype(f["DIST_MAX_z"],'f8')
    f["DIST_NBIN_z"]   = settype(f["DIST_NBIN_z"],'i4')

    f["DIST_MIN_vpa"]  = settype(f["DIST_MIN_vpa"],'f8')
    f["DIST_MAX_vpa"]  = settype(f["DIST_MAX_vpa"],'f8')
    f["DIST_NBIN_vpa"] = settype(f["DIST_NBIN_vpa"],'i4')

    f["DIST_MIN_vpe"]  = settype(f["DIST_MIN_vpe"],'f8')
    f["DIST_MAX_vpe"]  = settype(f["DIST_MAX_vpe"],'f8')
    f["DIST_NBIN_vpe"] = settype(f["DIST_NBIN_vpe"],'i4')

    f["DIST_MIN_t"]    = settype(f["DIST_MIN_t"],'f8')
    f["DIST_MAX_t"]    = settype(f["DIST_MAX_t"],'f8')
    f["DIST_NBIN_t"]   = settype(f["DIST_NBIN_t"],'i4')

    f["DIST_MIN_q"]    = settype(f["DIST_MIN_q"],'i4')
    f["DIST_MAX_q"]    = settype(f["DIST_MAX_q"],'i4')
    f["DIST_NBIN_q"]   = settype(f["DIST_NBIN_q"],'i4')
    
    f["DIST_MIN_vR"]    = settype(f["DIST_MIN_vR"],'f8')
    f["DIST_MAX_vR"]    = settype(f["DIST_MAX_vR"],'f8')
    f["DIST_NBIN_vR"]   = settype(f["DIST_NBIN_vR"],'i4')

    f["DIST_MIN_vphi"]  = settype(f["DIST_MIN_vphi"],'f8')
    f["DIST_MAX_vphi"]  = settype(f["DIST_MAX_vphi"],'f8')
    f["DIST_NBIN_vphi"] = settype(f["DIST_NBIN_vphi"],'i4')

    f["DIST_MIN_vz"]    = settype(f["DIST_MIN_vz"],'f8')
    f["DIST_MAX_vz"]    = settype(f["DIST_MAX_vz"],'f8')
    f["DIST_NBIN_vz"]   = settype(f["DIST_NBIN_vz"],'i4')

    ## Orbit writing specific options ##
    f["ENABLE_ORBITWRITE"]         = settype(f["ENABLE_ORBITWRITE"],'i4')
    f["ORBITWRITE_MODE"]           = settype(f["ORBITWRITE_MODE"],'i4')
    f["ORBITWRITE_NTOROIDALPLOTS"] = settype(f["ORBITWRITE_NTOROIDALPLOTS"],'i4')
    f["ORBITWRITE_TOROIDALANGLES"] = settype(f["ORBITWRITE_TOROIDALANGLES"],'f8')
    f["ORBITWRITE_NPOLOIDALPLOTS"] = settype(f["ORBITWRITE_NPOLOIDALPLOTS"],'i4')
    f["ORBITWRITE_POLOIDALANGLES"] = settype(f["ORBITWRITE_POLOIDALANGLES"],'f8')
    f["ORBITWRITE_INTERVAL"]       = settype(f["ORBITWRITE_INTERVAL"],'f8')
    f["ORBITWRITE_LASTNPOINTS"]    = settype(f["ORBITWRITE_LASTNPOINTS"],'i4')

    ## Debug options ##
    f["ENABLE_DEBUGDIST"] = settype(f["ENABLE_DEBUGDIST"],'i4')

    return f

def settype(data, type_):
    """
    Convert input into numpy array of specific type.
    """
    data = np.array(data)
    data = data.astype(type_)
    return data

