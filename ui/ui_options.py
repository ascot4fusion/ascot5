import h5py
import numpy as np


class ui_optionsIO:

    ##################################
    # SIMULATION MODE AND TIME STEPS #
    ##################################

    # SIM_MODE == 1, 2, or 3
    # - Simulation mode: 
    #   1 - Gyro-orbit 
    #   2 - Guiding center
    #   3 - Hybrid !!TODO Not implemented!!
    #   4 - Magnetic field lines !!TODO Not implemented!!
    #
    # ENABLE_ADAPTIVE == 0, 1, or 2
    # - Use adaptive time-step in guiding center and hybrid 
    #   (only in the guiding center part) modes
    #   0 - Use fixed time-step
    #   1 - Use adaptive time-step
    #   2 - Use acceleration (valid only in axisymmetric fields) !!TODO Not implemented!!
    #

    SIM_MODE        = 1
    ENABLE_ADAPTIVE = 0
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
    #   N is this parameter !!TODO not implemented!!
    #

    FIXEDSTEP_USE_USERDEFINED     = 1
    FIXEDSTEP_USERDEFINED         = 1e-8
    FIXEDSTEP_NSTEPS_PER_GYROTIME = 10


    ## Adaptive time-step specific options ##
    # ADAPTIVE_TOL_ORBIT == 1e-8
    # - Relative error tolerance for orbit following
    #
    # ADAPTIVE_TOL_CCOL == 1e-2
    # - Relative error tolerance for Coulomb collisions
    # 
    # ADAPTIVE_MAX_DRHO == 1.0 !!TODO Not implemented!!
    # - Maximum allowed change in rho during one time-step
    #
    # ADAPTIVE_MAX_DPHI == 1.0 !!TODO Not implemented!!
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
    # ENDCOND_CPUTIMELIM == 0 !!TODO Not implemented!!
    # - Terminate when marker specific cpu time is reached
    #
    # ENDCOND_RHOLIM == 0 !!TODO Not implemented!!
    # - Terminate when marker goes below or over user-specified 
    #   rho (or stellarator equivalent) values
    #
    # ENDCOND_ENERGYLIM  == 1
    # - Terminate when marker energy is below a user-specified value
    #
    # ENDCOND_WALLHIT == 1
    # - Terminate when marker impacts wall
    #
    # ENDCOND_MAXORBS == 0 !!TODO Not implemented!!
    # - Terminate when marker has completed user-specified number of
    #   toroidal or poloidal orbits
    #

    ENDCOND_SIMTIMELIM = 1
    ENDCOND_CPUTIMELIM = 0
    ENDCOND_RHOLIM     = 0
    ENDCOND_ENERGYLIM  = 0
    ENDCOND_WALLHIT    = 0
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

    ENDCOND_MAX_SIM_TIME             = 1.0e-5
    ENDCOND_MAX_CPU_TIME             = 1.0e5
    ENDCOND_MAX_RHO                  = 1.0
    ENDCOND_MIN_RHO                  = 0.0
    ENDCOND_MIN_ENERGY               = 10.0
    ENDCOND_MIN_ENERGY_TIMES_THERMAL = 2.0
    ENDCOND_MAX_TOROIDALORBS         = 1000
    ENDCOND_MAX_POLOIDALORBS         = 500


    ###########
    # PHYSICS #
    ###########

    ## Define what physics are included in simulations ##
    # ENABLE_ORBIT_FOLLOWING == 1
    # - Trace particles in an electromagnetic field
    #
    # ENABLE_COULOMB_COLLISIONS == 0
    # - Particles experience Coulomb collisions with background 
    #   plasma
    #

    ENABLE_ORBIT_FOLLOWING    = 1
    ENABLE_COULOMB_COLLISIONS = 1


    ###############
    # DIAGNOSTICS #
    ###############

    ## Distribution specific options ##
    # ENABLE_RZVparaVperp_DIST == 0 or 1
    # - Enables collection of marker histogram in coordinates
    #     R - major radius in [m]
    #     z - z-coordinate in [m]
    #     vpara - velocity component parallel to magnetic field [m/s]
    #     vperp - velocity component perpendicular to magnetic field [m/s]
    #
    # Distributions are defined with parameters
    # MIN_X - Smallest bin edge for coordinate X
    # MAX_X - Largest bin edge for coordinate X
    # BIN_X - Number of bins interval [MIN_X MAX_X] is divided into
    #

    ENABLE_RZVparaVperp_DIST    = 0

    DIST_RZVparaVperp_MIN_R     = 3.0
    DIST_RZVparaVperp_MAX_R     = 8.5
    DIST_RZVparaVperp_BIN_R     = 20

    DIST_RZVparaVperp_MIN_Z     = -4.25
    DIST_RZVparaVperp_MAX_Z     = 3.6
    DIST_RZVparaVperp_BIN_Z     = 40

    DIST_RZVparaVperp_MIN_VPARA = -1.5e7
    DIST_RZVparaVperp_MAX_VPARA = 1.5e7
    DIST_RZVparaVperp_BIN_VPARA = 40

    DIST_RZVparaVperp_MIN_VPERP = 0
    DIST_RZVparaVperp_MAX_VPERP = 1.5e7
    DIST_RZVparaVperp_BIN_VPERP = 20


    ## Orbit writing specific options ##
    # ENABLE_ORBITWRITE == 0 or 1 !!TODO Not implemented!!
    # - Write exact marker state
    #
    # ORBITWRITE_MODE == 0, 1, or 2
    # - When particle state is written
    #   0 - When particle crosses a plane (Poincare-plot)
    #   1 - Between given time intervals
    #   2 - Keep last positions stored and write them
    #       at the end of simulation
    #
    # ORBITWRITE_NTOROIDALPLOTS == 1
    # - Number of toroidal plots collected
    #
    # ORBITWRITE_TOROIDALANGLES == [0 180]
    # - Poloidal angles defining the toroidal planes where
    #   toroidal plots are collected
    #
    # ORBITWRITE_NPOLOIDALPLOTS == 1
    # - Number of poloidal plots collected
    #
    # ORBITWRITE_POLOIDALANGLES == [0 180]
    # - Toroidal angles defining the poloidal planes where
    #   poloidal plots are collected
    #
    # ORBITWRITE_INTERVAL == 1e-6
    # - Time interval [s] for writing particle state
    #
    # ORBITWRITE_LASTNPOINTS == 100
    # - Keep last N positions stored and write them when
    #   marker is terminated
    #

    ENABLE_ORBITWRITE         = 1
    ORBITWRITE_MODE           = 1
    ORBITWRITE_NTOROIDALPLOTS = 1
    ORBITWRITE_TOROIDALANGLES = np.array([0, 180])
    ORBITWRITE_NPOLOIDALPLOTS = 1
    ORBITWRITE_POLOIDALANGLES = np.array([0, 180])
    ORBITWRITE_INTERVAL       = 1e-6
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

##########################################################
## A python script for writing the options in hdf5 file ##
##    DO NOT MODIFY THIS IF YOU ARE NOT A DEVELOPER     ##


# Writes options to given hdf5-file
# c -options class (defined at the start of this file)
# fn -filename for hdf5-file
def writeHdf5(c,fn):
    f = h5py.File(fn, "a")
    o = f.create_group("options")

    o.create_dataset("SIM_MODE", data = c.SIM_MODE, dtype='i4')
    o.create_dataset("ENABLE_ADAPTIVE", data = c.ENABLE_ADAPTIVE, dtype='i4')
    o.create_dataset("RECORD_GO_AS_GC", data = c.RECORD_GO_AS_GC, dtype='i4')
    o.create_dataset("FIXEDSTEP_USE_USERDEFINED", data = c.FIXEDSTEP_USE_USERDEFINED, dtype='i4')
    o.create_dataset("FIXEDSTEP_USERDEFINED", data = c.FIXEDSTEP_USERDEFINED, dtype='f8')
    o.create_dataset("FIXEDSTEP_NSTEPS_PER_GYROTIME", data = c.FIXEDSTEP_NSTEPS_PER_GYROTIME, dtype='i4')
    o.create_dataset("ADAPTIVE_TOL_ORBIT", data = c.ADAPTIVE_TOL_ORBIT, dtype='f8')
    o.create_dataset("ADAPTIVE_TOL_CCOL", data = c.ADAPTIVE_TOL_CCOL, dtype='f8')
    o.create_dataset("ADAPTIVE_MAX_DRHO", data = c.ADAPTIVE_MAX_DRHO, dtype='f8')
    o.create_dataset("ADAPTIVE_MAX_DPHI", data = c.ADAPTIVE_MAX_DPHI, dtype='f8')
    o.create_dataset("ADAPTIVE_MAX_ACC", data = c.ADAPTIVE_MAX_ACC, dtype='f8')
    o.create_dataset("ENDCOND_SIMTIMELIM", data = c.ENDCOND_SIMTIMELIM, dtype='i4')
    o.create_dataset("ENDCOND_CPUTIMELIM", data = c.ENDCOND_CPUTIMELIM, dtype='i4')
    o.create_dataset("ENDCOND_RHOLIM", data = c.ENDCOND_RHOLIM, dtype='i4')
    o.create_dataset("ENDCOND_ENERGYLIM", data = c.ENDCOND_ENERGYLIM, dtype='i4')
    o.create_dataset("ENDCOND_WALLHIT", data = c.ENDCOND_WALLHIT, dtype='i4')
    o.create_dataset("ENDCOND_MAXORBS", data = c.ENDCOND_MAXORBS, dtype='i4')
    o.create_dataset("ENDCOND_MAX_SIM_TIME", data = c.ENDCOND_MAX_SIM_TIME, dtype='f8')
    o.create_dataset("ENDCOND_MAX_CPU_TIME", data = c.ENDCOND_MAX_CPU_TIME, dtype='f8')
    o.create_dataset("ENDCOND_MAX_RHO", data = c.ENDCOND_MAX_RHO, dtype='f8')
    o.create_dataset("ENDCOND_MIN_RHO", data = c.ENDCOND_MIN_RHO, dtype='f8')
    o.create_dataset("ENDCOND_MIN_ENERGY", data = c.ENDCOND_MIN_ENERGY, dtype='f8')
    o.create_dataset("ENDCOND_MIN_ENERGY_TIMES_THERMAL", data = c.ENDCOND_MIN_ENERGY_TIMES_THERMAL, dtype='f8')
    o.create_dataset("ENDCOND_MAX_TOROIDALORBS", data = c.ENDCOND_MAX_TOROIDALORBS, dtype='i4')
    o.create_dataset("ENDCOND_MAX_POLOIDALORBS", data = c.ENDCOND_MAX_POLOIDALORBS, dtype='i4')
    o.create_dataset("ENABLE_ORBIT_FOLLOWING", data = c.ENABLE_ORBIT_FOLLOWING, dtype='i4')
    o.create_dataset("ENABLE_COULOMB_COLLISIONS", data = c.ENABLE_COULOMB_COLLISIONS, dtype='i4')
    o.create_dataset("ENABLE_RZVparaVperp_DIST", data = c.ENABLE_RZVparaVperp_DIST, dtype='i4')
    o.create_dataset("DIST_RZVparaVperp_MIN_R", data = c.DIST_RZVparaVperp_MIN_R, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_MAX_R", data = c.DIST_RZVparaVperp_MAX_R, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_BIN_R", data = c.DIST_RZVparaVperp_BIN_R, dtype='i4')
    o.create_dataset("DIST_RZVparaVperp_MIN_Z", data = c.DIST_RZVparaVperp_MIN_Z, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_MAX_Z", data = c.DIST_RZVparaVperp_MAX_Z, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_BIN_Z", data = c.DIST_RZVparaVperp_BIN_Z, dtype='i4')
    o.create_dataset("DIST_RZVparaVperp_MIN_VPARA", data = c.DIST_RZVparaVperp_MIN_VPARA, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_MAX_VPARA", data = c.DIST_RZVparaVperp_MAX_VPARA, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_BIN_VPARA", data = c.DIST_RZVparaVperp_BIN_VPARA, dtype='i4')
    o.create_dataset("DIST_RZVparaVperp_MIN_VPERP", data = c.DIST_RZVparaVperp_MIN_VPERP, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_MAX_VPERP", data = c.DIST_RZVparaVperp_MAX_VPERP, dtype='f8')
    o.create_dataset("DIST_RZVparaVperp_BIN_VPERP", data = c.DIST_RZVparaVperp_BIN_VPERP, dtype='i4')
    o.create_dataset("ENABLE_ORBITWRITE", data = c.ENABLE_ORBITWRITE, dtype='i4')
    o.create_dataset("ORBITWRITE_MODE", data = c.ORBITWRITE_MODE, dtype='i4')
    o.create_dataset("ORBITWRITE_NTOROIDALPLOTS", data = c.ORBITWRITE_NTOROIDALPLOTS, dtype='i4')
    o.create_dataset("ORBITWRITE_TOROIDALANGLES", data = c.ORBITWRITE_TOROIDALANGLES, dtype='f8')
    o.create_dataset("ORBITWRITE_NPOLOIDALPLOTS", data = c.ORBITWRITE_NPOLOIDALPLOTS, dtype='i4')
    o.create_dataset("ORBITWRITE_POLOIDALANGLES", data = c.ORBITWRITE_POLOIDALANGLES, dtype='f8')
    o.create_dataset("ORBITWRITE_INTERVAL", data = c.ORBITWRITE_INTERVAL, dtype='f8')
    o.create_dataset("ORBITWRITE_LASTNPOINTS", data = c.ORBITWRITE_LASTNPOINTS, dtype='i4')
    o.create_dataset("ENABLE_DEBUGDIST", data = c.ENABLE_DEBUGDIST, dtype='i4')


    f.close()

if __name__ == "__main__":
    writeHdf5(ui_optionsIO(), "ascot.h5")
