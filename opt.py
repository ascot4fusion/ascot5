"""
Human readable options file for generating ASCOT5 options.

Users can modify the values here in class opt.opt which are written to the HDF5
file by calling a5setoptions (in a5py).

Developers need to define type for the new options parameter in generateopt()
each time one is added here.

File: opt.py
"""
import h5py
import numpy as np


class opt():
    """
    Class for storing option parameters.

    Feel free to change these values. Use generateopt() to generate dictionary
    from this class which can then be written to HDF5 file.
    """

    #**************************************************************************#
    #*                      SIMULATION MODE AND TIME-STEP                     *#
    #*                                                                        *#
    #**************************************************************************#

    ## Simulation mode (1, 2, 3, 4)
    #
    # - 1 Gyro-orbit
    # - 2 Guiding center
    # - 3 Hybrid
    # - 4 Magnetic field lines
    #
    SIM_MODE = 2


    ## Use adaptive time-step (0, 1)
    #
    # This option is used only if opt.opt.SIM_MODE = 2 or 3. Gyro-orbit
    # simulations are always done with fixed time-step and magnetic field line
    # simulations with adaptive time-step.
    #
    # - 0 Use fixed time-step
    # - 1 Use adaptive time-step
    #
    ENABLE_ADAPTIVE = 1


    ## Record GOs as GCs in diagnostics (0, 1)
    #
    # - 0 Record GOs as GOs
    # - 1 Record GOs as GCs
    #
    RECORD_GO_AS_GC = 0


    ## Define fixed time-step value explicitly (0,1)
    #
    # Note: The adaptive scheme uses fixed time-step value as an initial step.
    #
    # - 0 Calculate time-step from opt.opt.FIXEDSTEP_NSTEPS_PER_GYROTIME
    # - 1 Use opt.opt.FIXEDSTEP_USERDEFINED as a time-step
    #
    FIXEDSTEP_USE_USERDEFINED = 1


    ## User-defined time-step [s]
    FIXEDSTEP_USERDEFINED = 1e-10


    ## Time-step is 2pi / ( gyrofrequency * N ) where N is this parameter
    FIXEDSTEP_NSTEPS_PER_GYROTIME = 50


    ## Relative error tolerance for orbit following in adaptive scheme
    ADAPTIVE_TOL_ORBIT = 1e-8


    ## Relative error tolerance for Coulomb collisions in adaptive scheme
    ADAPTIVE_TOL_CCOL = 1e-2


    ## Maximum allowed change in rho during one time-step in adaptive scheme
    ADAPTIVE_MAX_DRHO = 1.0


    ## Maximum allowed change in phi during one time-step in adaptive scheme
    ADAPTIVE_MAX_DPHI = 1.0


    #**************************************************************************#
    #*                             END CONDITIONS                             *#
    #*                                                                        *#
    #**************************************************************************#

    ## Terminate when marker's clock ("laboratory") time reaches a limit
    #
    # The limit is set by opt.opt.ENDCOND_MAX_SIM_TIME
    #
    ENDCOND_SIMTIMELIM = 1


    ## Terminate when marker's cpu time reaches a limit
    #
    # The limit is set by opt.opt.ENDCOND_MAX_CPU_TIME
    #
    ENDCOND_CPUTIMELIM = 1


    ## Terminate if marker goes outside given rho boundaries
    #
    # rho boundaries are defined by opt.opt.ENDCOND_MAX_RHO and
    # opt.opt.ENDCOND_MAX_RHO.
    #
    ENDCOND_RHOLIM = 1


    ## Terminate when marker energy is below a user-specified value
    #
    # The user specified values are opt.opt.ENDCOND_MIN_ENERGY and
    # opt.opt.ENDCOND_MIN_ENERGY_TIMES_THERMAL. Marker is terminated when either
    # of these limits is reached.
    #
    ENDCOND_ENERGYLIM = 1


    ## Terminate when marker impacts wall
    ENDCOND_WALLHIT = 0


    ## Terminate when marker has completed user-specified number of orbits
    #
    # Limit opt.opt.ENDCOND_MAX_TOROIDALORBS is used for a number of toroidal
    # and opt.opt.ENDCOND_MAX_POLOIDALORBS for poloidal orbits. Marker is
    # terminated when either of these limits is reached.
    #
    ENDCOND_MAXORBS = 1


    ## Maximum simulation time [s]
    ENDCOND_MAX_SIM_TIME = 1.0e-5


    ## Maximum cpu time [s]
    ENDCOND_MAX_CPU_TIME = 1.0e4


    ## Maximum rho value
    ENDCOND_MAX_RHO = 1.0


    ## Minimum rho value
    ENDCOND_MIN_RHO = 0.0


    ## Minimum energy [eV]
    ENDCOND_MIN_ENERGY = 1.0e3


    ## Minimum energy limit local electron thermal energy times this value
    ENDCOND_MIN_ENERGY_TIMES_THERMAL = 2.0


    ## Maximum number of toroidal orbits
    ENDCOND_MAX_TOROIDALORBS = 3


    ## Maximum number of poloidal orbits
    ENDCOND_MAX_POLOIDALORBS = 500

    #**************************************************************************#
    #*                               PHYSICS                                  *#
    #*                                                                        *#
    #**************************************************************************#

    ## Trace markers in an electromagnetic field
    ENABLE_ORBIT_FOLLOWING    = 1

    ## Markers experience Coulomb collisions with background plasma
    ENABLE_COULOMB_COLLISIONS = 1

    ## Markers interact with MHD modes
    ENABLE_MHD = 1

    ## Disable first order guiding center transformation in velocity space
    DISABLE_FIRSTORDER_GCTRANS = 0

    ## Disable first order guiding center transformation in velocity space
    DISABLE_ENERGY_CCOLL = 0

    ## Disable first order guiding center transformation in velocity space
    DISABLE_PITCH_CCOLL  = 0

    ## Disable first order guiding center transformation in velocity space
    DISABLE_GCDIFF_CCOLL = 0

    #**************************************************************************#
    #*                            DISTRIBUTIONS                               *#
    #*                                                                        *#
    #**************************************************************************#

    ## Collect distribution histogram in [R, phi, z, vpa, vpe, t, q]
    #
    # The coordinates are
    # - R   major radius
    # - phi toroidal angle
    # - z   z-coordinate
    # - vpa velocity component parallel to magnetic field
    # - vpe velocity component perpendicular to magnetic field
    # - t   time
    # - q   charge
    #
    ENABLE_DIST_5D = 1


    ## Collect distribution histogram in [R, phi, z, vR, vphi, vz, t, q]
    #
    # The coordinates are
    # - R    major radius
    # - phi  toroidal angle
    # - z    z-coordinate
    # - vR   velocity R-component
    # - vphi velocity phi-component
    # - vz   velocity z-component
    # - t    time
    # - q    charge
    #
    ENABLE_DIST_6D = 0


    ## Collect distribution histogram in [rho, pol, phi, vpa, vpe, t, q]
    #
    # The coordinates are
    # - rho  flux surface
    # - pol  poloidal angle
    # - phi  toroidal angle
    # - z    z-coordinate
    # - vpa  velocity component parallel to magnetic field
    # - vpe  velocity component perpendicular to magnetic field
    # - t    time
    # - q    charge
    #
    ENABLE_DIST_rho5D = 1


    ## Collect distribution histogram in [rho, pol, phi, vR, vphi, vz, t, q]
    #
    # The coordinates are
    # - rho  flux surface
    # - pol  poloidal angle
    # - phi  toroidal angle
    # - z    z-coordinate
    # - vR   velocity R-component
    # - vphi velocity phi-component
    # - vz   velocity z-component
    # - t    time
    # - q    charge
    #
    ENABLE_DIST_rho6D = 0


    ## Minimum bin edge for major R coordinate [m]
    DIST_MIN_R    = 4.5
    ## Maximum bin edge for R coordinate [m]
    DIST_MAX_R    = 6.5
    ## Number of bins the interval [DIST_MIN_R, DIST_MAX_R] is divided to
    DIST_NBIN_R   = 40


    ## Minimum bin edge for phi coordinate [deg]
    DIST_MIN_phi  = 0
    ## Maximum bin edge for phi coordinate [deg]
    DIST_MAX_phi  = 360
    ## Number of bins the interval [DIST_MIN_phi, DIST_MAX_phi] is divided to
    DIST_NBIN_phi = 40


    ## Minimum bin edge for z coordinate [m]
    DIST_MIN_z    = -1.45
    ## Maximum bin edge for z coordinate [m]
    DIST_MAX_z    = 1.45
    ## Number of bins the interval [DIST_MIN_z, DIST_MAX_z] is divided to
    DIST_NBIN_z   = 40


    ## Minimum bin edge for rho coordinate
    DIST_MIN_rho  = 0
    ## Maximum bin edge for rho coordinate
    DIST_MAX_rho  = 1
    ## Number of bins the interval [DIST_MIN_rho, DIST_MAX_rho] is divided to
    DIST_NBIN_rho = 40


    ## Minimum bin edge for pol coordinate [deg]
    DIST_MIN_pol  = 0
    ## Maximum bin edge for pol coordinate [deg]
    DIST_MAX_pol  = 360
    ## Number of bins the interval [DIST_MIN_pol, DIST_MAX_pol] is divided to
    DIST_NBIN_pol = 40


    ## Minimum bin edge for vpa coordinate [m/s]
    DIST_MIN_vpa  = -1.5e7
    ## Maximum bin edge for vpa coordinate [m/s]
    DIST_MAX_vpa  = 1.5e7
    ## Number of bins the interval [DIST_MIN_vpa, DIST_MAX_vpa] is divided to
    DIST_NBIN_vpa = 1


    ## Minimum bin edge for vpe coordinate [m/s]
    DIST_MIN_vpe  = 0
    ## Maximum bin edge for vpe coordinate [m/s]
    DIST_MAX_vpe  = 1.5e7
    ## Number of bins the interval [DIST_MIN_vpe, DIST_MAX_vpe] is divided to
    DIST_NBIN_vpe = 1


    ## Minimum bin edge for vR coordinate [m/s]
    DIST_MIN_vR    = -1.5e7
    ## Maximum bin edge for vR coordinate [m/s]
    DIST_MAX_vR    = 1.5e7
    ## Number of bins the interval [DIST_MIN_vR, DIST_MAX_vR] is divided to
    DIST_NBIN_vR   = 1


    ## Minimum bin edge for vphi coordinate [m/s]
    DIST_MIN_vphi  = -1.5e7
    ## Maximum bin edge for vphi coordinate [m/s]
    DIST_MAX_vphi  = 1.5e7
    ## Number of bins the interval [DIST_MIN_vphi, DIST_MAX_vphi] is divided to
    DIST_NBIN_vphi = 1


    ## Minimum bin edge for vz coordinate [m/s]
    DIST_MIN_vz    = -1.5e7
    ## Maximum bin edge for vz coordinate [m/s]
    DIST_MAX_vz    = 1.5e7
    ## Number of bins the interval [DIST_MIN_vz, DIST_MAX_vz] is divided to
    DIST_NBIN_vz   = 1


    ## Minimum bin edge for t coordinate [s]
    DIST_MIN_t    = 0
    ## Maximum bin edge for t coordinate [s]
    DIST_MAX_t    = 100
    ## Number of bins the interval [DIST_MIN_t, DIST_MAX_t] is divided to
    DIST_NBIN_t   = 1


    ## Minimum bin edge for q coordinate [e]
    DIST_MIN_q    = -100
    ## Maximum bin edge for q coordinate [e]
    DIST_MAX_q    = 100
    ## Number of bins the interval [DIST_MIN_q, DIST_MAX_q] is divided to
    DIST_NBIN_q   = 1


    #**************************************************************************#
    #*                             ORBIT WRITE                                *#
    #*                                                                        *#
    #**************************************************************************#

    ## Enable diagnostics that store marker orbit
    #
    # - 0 Marker orbit diagnostics are not collected
    # - 1 Marker orbit diagnostics are collected
    #
    ENABLE_ORBITWRITE         = 1


    ## What kind of marker orbit diagnostics are collected
    #
    # These are only used if opt.opt.ENABLE_ORBITWRITE is active.
    # - 0 When marker crosses a plane (Poincare-plot)
    # - 1 Between given time intervals
    #
    ORBITWRITE_MODE           = 1

    ## Maximum number of points (per marker) to be written
    #
    # If this number is exceeded when marker is being simulated, the oldest
    # points will be replaced as long as the simulation continues. Thus, this
    # parameter is effectively the number of marker's last positions that are
    # stored.
    #
    ORBITWRITE_MAXPOINTS    = 100


    ## Poloidal angles of the toroidal planes where toroidal plots are collected
    #
    # Used when opt.opt.ENABLE_ORBITWRITE = 1 and opt.opt.ORBITWRITE_MODE = 0.
    #
    ORBITWRITE_TOROIDALANGLES = [0, 180]


    ## Toroidal angles of the poloidal planes where poloidal plots are collected
    #
    # Used when opt.opt.ENABLE_ORBITWRITE = 1 and opt.opt.ORBITWRITE_MODE = 0.
    #
    ORBITWRITE_POLOIDALANGLES = [0, 180]


    ## Time interval for writing marker state [s]
    #
    # Used when opt.opt.ENABLE_ORBITWRITE = 1 and opt.opt.ORBITWRITE_MODE = 1.
    #
    ORBITWRITE_INTERVAL       = 1e-9


def generateopt():
    """
    Converts options class to python dictionary.

    The options to be converted is the one specified by class opt.opt. All
    fields in this class are converted to numpy arrays of specific type
    (otherwise h5py cannot process them correctly) before they are stored
    into the dictionary. The key names are the same as field names in opt.opt
    class.

    Returns:
        Options as a python dictionary storing parameters as numpy arrays
    """

    # This ugly thing converts class variables into dictionary
    f = {key:value for key, value in opt.__dict__.items() if
         not key.startswith('__') and not callable(key)}

    f = settypes(f)
    return f

def settypes(f):
    """
    Convert option parameters in dictionary to correct format

    This is a helper routine for generateopt() but can also be called if the
    option parameters in dictionary were modified after the options dictionary
    was produced by generateopt() (if in case a plain number was assigned
    instead of proper numpy array). This function makes sure all parameters are
    numpy arrays of correct type.

    Args:
        f: Option dictionary
    Returns:
        Same dictionary but with all parameters ensured to be in correct format
    """

    for i in f:
        if type(f[i]) is not np.array:
            f[i] = np.array(f[i])

        f[i] =f[i].astype("f8")

    return f
