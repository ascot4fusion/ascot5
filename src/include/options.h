/**
 * @file options.h
 * Simulation options.
 */
#ifndef OPTIONS_H
#define OPTIONS_H

#include "defines.h"

#define MAXPOINCARE 32

/**
 * Simulation modes.
 *
 * These enums are used to determine which simulation mode will be executed.
 */
enum SIMULATION_MODE
{
    /** Models markers as particles using MarkerGyroOrbit struct and
        simulate_fo_fixed.c simulation loop .                                  */
    simulate_mode_fo = 1,
    /** Models markers as guiding centers using MarkerGuidingCenter struct and
        simulate_gc_fixed.c or simulate_gc_adaptive.c simulation loops.      */
    simulate_mode_gc = 2,
    /** Models markers first like using simulate_mode_gc. Additional end
        condition is used for markers that get close to wall. After all
        markers are finished, simulation for markers that were close to the
        wall is continued with using simulate_mode_fo mode.                 */
    simulate_mode_hybrid = 3,
    /** Models markers as field lines using MarkerFieldLine struct and
        simulate_ml_adaptive.c simulation loop.                              */
    simulate_mode_ml = 4
};

/**
 * Simulation parameters.
 */
typedef struct
{
    /**
     * Indicates if we are tracing gyro orbits, guiding centers, or field lines
     * or using the hybrid mode.
     */
    int simulation_mode;

    /** Is adaptive time-step used. */
    int enable_adaptive;

    /** Which record mode is used. */
    int record_mode;

    /** User defined value for the (initial) time-step. */
    int use_explicit_fixedstep;

    /** Time-step is gyrotime divided by this value. */
    int gyrodefined_fixedstep;

    /** User defined time-step value. */
    real explicit_fixedstep;

    /** Tolerance for relative error in orbit-following. */
    real adaptive_tolerance_orbit;

    /** Tolerance for relative error in Coulomb collisions. */
    real adaptive_tolerance_collisions;

    /**
     * Maximum rho distance marker is allowed to travel during single adaptive
     * time-step.
     */
    real adaptive_max_drho;

    /**
     * Maximum phi distance marker is allowed to travel during single adaptive
     * time-step.
     */
    real adaptive_max_dphi;

    /** Toggle orbit-following. */
    int enable_orbit_following;

    /** Toggle Coulomb collisions. */
    int enable_coulomb_collisions;

    /** Toggle MHD. */
    int enable_mhd;

    /** Toggle atomic reactions. */
    int enable_atomic;

    /** Toggle RFOF. */
    int enable_icrh;

    /** Toggle radiation losses. */
    int enable_aldforce;

    /** Disables first order velocity terms in guiding center transformation. */
    int disable_first_order_gctransformation;

    /** Toggle GC energy operator. */
    int disable_ccoll_gcenergy;

    /** Toggle GC pitch operator. */
    int disable_ccoll_gcpitch;

    /** Toggle GC spatial operator. */
    int disable_ccoll_gcspatial;

    /** Set time running backwards in simulation. */
    int reverse_time;

    /** Bit array indicating active end conditions. */
    int endcond_active;

    /** Both toroidal and poloidal limits has to be met if this is true. */
    int require_both_tor_and_pol;

    /** Simulation time limit [s]. */
    real lab_time_limit;

    /** Maximum simulation duration [s]. */
    real max_mileage;

    /** Maximum wall-clock time [s]. */
    real max_real_time;

    /** Minimum and maximum rho limit. */
    real rho_coordinate_limits[2];

    /** Fixed minimum kinetic energy limit [J]. */
    real min_energy;

    /** Minimum energy limit is this times local thermal energy. */
    real min_local_thermal_energy;

    /** Maximum limit for toroidal distance [rad]. */
    real max_number_of_toroidal_orbits;

    /** Maximum limit for poloidal distance [rad]. */
    real max_number_of_poloidal_orbits;

    /** Flag for collecting orbit data. */
    int collect_orbit;

    /** Collect 5D distribution. */
    int collect_dist5d;

    /** Collect 6D distribution. */
    int collect_dist6d;

    /** Collect 5D rho distribution. */
    int collect_dist5drho;

    /** Collect 6D rho distribution. */
    int collect_dist6drho;

    /** Collect COM distribution. */
    int collect_distcom;

    /** Collect transport coefficients. */
    int collect_transport_coefficient;

    /** Record orbit every crossing. */
    int poincare;

    /** Number of points in record. */
    size_t number_of_points_per_marker;

    /** Number of toroidal Poincaré planes. */
    size_t ntoroidalplots;

    /** Number of poloidal Poincaré planes. */
    size_t npoloidalplots;

    /** Number of radial Poincaré surfaces. */
    size_t nradialplots;

    /** Interval for recording orbit. */
    real interval;

    /** Toroidal plane angles. */
    real toroidal_angles[MAXPOINCARE];

    /** Poloidal plane angles. */
    real poloidal_angles[MAXPOINCARE];

    /** Radial plane angles. */
    real radial_distances[MAXPOINCARE];

    /**
     * How many data points are averaged before calculating a transport
     * coefficient.
     */
    size_t number_of_points_to_average;

    /** Evaluate transport coefficients as a function of rho or R. */
    int record_rho_instead_of_r;

    /**
     * Margin before making a new record in transport coefficient evaluation.
     */
    real margin;

    /** Number of bins in R abscissa. */
    int r_bins;

    /** Number of bins in phi abscissa. */
    int phi_bins;

    /** Number of bins in z abscissa. */
    int z_bins;

    /** Number of bins in ppara abscissa. */
    int ppara_bins;

    /** Number of bins in pperp abscissa. */
    int pperp_bins;

    /** Number of bins in time abscissa. */
    int time_bins;

    /** Number of bins in charge abscissa. */
    int charge_bins;

    /** Min/max values in R abscissa. */
    real r_interval[2];

    /** Min/max values in phi abscissa. */
    real phi_interval[2];

    /** Min/max values in z abscissa. */
    real z_interval[2];

    /** Min/max values in ppara abscissa. */
    real ppara_interval[2];

    /** Min/max values in pperp abscissa. */
    real pperp_interval[2];

    /** Min/max values in time abscissa. */
    real time_interval[2];

    /** Min/max values in charge abscissa. */
    real charge_interval[2];

} Options;

#endif
