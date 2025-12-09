/**
 * @file simulate.h
 * @brief Header file for simulate.c
 *
 * Contains declarations of simulation_offload_data and simulation_data structs.
 * Also simulation mode enums are declared here.
 */
#ifndef SIMULATE_H
#define SIMULATE_H

#include "ascot5.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "boozer.h"
#include "mhd.h"
#include "asigma.h"
#include "nbi.h"
#include "diag.h"
#include "random.h"
#include "simulate/mccc/mccc.h"
#include "rfof.h"

/**
 * @brief Simulaton modes
 *
 * These enums are used to determine which simulation mode will be executed.
 */
enum SIMULATION_MODE {
    /** Models markers as particles using particle_simd_fo struct and
        simulate_fo_fixed.c simulation loop                                 */
    simulate_mode_fo = 1,
    /** Models markers as guiding centers using particle_simd_gc struct and
        simulate_gc_fixed.c or simulate_gc_adaptive.c simulation loops      */
    simulate_mode_gc = 2,
    /** Models markers first like using simulate_mode_gc. Additional end
        condition is used for markers that get close to wall. After all
        markers are finished, simulation for markers that were close to the
        wall is continued with using simulate_mode_fo mode                  */
    simulate_mode_hybrid = 3,
    /** Models markers as field lines using particle_simd_ml struct and
        simulate_ml_adaptive.c simulation loop                              */
    simulate_mode_ml = 4
};

/**
 * @brief Simulation data struct
 *
 * This structure holds all data required to simulate markers except the
 * markers themselves. Options, initialized input data, etc. are all stored
 * here.
 *
 * Even though most fields are identical to sim_offload_data, a separate struct
 * is required as pointers cannot be offloaded.
 */
typedef struct {
    /* Input and diagnostic interfaces */
    B_field_data B_data;       /**< Magnetic field interface                  */
    E_field_data E_data;       /**< Electric field interface                  */
    plasma_data plasma_data;   /**< Plasma data interface                     */
    neutral_data neutral_data; /**< Neutral data interface                    */
    wall_data wall_data;       /**< Wall data interface                       */
    boozer_data boozer_data;   /**< Boozer data interface                     */
    mhd_data mhd_data;         /**< MHD data interface                        */
    asigma_data asigma_data;   /**< Atomic sigma data interface               */
    nbi_data nbi_data;         /**< Neutral beam injection data interface     */
    diag_data diag_data;       /**< Diagnostics data interface                */
    rfof_data rfof_data;       /**< Void pointers to ICRH wave field and input
                                    parameters Fortran structs.               */

    /* Metadata */
    random_data random_data;   /**< Random number generator                   */
    mccc_data mccc_data;       /**< Tabulated special functions and collision
                                    operator parameters                       */
    /* Options - general */
    int sim_mode;        /**< Which simulation mode is used                   */
    int enable_ada;      /**< Is adaptive time-step used                      */
    int record_mode;     /**< Which record mode is used                       */

    /* Options - fixed time-step */
    int fix_usrdef_use;    /**< Use user defined value for (initial) time-step*/
    real fix_usrdef_val;   /**< User defined time-step value                  */
    int fix_gyrodef_nstep; /**< Time-step = gyrotime/fix_stepsPerGO if not
                                explicitly user defined                       */

    /* Options - adaptive time-step */
    real ada_tol_orbfol;       /**< Tolerance for relative error in
                                    orbit-following                           */
    real ada_tol_clmbcol;      /**< Tolerance for relative error in Coulomb
                                    collisions                                */
    real ada_max_drho;         /**< Maximum rho distance marker is allowed to
                                    travel during single adaptive time-step   */
    real ada_max_dphi;         /**< Maximum phi distance marker is allowed to
                                    travel during single adaptive time-step   */

    /* Options - physics */
    int enable_orbfol;         /**< Is orbit-following enabled                */
    int enable_clmbcol;        /**< Are Coulomb collisions enabled            */
    int enable_mhd;            /**< Are MHD modes enabled                     */
    int enable_atomic;         /**< Are atomic reactions enabled              */
    int enable_icrh;           /**< Is RFOF enabled                           */
    int enable_aldforce;       /**< Is radiation reaction force enabled       */
    int enable_flr_losses;     /**< Are finite Larmor radius losses enabled   */
    int disable_gctransform;   /**< Disables first order velocity terms in
                                    guiding center transformation             */
    int disable_energyccoll;   /**< Disables energy component from Coulomb
                                    collisions */
    int disable_pitchccoll;    /**< Disables pitch component from Coulomb
                                    collisions */
    int disable_gcdiffccoll;   /**< Disables guiding center spatial diffusion
                                    from Coulomb collisions */
    int reverse_time;          /**< Set time running backwards in simulation  */

    /* Options - end conditions */
    int endcond_active;        /**< Bit array notating active end conditions  */
    real endcond_lim_simtime;  /**< Simulation time limit [s]                 */
    real endcond_max_mileage;  /**< Maximum simulation duration [s]           */
    real endcond_max_cputime;  /**< Maximum wall-clock time [s]               */
    real endcond_min_rho;      /**< Minimum rho limit                         */
    real endcond_max_rho;      /**< Maximum rho limit                         */
    real endcond_min_ekin;     /**< Fixed minimum kinetic energy limit [J]    */
    real endcond_min_thermal;  /**< Thermal minimum energy limit is this
                                    parameter times local thermal energy      */
    real endcond_max_tororb;   /**< Maximum limit for toroidal distance [rad] */
    real endcond_max_polorb;   /**< Maximum limit for poloidal distance [rad] */
    int endcond_torandpol;     /**< Flag whether both tor and pol must be met */

    /* Metadata */
    char hdf5_in[256];     /**< Name of the input HDF5 file  */
    char hdf5_out[256];    /**< Name of the output HDF5 file */
    char qid[256];         /**< QID of current run           */
    char description[256]; /**< Current run's description    */

    int mpi_root; /**< Rank of the root process      */
    int mpi_rank; /**< Rank of this MPI process      */
    int mpi_size; /**< Total number of MPI processes */

    /* QIDs for inputs if the active inputs are not used */
    char qid_options[256]; /**< Options QID if active not used */
    char qid_bfield[256];  /**< Bfield QID if active not used  */
    char qid_efield[256];  /**< Efield QID if active not used  */
    char qid_marker[256];  /**< Marker QID if active not used  */
    char qid_wall[256];    /**< Wall QID if active not used    */
    char qid_plasma[256];  /**< Plasma QID if active not used  */
    char qid_neutral[256]; /**< Neutral QID if active not used */
    char qid_boozer[256];  /**< Boozer QID if active not used  */
    char qid_mhd[256];     /**< MHD QID if active not used     */
    char qid_asigma[256];  /**< Asigma QID if active not used  */
    char qid_nbi[256];     /**< NBI QID if active not used     */

} sim_data;

void simulate_init(sim_data* sim);

void simulate(int n_particles, particle_state* p, sim_data* sim);

#endif
