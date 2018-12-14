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
#include "diag.h"
#include "offload.h"
#include "random.h"

/**
 * @brief Simulaton modes
 *
 * These enums are used to determine which simulation mode will be executed.
 *
 * @todo Make these real enums.
 */
enum {
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
 * @brief Simulation offload struct
 *
 * This structure holds necessary data to initialize the simulation data struct
 * target. Any IO related data (input filenames etc.) are also stored here (but
 * not in the simulation data struct as these are not needed on target).
 */
typedef struct {
    /* Input and diagnostic interface offload data */
    B_field_offload_data B_offload_data;       /**< Magnetic field offload data */
    E_field_offload_data E_offload_data;       /**< Electric field offload data */
    plasma_offload_data plasma_offload_data;   /**< Plasma offload data         */
    neutral_offload_data neutral_offload_data; /**< Neutral offload data        */
    wall_offload_data wall_offload_data;       /**< Wall offload data           */
    diag_offload_data diag_offload_data;       /**< Diagnostics offload data    */

    /* Options - general */
    int sim_mode;        /**< Which simulation mode is used                   */
    int enable_ada;      /**< Is adaptive time-step used                      */
    int record_GOasGC;   /**< Is particles GC coordinates used in diagnostics */

    /* Options - fixed time-step */
    int fix_usrdef_use;  /**< Use user defined value for (initial) time-step  */
    real fix_usrdef_val; /**< User defined time-step value                    */
    int fix_stepsPerGO;  /**< Time-step = gyrotime/fix_stepsPerGO if not
                              explicitly user defined                         */

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
    int disable_gctransform;   /**< Disables first order velocity terms in
                                    guiding center transformation             */
    int disable_energyccoll;   /**< Disables energy component from Coulomb
                                    collisions */
    int disable_pitchccoll;    /**< Disables pitch component from Coulomb
                                    collisions */
    int disable_gcdiffccoll;   /**< Disables guiding center spatial diffusion
                                    from Coulomb collisions */

    /* Options - end conditions */
    int endcond_active;        /**< Bit array notating active end conditions  */
    real endcond_maxSimTime;   /**< Maximum simulation time [s]               */
    real endcond_maxCpuTime;   /**< Maximum wall-clock time [s]               */
    real endcond_minRho;       /**< Minimum rho limit                         */
    real endcond_maxRho;       /**< Maximum rho limit                         */
    real endcond_minEkin;      /**< Fixed minimum kinetic energy limit [J]    */
    real endcond_minEkinPerTe; /**< Thermal minimum energy limit is this
                                    parameter times local thermal energy      */
    real endcond_maxTorOrb;    /**< Maximum limit for toroidal distance [rad] */
    real endcond_maxPolOrb;    /**< Maximum limit for poloidal distance [rad] */

    /* Metadata */
    char hdf5_in[256];
    char hdf5_out[256];
    char outfn[256];
    char qid[256];
    char description[256];

    int mpi_rank;
    int mpi_size;

} sim_offload_data;

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
    diag_data diag_data;       /**< Diagnostics data interface                */

    /* Options - general */
    int sim_mode;        /**< Which simulation mode is used                   */
    int enable_ada;      /**< Is adaptive time-step used                      */
    int record_GOasGC;   /**< Is particles GC coordinates used in diagnostics */

    /* Options - fixed time-step */
    int fix_usrdef_use;  /**< Use user defined value for (initial) time-step  */
    real fix_usrdef_val; /**< User defined time-step value                    */
    int fix_stepsPerGO;  /**< Time-step = gyrotime/fix_stepsPerGO if not
                              explicitly user defined                         */

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
    int disable_gctransform;   /**< Disables first order velocity terms in
                                    guiding center transformation             */
    int disable_energyccoll;   /**< Disables energy component from Coulomb
                                    collisions */
    int disable_pitchccoll;    /**< Disables pitch component from Coulomb
                                    collisions */
    int disable_gcdiffccoll;   /**< Disables guiding center spatial diffusion
                                    from Coulomb collisions */

    /* Options - end conditions */
    int endcond_active;        /**< Bit array notating active end conditions  */
    real endcond_maxSimTime;   /**< Maximum simulation time [s]               */
    real endcond_maxCpuTime;   /**< Maximum wall-clock time [s]               */
    real endcond_minRho;       /**< Minimum rho limit                         */
    real endcond_maxRho;       /**< Maximum rho limit                         */
    real endcond_minEkin;      /**< Fixed minimum kinetic energy limit [J]    */
    real endcond_minEkinPerTe; /**< Thermal minimum energy limit is this
                                    parameter times local thermal energy      */
    real endcond_maxTorOrb;    /**< Maximum limit for toroidal distance [rad] */
    real endcond_maxPolOrb;    /**< Maximum limit for poloidal distance [rad] */

    /* Metadata */
    random_data random_data;   /**< Random number generator                   */
    real* coldata;             /**< Look-up tables for collision operator if
                                    collision coefficients are interpolated
                                    and not calculated run-time               */
} sim_data;

#pragma omp declare target
void simulate(int id, int n_particles, particle_state* p,
              sim_offload_data* sim_offload,
              offload_package* offload_data,
              real* offload_array, real* diag_offload_array);
#pragma omp end declare target

#endif
