/**
 * @file simulate.h
 * Contains declarations of simulation_offload_data and simulation_data structs.
 * Also simulation mode enums are declared here.
 */
#ifndef SIMULATE_H
#define SIMULATE_H

#include "bfield.h"
#include "efield.h"
#include "defines.h"
#include "atomic.h"
#include "boozer.h"
#include "diag.h"
#include "mhd.h"
#include "nbi.h"
#include "neutral.h"
#include "options.h"
#include "plasma.h"
#include "random.h"
#include "rfof.h"
#include "wall.h"
#include "diag.h"

/**
 * Simulation data struct.
 *
 * This structure holds all data required to simulate markers except the
 * markers themselves. Options, initialized input data, etc. are all stored
 * here.
 *
 * Even though most fields are identical to sim_offload_data, a separate struct
 * is required as pointers cannot be offloaded.
 */
typedef struct
{
    Bfield bfield;       /**< Magnetic field interface.                 */
    Efield efield;       /**< Electric field interface.                 */
    Plasma plasma;   /**< Plasma data interface.                    */
    Neutral neutral; /**< Neutral data interface.                   */
    Wall wall;       /**< Wall data interface.                      */
    Boozer *boozer;  /**< Boozer data interface.                    */
    Mhd mhd;         /**< MHD data interface.                       */
    Atomic atomic;   /**< Atomic sigma data interface.              */
    Nbi nbi;         /**< Neutral beam injection data interface.    */
    Rfof rfof;       /**< Void pointers to ICRH wave field and input
                                    parameters Fortran structs.               */
    Diagnostics *diagnostics;  /**< Diagnostics interface.                    */
    random_data *random_data;  /**< Random number generator.                  */
    mccc_data *mccc_data;      /**< Tabulated special functions and collision
                                    operator parameters.                      */
    sim_parameters *params;    /**< Simulation parameters.                    */

} sim_data;

void simulate(int nmarkers, particle_state *p, sim_data *sim);

/**
 * Simulates particles using fixed time-step
 *
 * The simulation includes:
 * - orbit-following with Volume-Preserving Algorithm
 * - Coulomb collisions with Euler-Maruyama method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined: either a directly given fixed value
 * or a given fraction of gyrotime.
 *
 * @param pq particles to be simulated
 * @param sim simulation data struct
 * @param mrk_array_size size of particle arrays
 */
void simulate_fo_fixed(particle_queue* pq, sim_data* sim, int mrk_array_size);

/**
 * Simulates guiding centers using adaptive time-step
 *
 * The simulation is carried until all marker have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
 *
 * @param pq particles to be simulated
 * @param sim simulation data
 */
void simulate_gc_adaptive(particle_queue* pq, sim_data* sim);

/**
 * Simulates guiding centers using fixed time-step
 *
 * The simulation includes:
 * - orbit-following with RK4 method
 * - Coulomb collisions with Euler-Maruyama method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined either directly or as a fraction of
 * gyrotime.
 *
 * @param pq particles to be simulated
 * @param sim simulation data
 */
void simulate_gc_fixed(particle_queue* pq, sim_data* sim);

/**
 * Simulates magnetic field-lines using adaptive time-step.
 *
 * The simulation includes:
 * - orbit-following with Cash-Karp method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
 *
 * Note that even though we might refer the integration time-step
 * as "time", in reality we are integrating over distance. The time
 * step is therefore step in meters marker orbit is intgerated. Marker does have
 * time in its field, but it is the global time and that is not being changed
 * during the simulation.
 *
 * @param pq field lines to be simulated
 * @param sim simulation data struct
 */
void simulate_fl_adaptive(particle_queue* pq, sim_data* sim);

#endif
