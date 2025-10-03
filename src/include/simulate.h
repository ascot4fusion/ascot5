/**
 * @file simulate.h
 * Contains declarations of simulation_offload_data and simulation_data structs.
 * Also simulation mode enums are declared here.
 */
#ifndef SIMULATE_H
#define SIMULATE_H

#include "B_field.h"
#include "E_field.h"
#include "ascot5.h"
#include "asigma.h"
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
    B_field_data B_data;       /**< Magnetic field interface.                 */
    E_field_data E_data;       /**< Electric field interface.                 */
    plasma_data plasma_data;   /**< Plasma data interface.                    */
    neutral_data neutral_data; /**< Neutral data interface.                   */
    wall_data wall_data;       /**< Wall data interface.                      */
    boozer_data *boozer_data;  /**< Boozer data interface.                    */
    mhd_data mhd_data;         /**< MHD data interface.                       */
    asigma_data asigma_data;   /**< Atomic sigma data interface.              */
    nbi_data nbi_data;         /**< Neutral beam injection data interface.    */
    rfof_data rfof_data;       /**< Void pointers to ICRH wave field and input
                                    parameters Fortran structs.               */
    Diagnostics *diagnostics;  /**< Diagnostics interface.                    */
    random_data *random_data;  /**< Random number generator.                  */
    mccc_data *mccc_data;      /**< Tabulated special functions and collision
                                    operator parameters.                      */
    sim_parameters *params;    /**< Simulation parameters.                    */

} sim_data;

void simulate(int nmarkers, particle_state *p, sim_data *sim);

#endif
