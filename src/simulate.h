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
#include "options.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "boozer.h"
#include "mhd.h"
#include "asigma.h"
#include "nbi.h"
#include "random.h"
#include "simulate/mccc/mccc.h"
#include "rfof.h"
#include "diag/dist_5D.h"
#include "diag/dist_6D.h"
#include "diag/dist_rho5D.h"
#include "diag/dist_rho6D.h"
#include "diag/dist_com.h"
#include "diag/diag_orb.h"
#include "diag/diag_transcoef.h"

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
    /* Input */
    B_field_data B_data;       /**< Magnetic field interface                  */
    E_field_data E_data;       /**< Electric field interface                  */
    plasma_data plasma_data;   /**< Plasma data interface                     */
    neutral_data neutral_data; /**< Neutral data interface                    */
    wall_data wall_data;       /**< Wall data interface                       */
    boozer_data* boozer_data;  /**< Boozer data interface                     */
    mhd_data mhd_data;         /**< MHD data interface                        */
    asigma_data asigma_data;   /**< Atomic sigma data interface               */
    nbi_data nbi_data;         /**< Neutral beam injection data interface     */
    rfof_data rfof_data;       /**< Void pointers to ICRH wave field and input
                                    parameters Fortran structs.               */

    /* Diagnostics */
    diag_orb_data* orbit;     /**< Orbit diagnostics data                   */
    dist_5D_data* dist5d;       /**< 5D distribution diagnostics data         */
    dist_6D_data* dist6d;       /**< 6D distribution diagnostics data         */
    dist_rho5D_data* dist5drho; /**< 5D rho distribution diagnosticsd data    */
    dist_rho6D_data* dist6drho; /**< 6D rho distribution diagnostics data     */
    dist_COM_data* distcom;     /**< COM distribution diagnostics data        */
    diag_transcoef_data* transport_coefficient; /**< Transp. Coef. diagnostics data       */

    /* Metadata */
    random_data* random_data;  /**< Random number generator                   */
    mccc_data* mccc_data;      /**< Tabulated special functions and collision
                                    operator parameters                       */
    sim_parameters* params;    /**< Simulation parameters                     */

} sim_data;

void simulate_init(sim_data* sim, int nmarkers);

void simulate(int nmarkers, particle_state* p, sim_data* sim);

#endif
