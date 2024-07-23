/**
 * @file rfof.h
 * @brief Contains the functions to be called from the simulation loop when
 * using ICRH.
**/
#include <stdlib.h>
#include "ascot5.h"

#ifndef RFOF_H
#define RFOF_H

/**
 * @brief Reusable struct for storing marker specific data during the simulation
 * loop.
 *
 * The data in this struct is altered during the simulation. Only pointers are
 * stored as the actual data is stored on the Fortran side.
 */
typedef struct rfof_marker {
    void* p[NSIMD]; /**< The marker struct in a format required by librfof */
    void* history_array[NSIMD]; /**< Stores values of the resonance function for
                                     estimating the next time-step */
    void* diag_array[NSIMD]; /**< C equivalents of Fortran diagnostics pointers
                                  which are required but unused at the moment */
    int nrow[NSIMD]; /**< Number of rows in an resonance history matrix    */
    int ncol[NSIMD]; /**< Number of columns in an resonance history matrix */

} rfof_marker;

/** @brief RFOF simulation input data
 *
 * Immutable input data shared between all markers. The actual data is stored
 * in the Fortran side and this struct only stores the pointers.
*/
typedef struct {
    void* rfof_input_params; /**< Pointer to rfof_input_param struct on
                                  the fortran side                 */
    void* rfglobal;          /**< Wave field; same for all markers */
} rfof_data;

void rfof_init_offload(rfof_data* rfof_data);

void rfof_init(rfof_data* rfof, rfof_data* rfof_offload_data);

void rfof_free_offload(rfof_data* rfof);

void rfof_set_marker_manually(void** cptr_marker, int* id,
    real* weight, real* R, real* phi, real* z, real* psi, real* charge,
    real* mass, real* Ekin, real* vnorm, real* mu, real* Pphi,
    real* vpar, real* vperp, real* gyrof, real* vdriftRho, real* acceleration,
    int* is_accelerated, int* is_already_allocated);

void rfof_set_up(rfof_marker* rfof_mrk, rfof_data* rfof_data);

void rfof_tear_down(rfof_marker* rfof_mrk);

void rfof_clear_history(rfof_marker* rfof_mrk, int imrk);

void rfof_resonance_check_and_kick_gc(
    particle_simd_gc* p, real* hin, real* hout_rfof, rfof_marker* rfof_mrk,
    rfof_data* rfof_data, B_field_data* Bdata);

void rfof_eval_rf_wave(
    real* e_plus_real, real* e_minus_real, real* e_plus_imag,
    real* e_minus_imag, real R, real z, rfof_data* rfof);
void rfof_eval_resonance_function(
    real* omega_res, int* nharm, rfof_marker* rfof_mrk, rfof_data* rfof);
#endif
