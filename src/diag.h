/**
 * @file diag.h
 * @brief Header file for diag.c
 */
#ifndef DIAG_H
#define DIAG_H
#include "ascot5.h"
#include "particle.h"
#include "B_field.h"
#include "diag/dist_5D.h"
#include "diag/dist_6D.h"
#include "diag/dist_rho5D.h"
#include "diag/dist_rho6D.h"
#include "diag/dist_com.h"
#include "diag/diag_orb.h"
#include "diag/diag_transcoef.h"

/**
 * @brief Diagnostics data struct
 */
typedef struct {
    int diagorb_collect;   /**< Flag for collecting orbit data               */
    int dist5D_collect;    /**< Flag for collecting 5D distribution          */
    int dist6D_collect;    /**< Flag for collecting 6D distribution          */
    int distrho5D_collect; /**< Flag for collecting 5D rho distribution      */
    int distrho6D_collect; /**< Flag for collecting 6D rho distribution      */
    int distCOM_collect;   /**< Flag for collecting COM distribution         */
    int diagtrcof_collect; /**< Flag for collecting transport coefficients   */

    diag_orb_data diagorb;     /**< Orbit diagnostics data                   */
    dist_5D_data dist5D;       /**< 5D distribution diagnostics data         */
    dist_6D_data dist6D;       /**< 6D distribution diagnostics data         */
    dist_rho5D_data distrho5D; /**< 5D rho distribution diagnosticsd data    */
    dist_rho6D_data distrho6D; /**< 6D rho distribution diagnostics data     */
    dist_COM_data distCOM;     /**< COM distribution diagnostics data        */
    diag_transcoef_data diagtrcof; /**< Transp. Coef. diagnostics data       */

} diag_data;

int diag_init(diag_data* data, int Nmrk);
void diag_free(diag_data* data);
void diag_offload(diag_data* data);

void diag_sum(diag_data* data1, diag_data* data2);

void diag_update_fo(diag_data* data, B_field_data* Bdata, particle_simd_fo* p_f,
                    particle_simd_fo* p_i, int n_running_ref);

void diag_update_gc(diag_data* data, B_field_data* Bdata, particle_simd_gc* p_f,
                    particle_simd_gc* p_i);

void diag_update_ml(diag_data* data, particle_simd_ml* p_f,
                    particle_simd_ml* p_i);


#endif
