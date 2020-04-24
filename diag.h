/**
 * @file diag.h
 * @brief Header file for diag.c
 */
#ifndef DIAG_H
#define DIAG_H
#include "ascot5.h"
#include "particle.h"
#include "diag/dist_5D.h"
#include "diag/dist_6D.h"
#include "diag/dist_rho5D.h"
#include "diag/dist_rho6D.h"
#include "diag/diag_orb.h"

/**
 * @brief Diagnostics offload data struct
 */
typedef struct {
    int diagorb_collect;   /**< Flag for collecting orbit data               */
    int dist5D_collect;    /**< Flag for collecting 5D distribution          */
    int dist6D_collect;    /**< Flag for collecting 6D distribution          */
    int distrho5D_collect; /**< Flag for collecting 5D rho distribution      */
    int distrho6D_collect; /**< Flag for collecting 6D rho distribution      */

    diag_orb_offload_data diagorb;     /**< Orbit offload data               */
    dist_5D_offload_data dist5D;       /**< 5D distribution offload data     */
    dist_6D_offload_data dist6D;       /**< 6D distribution offload data     */
    dist_rho5D_offload_data distrho5D; /**< 5D rho distribution offload data */
    dist_rho6D_offload_data distrho6D; /**< 6D rho distribution offload data */

    int offload_dist5D_index;    /**< Index for 5D dist in offload array     */
    int offload_dist6D_index;    /**< Index for 5D dist in offload array     */
    int offload_distrho5D_index; /**< Index for 5D dist in offload array     */
    int offload_distrho6D_index; /**< Index for 6D rho dist in offload array */
    int offload_diagorb_index;   /**< Index for orbit data in offload array  */
    int offload_array_length;    /**< Number of elements in offload_array    */

} diag_offload_data;

/**
 * @brief Diagnostics data struct
 */
typedef struct {
    int diagorb_collect;   /**< Flag for collecting orbit data               */
    int dist5D_collect;    /**< Flag for collecting 5D distribution          */
    int dist6D_collect;    /**< Flag for collecting 6D distribution          */
    int distrho5D_collect; /**< Flag for collecting 5D rho distribution      */
    int distrho6D_collect; /**< Flag for collecting 6D rho distribution      */

    diag_orb_data diagorb;     /**< Orbit diagnostics data                   */
    dist_5D_data dist5D;       /**< 5D distribution diagnostics data         */
    dist_6D_data dist6D;       /**< 6D distribution diagnostics data         */
    dist_rho5D_data distrho5D; /**< 5D rho distribution diagnosticsd data    */
    dist_rho6D_data distrho6D; /**< 6D rho distribution diagnostics data     */

} diag_data;

int diag_init_offload(diag_offload_data* data, real** offload_array, int Nmrk);

void diag_free_offload(diag_offload_data* data, real** offload_array);

void diag_sum(diag_offload_data* data, real* array1, real* array2);

#pragma omp declare target
void diag_init(diag_data* data, diag_offload_data* offload_data,
               real* offload_array);

void diag_free(diag_data* data);

void diag_update_fo(diag_data* data, particle_simd_fo* p_f,
                    particle_simd_fo* p_i);

void diag_update_gc(diag_data* data, particle_simd_gc* p_f,
                    particle_simd_gc* p_i);

void diag_update_ml(diag_data* data, particle_simd_ml* p_f,
                    particle_simd_ml* p_i);

#pragma omp end declare target

#endif
