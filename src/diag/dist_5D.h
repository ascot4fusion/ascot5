/**
 * @file dist_5D.h
 * @brief Header file for dist_5D.c
 */
#ifndef DIST_5D_H
#define DIST_5D_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"
#include "../options.h"

/**
 * @brief Histogram parameters
 */
typedef struct {
    size_t step_1;    /**< step for 2nd fastest running index   */
    size_t step_2;    /**< step for 3rd fastest running index   */
    size_t step_3;    /**< step for 4th fastest running index   */
    size_t step_4;    /**< step for 5th fastest running index   */
    size_t step_5;    /**< step for 6th fastest running index   */
    size_t step_6;    /**< step for 7th fastest running index   */

    real* histogram;  /**< pointer to start of histogram array */
} dist_5D_data;

size_t dist_5D_index(int i_r, int i_phi, int i_z, int i_ppara, int i_pperp,
                     int i_time, int i_q, size_t step_6, size_t step_5,
                     size_t step_4, size_t step_3, size_t step_2,
                     size_t step_1);
int dist_5D_init(dist_5D_data* data, sim_parameters* params);
void dist_5D_free(dist_5D_data* data);
void dist_5D_offload(dist_5D_data* data);
void dist_5D_update_fo(dist_5D_data* dist, sim_parameters* params,
                       particle_simd_fo* p_f, particle_simd_fo* p_i);
void dist_5D_update_gc(dist_5D_data* dist, sim_parameters* params,
                       particle_simd_gc* p_f, particle_simd_gc* p_i);

#endif
