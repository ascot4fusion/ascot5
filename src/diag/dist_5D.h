/**
 * @file dist_5D.h
 * @brief Header file for dist_5D.c
 */
#ifndef DIST_5D_H
#define DIST_5D_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"

/**
 * @brief Histogram parameters that will be offloaded to target
 */
typedef struct {
    int n_r;          /**< number of r bins                     */
    real min_r;       /**< value of lowest r bin                */
    real max_r;       /**< value of highest r bin               */

    int n_phi;        /**< number of phi bins                   */
    real min_phi;     /**< value of lowest phi bin              */
    real max_phi;     /**< value of highest phi bin             */

    int n_z;          /**< number of z bins                     */
    real min_z;       /**< value of lowest z bin                */
    real max_z;       /**< value of highest z bin               */

    int n_ppara;      /**< number of p_parallel bins            */
    real min_ppara;   /**< value of lowest p_parallel bin       */
    real max_ppara;   /**< value of highest p_parallel bin      */

    int n_pperp;      /**< number of p_perpendicular bins       */
    real min_pperp;   /**< value of lowest p_perpendicular bin  */
    real max_pperp;   /**< value of highest p_perpendicular bin */

    int n_time;       /**< number of time bins                  */
    real min_time;    /**< value of lowest time bin             */
    real max_time;    /**< value of highest time bin            */

    int n_q;          /**< number of charge bins                */
    real min_q;       /**< value of lowest charge bin           */
    real max_q;       /**< value of highest charge bin          */
} dist_5D_offload_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct {
    int n_r;          /**< number of r bins                     */
    real min_r;       /**< value of lowest r bin                */
    real max_r;       /**< value of highest r bin               */

    int n_phi;        /**< number of r bins                     */
    real min_phi;     /**< value of lowest r bin                */
    real max_phi;     /**< value of highest r bin               */

    int n_z;          /**< number of z bins                     */
    real min_z;       /**< value of lowest z bin                */
    real max_z;       /**< value of highest z bin               */

    int n_ppara;      /**< number of p_parallel bins            */
    real min_ppara;   /**< value of lowest p_parallel bin       */
    real max_ppara;   /**< value of highest p_parallel bin      */

    int n_pperp;      /**< number of p_perpendicular bins       */
    real min_pperp;   /**< value of lowest p_perpendicular bin  */
    real max_pperp;   /**< value of highest p_perpendicular bin */

    int n_time;       /**< number of r bins                     */
    real min_time;    /**< value of lowest r bin                */
    real max_time;    /**< value of highest r bin               */

    int n_q;          /**< number of r bins                     */
    real min_q;       /**< value of lowest r bin                */
    real max_q;       /**< value of highest r bin               */

    size_t step_1;    /**< step for 2nd fastest running index   */
    size_t step_2;    /**< step for 3rd fastest running index   */
    size_t step_3;    /**< step for 4th fastest running index   */
    size_t step_4;    /**< step for 5th fastest running index   */
    size_t step_5;    /**< step for 6th fastest running index   */
    size_t step_6;    /**< step for 7th fastest running index   */

    real* histogram;  /**< pointer to start of histogram array */
} dist_5D_data;

#pragma omp declare target
size_t dist_5D_index(int i_r, int i_phi, int i_z, int i_ppara, int i_pperp,
                     int i_time, int i_q, size_t step_6, size_t step_5,
                     size_t step_4, size_t step_3, size_t step_2,
                     size_t step_1);
void dist_5D_init(dist_5D_data* dist_data,
                  dist_5D_offload_data* offload_data,
                  real* offload_array);
void dist_5D_update_fo(dist_5D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i);
void dist_5D_update_gc(dist_5D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i);
#pragma omp end declare target

#endif
