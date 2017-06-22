/**
 * @file distributions.h
 * @brief Header file for distributions.c
 */
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "ascot5.h"
#include "particle.h"

/**
 * @brief Histogram parameters that will be offloaded to target
 */
typedef struct {
    int n_r;          /**< number of r bins */ 
    real min_r;       /**< value of lowest r bin */
    real max_r;       /**< value of highest r bin */
    int n_z;          /**< number of z bins */
    real min_z;       /**< value of lowest z bin */
    real max_z;       /**< value of highest z bin */
    int n_vpara;      /**< number of v_parallel bins */
    real min_vpara;   /**< value of lowest v_parallel bin */
    real max_vpara;   /**< value of highest v_parallel bin */
    int n_vperp;      /**< number of v_perpendicular bins */
    real min_vperp;   /**< value of lowest v_perpendicular bin */
    real max_vperp;   /**< value of highest v_perpendicular bin */
    int offload_array_length; /**< number of elements in offload_array */
} dist_rzvv_offload_data;

typedef struct {
    int n_r;          /**< number of r bins */ 
    real min_r;       /**< value of lowest r bin */
    real max_r;       /**< value of highest r bin */
    int n_z;          /**< number of z bins */
    real min_z;       /**< value of lowest z bin */
    real max_z;       /**< value of highest z bin */
    int n_vpara;      /**< number of v_parallel bins */
    real min_vpara;   /**< value of lowest v_parallel bin */
    real max_vpara;   /**< value of highest v_parallel bin */
    int n_vperp;      /**< number of v_perpendicular bins */
    real min_vperp;   /**< value of lowest v_perpendicular bin */
    real max_vperp;   /**< value of highest v_perpendicular bin */
    real* histogram;  /**< pointer to start of histogram array */
} dist_rzvv_data;

void dist_rzvv_init_offload(dist_rzvv_offload_data* offload_data,
                            real** offload_array);
void dist_rzvv_free_offload(dist_rzvv_offload_data* offload_data,
                            real** offload_array);

void dist_rzvv_print_rz(dist_rzvv_offload_data* dist, real* histogram);
void dist_rzvv_print_vv(dist_rzvv_offload_data* dist, real* histogram);

void dist_rzvv_sum(dist_rzvv_offload_data* dist1, real* array1, real* array2);

#pragma omp declare target
void dist_rzvv_init(dist_rzvv_data* dist_data,
                    dist_rzvv_offload_data* offload_data,
                    real* offload_array);
void dist_rzvv_update_fo(dist_rzvv_data* dist, particle_simd_fo* p_f, particle_simd_fo* p_i);
void dist_rzvv_update_gc(dist_rzvv_data* dist, particle_simd_gc* p_f, particle_simd_gc* p_i);
#pragma omp end declare target

#endif
