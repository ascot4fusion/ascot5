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
    
    int n_phi;        /**< number of phi bins */
    real min_phi;     /**< value of lowest phi bin */
    real max_phi;     /**< value of highest phi bin */
    
    int n_z;          /**< number of v_parallel bins */
    real min_z;       /**< value of lowest z bin */
    real max_z;       /**< value of highest z bin */
    
    int n_vpara;      /**< number of v_parallel bins */
    real min_vpara;   /**< value of lowest v_parallel bin */
    real max_vpara;   /**< value of highest v_parallel bin */
    
    int n_vperp;      /**< number of v_perpendicular bins */
    real min_vperp;   /**< value of lowest v_perpendicular bin */
    real max_vperp;   /**< value of highest v_perpendicular bin */
    
    int n_time;       /**< number of time bins */ 
    real min_time;    /**< value of lowest time bin */
    real max_time;    /**< value of highest time bin */

    int n_q;          /**< number of charge bins */ 
    real min_q;       /**< value of lowest charge bin */
    real max_q;       /**< value of highest charge bin */
} dist_rzvv_offload_data;

typedef struct {
    int n_r;          /**< number of r bins */ 
    real min_r;       /**< value of lowest r bin */
    real max_r;       /**< value of highest r bin */

    int n_phi;        /**< number of r bins */ 
    real min_phi;     /**< value of lowest r bin */
    real max_phi;     /**< value of highest r bin */
    
    int n_z;          /**< number of z bins */
    real min_z;       /**< value of lowest z bin */
    real max_z;       /**< value of highest z bin */
    
    int n_vpara;      /**< number of v_parallel bins */
    real min_vpara;   /**< value of lowest v_parallel bin */
    real max_vpara;   /**< value of highest v_parallel bin */
    
    int n_vperp;      /**< number of v_perpendicular bins */
    real min_vperp;   /**< value of lowest v_perpendicular bin */
    real max_vperp;   /**< value of highest v_perpendicular bin */

    int n_time;       /**< number of r bins */ 
    real min_time;    /**< value of lowest r bin */
    real max_time;    /**< value of highest r bin */

    int n_q;          /**< number of r bins */ 
    real min_q;       /**< value of lowest r bin */
    real max_q;       /**< value of highest r bin */
    
    real* histogram;  /**< pointer to start of histogram array */
} dist_rzvv_data;

void dist_rzvv_print_rz(dist_rzvv_offload_data* dist, real* histogram);
void dist_rzvv_print_vv(dist_rzvv_offload_data* dist, real* histogram);

void dist_rzvv_sum(int start, int stop, real* array1, real* array2);

#pragma omp declare target
void dist_rzvv_init(dist_rzvv_data* dist_data,
                    dist_rzvv_offload_data* offload_data,
                    real* offload_array);
void dist_rzvv_update_fo(dist_rzvv_data* dist, particle_simd_fo* p_f, particle_simd_fo* p_i);
void dist_rzvv_update_gc(dist_rzvv_data* dist, particle_simd_gc* p_f, particle_simd_gc* p_i);
#pragma omp end declare target

#endif
