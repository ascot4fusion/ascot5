/**
 * @file dist_rho5D.h
 * @brief Header file for dist_rho5D.c
 */
#ifndef DIST_RHO5D_H
#define DIST_RHO5D_H

#include "ascot5.h"
#include "particle.h"

/**
 * @brief Histogram parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;          /**< number of rho bins */
    real min_rho;       /**< value of lowest rho bin */
    real max_rho;       /**< value of highest rho bin */
    
	int n_pol;          /**< number of poloidal angle bins */
    real min_pol;       /**< value of lowest pol bin */
    real max_pol;       /**< value of highest pol bin */
	
	int n_phi;          /**< number of phi bins */
    real min_phi;       /**< value of lowest phi bin */
    real max_phi;       /**< value of highest phi bin */
    
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
} dist_rho5D_offload_data;

typedef struct {
    int n_rho;          /**< number of rho bins */
    real min_rho;       /**< value of lowest rho bin */
    real max_rho;       /**< value of highest rho bin */
    
	int n_pol;          /**< number of poloidal angle bins */
    real min_pol;       /**< value of lowest pol bin */
    real max_pol;       /**< value of highest pol bin */
	
	int n_phi;          /**< number of phi bins */
    real min_phi;       /**< value of lowest phi bin */
    real max_phi;       /**< value of highest phi bin */
    
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
} dist_rho5D_data;

void dist_rho5D_sum(int start, int stop, real* array1, real* array2);

#pragma omp declare target
void dist_rho5D_init(dist_rho5D_data* dist_data,
                  dist_rho5D_offload_data* offload_data,
                  real* offload_array);
void dist_rho5D_update_fo(dist_rho5D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i);
void dist_rho5D_update_gc(dist_rho5D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i);
#pragma omp end declare target

#endif
