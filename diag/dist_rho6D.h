/**
 * @file dist_rho6D.h
 * @brief Header file for dist_rho6D.c
 */
#ifndef DIST_RHO6D_H
#define DIST_RHO6D_H

#include "../ascot5.h"
#include "../particle.h"

/**
 * @brief Histogram parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;        /**< number of rho bins */
    real min_rho;     /**< value of lowest rho bin */
    real max_rho;     /**< value of highest rho bin */

    int n_pol;        /**< number of poloidal angle bins */
    real min_pol;     /**< value of lowest pol bin */
    real max_pol;     /**< value of highest pol bin */

    int n_phi;        /**< number of phi bins */
    real min_phi;     /**< value of lowest phi bin */
    real max_phi;     /**< value of highest phi bin */

    int n_vr;         /**< number of v_r bins */
    real min_vr;      /**< value of lowest v_r bin */
    real max_vr;      /**< value of highest v_r bin */

    int n_vphi;       /**< number of v_phi bins */
    real min_vphi;    /**< value of lowest v_phi bin */
    real max_vphi;    /**< value of highest v_phi bin */

    int n_vz;         /**< number of v_z bins */
    real min_vz;      /**< value of lowest v_z bin */
    real max_vz;      /**< value of highest v_z bin */

    int n_time;       /**< number of time bins */
    real min_time;    /**< value of lowest time bin */
    real max_time;    /**< value of highest time bin */

    int n_q;          /**< number of charge bins */
    real min_q;       /**< value of lowest charge bin */
    real max_q;       /**< value of highest charge bin */
} dist_rho6D_offload_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct {
    int n_rho;        /**< number of rho bins */
    real min_rho;     /**< value of lowest rho bin */
    real max_rho;     /**< value of highest rho bin */

    int n_pol;        /**< number of poloidal angle bins */
    real min_pol;     /**< value of lowest pol bin */
    real max_pol;     /**< value of highest pol bin */

    int n_phi;        /**< number of phi bins */
    real min_phi;     /**< value of lowest phi bin */
    real max_phi;     /**< value of highest phi bin */

    int n_vr;         /**< number of v_r bins */
    real min_vr;      /**< value of lowest v_r bin */
    real max_vr;      /**< value of highest v_r bin */

    int n_vphi;       /**< number of v_phi bins */
    real min_vphi;    /**< value of lowest v_phi bin */
    real max_vphi;    /**< value of highest v_phi bin */

    int n_vz;         /**< number of v_z bins */
    real min_vz;      /**< value of lowest v_z bin */
    real max_vz;      /**< value of highest v_z bin */

    int n_time;       /**< number of time bins */
    real min_time;    /**< value of lowest time bin */
    real max_time;    /**< value of highest time bin */

    int n_q;          /**< number of charge bins */
    real min_q;       /**< value of lowest charge bin */
    real max_q;       /**< value of highest charge bin */

    real* histogram;  /**< pointer to start of histogram array */
} dist_rho6D_data;

#pragma omp declare target
void dist_rho6D_init(dist_rho6D_data* dist_data, dist_rho6D_offload_data* offload_data,
                     real* offload_array);
void dist_rho6D_update_fo(dist_rho6D_data* dist, particle_simd_fo* p_f,
                          particle_simd_fo* p_i);
void dist_rho6D_update_gc(dist_rho6D_data* dist, particle_simd_gc* p_f,
                          particle_simd_gc* p_i);
#pragma omp end declare target

#endif
