/**
 * @file dist_rho6D.h
 * @brief Header file for dist_rho6D.c
 */
#ifndef DIST_RHO6D_H
#define DIST_RHO6D_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"

/**
 * @brief Histogram parameters that will be offloaded to target
 */
typedef struct {
    int n_rho;        /**< number of rho bins            */
    real min_rho;     /**< value of lowest rho bin       */
    real max_rho;     /**< value of highest rho bin      */

    int n_theta;      /**< number of poloidal angle bins */
    real min_theta;   /**< value of lowest theta bin     */
    real max_theta;   /**< value of highest theta bin    */

    int n_phi;        /**< number of phi bins            */
    real min_phi;     /**< value of lowest phi bin       */
    real max_phi;     /**< value of highest phi bin      */

    int n_pr;         /**< number of p_r bins            */
    real min_pr;      /**< value of lowest p_r bin       */
    real max_pr;      /**< value of highest p_r bin      */

    int n_pphi;       /**< number of p_phi bins          */
    real min_pphi;    /**< value of lowest p_phi bin     */
    real max_pphi;    /**< value of highest p_phi bin    */

    int n_pz;         /**< number of p_z bins            */
    real min_pz;      /**< value of lowest p_z bin       */
    real max_pz;      /**< value of highest p_z bin      */

    int n_time;       /**< number of time bins           */
    real min_time;    /**< value of lowest time bin      */
    real max_time;    /**< value of highest time bin     */

    int n_q;          /**< number of charge bins         */
    real min_q;       /**< value of lowest charge bin    */
    real max_q;       /**< value of highest charge bin   */
} dist_rho6D_offload_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct {
    int n_rho;        /**< number of rho bins            */
    real min_rho;     /**< value of lowest rho bin       */
    real max_rho;     /**< value of highest rho bin      */

    int n_theta;      /**< number of poloidal angle bins */
    real min_theta;   /**< value of lowest theta bin     */
    real max_theta;   /**< value of highest theta bin    */

    int n_phi;        /**< number of phi bins            */
    real min_phi;     /**< value of lowest phi bin       */
    real max_phi;     /**< value of highest phi bin      */

    int n_pr;         /**< number of p_r bins            */
    real min_pr;      /**< value of lowest p_r bin       */
    real max_pr;      /**< value of highest p_r bin      */

    int n_pphi;       /**< number of p_phi bins          */
    real min_pphi;    /**< value of lowest p_phi bin     */
    real max_pphi;    /**< value of highest p_phi bin    */

    int n_pz;         /**< number of p_z bins            */
    real min_pz;      /**< value of lowest p_z bin       */
    real max_pz;      /**< value of highest p_z bin      */

    int n_time;       /**< number of time bins           */
    real min_time;    /**< value of lowest time bin      */
    real max_time;    /**< value of highest time bin     */

    int n_q;          /**< number of charge bins         */
    real min_q;       /**< value of lowest charge bin    */
    real max_q;       /**< value of highest charge bin   */

    size_t step_1;    /**< step for 2nd fastest running index   */
    size_t step_2;    /**< step for 3rd fastest running index   */
    size_t step_3;    /**< step for 4th fastest running index   */
    size_t step_4;    /**< step for 5th fastest running index   */
    size_t step_5;    /**< step for 6th fastest running index   */
    size_t step_6;    /**< step for 7th fastest running index   */
    size_t step_7;    /**< step for 8th fastest running index   */

    real* histogram;  /**< pointer to start of histogram array */
} dist_rho6D_data;

#pragma omp declare target
void dist_rho6D_init(dist_rho6D_data* dist_data, dist_rho6D_offload_data* offload_data,
                     real* offload_array);
void dist_rho6D_update_fo(dist_rho6D_data* dist, particle_simd_fo* p_f,
                          particle_simd_fo* p_i, particle_loc* p_loc);
void dist_rho6D_update_gc(dist_rho6D_data* dist, particle_simd_gc* p_f,
                          particle_simd_gc* p_i);
#pragma omp end declare target

#endif
