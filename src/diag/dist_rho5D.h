/**
 * @file dist_rho5D.h
 * @brief Header file for dist_rho5D.c
 */
#ifndef DIST_RHO5D_H
#define DIST_RHO5D_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"

/**
 * @brief Histogram parameters
 */
typedef struct {
    int n_rho;        /**< number of rho bins                   */
    real min_rho;     /**< value of lowest rho bin              */
    real max_rho;     /**< value of highest rho bin             */

    int n_theta;      /**< number of poloidal angle bins        */
    real min_theta;   /**< value of lowest pol bin              */
    real max_theta;   /**< value of highest pol bin             */

    int n_phi;        /**< number of phi bins                   */
    real min_phi;     /**< value of lowest phi bin              */
    real max_phi;     /**< value of highest phi bin             */

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

    size_t step_1;    /**< step for 2nd fastest running index   */
    size_t step_2;    /**< step for 3rd fastest running index   */
    size_t step_3;    /**< step for 4th fastest running index   */
    size_t step_4;    /**< step for 5th fastest running index   */
    size_t step_5;    /**< step for 6th fastest running index   */
    size_t step_6;    /**< step for 7th fastest running index   */

    real* histogram;  /**< pointer to start of histogram array */
} dist_rho5D_data;

int dist_rho5D_init(dist_rho5D_data* data);
void dist_rho5D_free(dist_rho5D_data* data);
void dist_rho5D_update_fo(dist_rho5D_data* dist, particle_simd_fo* p_f,
                          particle_simd_fo* p_i);
void dist_rho5D_update_gc(dist_rho5D_data* dist, particle_simd_gc* p_f,
                          particle_simd_gc* p_i);

#endif
