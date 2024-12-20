/**
 * @file dist_6D.h
 * @brief Header file for dist_6D.c
 */
#ifndef DIST_6D_H
#define DIST_6D_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"

/**
 * @brief Histogram parameters on target
 */
typedef struct {
    int n_r;          /**< number of r bins           */
    real min_r;       /**< value of lowest r bin      */
    real max_r;       /**< value of highest r bin     */

    int n_phi;        /**< number of r bins           */
    real min_phi;     /**< value of lowest r bin      */
    real max_phi;     /**< value of highest r bin     */

    int n_z;          /**< number of z bins           */
    real min_z;       /**< value of lowest z bin      */
    real max_z;       /**< value of highest z bin     */

    int n_pr;         /**< number of p_r bins         */
    real min_pr;      /**< value of lowest p_r bin    */
    real max_pr;      /**< value of highest p_r bin   */

    int n_pphi;       /**< number of p_phi bins       */
    real min_pphi;    /**< value of lowest p_phi bin  */
    real max_pphi;    /**< value of highest p_phi bin */

    int n_pz;         /**< number of p_z bins         */
    real min_pz;      /**< value of lowest p_z bin    */
    real max_pz;      /**< value of highest p_z bin   */

    int n_time;       /**< number of r bins           */
    real min_time;    /**< value of lowest r bin      */
    real max_time;    /**< value of highest r bin     */

    int n_q;          /**< number of r bins           */
    real min_q;       /**< value of lowest r bin      */
    real max_q;       /**< value of highest r bin     */

    size_t step_1;    /**< step for 2nd fastest running index   */
    size_t step_2;    /**< step for 3rd fastest running index   */
    size_t step_3;    /**< step for 4th fastest running index   */
    size_t step_4;    /**< step for 5th fastest running index   */
    size_t step_5;    /**< step for 6th fastest running index   */
    size_t step_6;    /**< step for 7th fastest running index   */
    size_t step_7;    /**< step for 8th fastest running index   */

    real* histogram;  /**< pointer to start of histogram array */
} dist_6D_data;

int dist_6D_init(dist_6D_data* data);
void dist_6D_free(dist_6D_data* data);
void dist_6D_offload(dist_6D_data* data);
void dist_6D_update_fo(dist_6D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i);
void dist_6D_update_gc(dist_6D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i);

#endif
