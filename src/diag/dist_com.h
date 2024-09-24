/**
 * @file dist_com.h
 * @brief Header file for dist_com.c
 */
#ifndef DIST_COM_H
#define DIST_COM_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"
#include "../B_field.h"

/**
 * @brief Histogram parameters that will be offloaded to target
 */
typedef struct {
    int n_mu;            /**< number of mu bins                     */
    real min_mu;         /**< value of lowest mu bin                */
    real max_mu;         /**< value of highest ,u bin               */

    int n_Ekin;          /**< number of Ekin bins                   */
    real min_Ekin;       /**< value of lowest Ekin bin              */
    real max_Ekin;       /**< value of highest Ekin bin             */

    int n_Ptor;          /**< number of Ptor bins                   */
    real min_Ptor;       /**< value of lowest Ptor bin              */
    real max_Ptor;       /**< value of highest Ptor bin             */
} dist_COM_offload_data;

/**
 * @brief Histogram parameters on target
 */
typedef struct {
    int n_mu;            /**< number of mu bins                     */
    real min_mu;         /**< value of lowest mu bin                */
    real max_mu;         /**< value of highest ,u bin               */

    int n_Ekin;          /**< number of Ekin bins                   */
    real min_Ekin;       /**< value of lowest Ekin bin              */
    real max_Ekin;       /**< value of highest Ekin bin             */

    int n_Ptor;          /**< number of Ptor bins                   */
    real min_Ptor;       /**< value of lowest Ptor bin              */
    real max_Ptor;       /**< value of highest Ptor bin             */

    size_t step_1;       /**< step for 2nd fastest running index    */
    size_t step_2;       /**< step for 3rd fastest running index    */

    real* histogram;  /**< pointer to start of histogram array */
} dist_COM_data;

#pragma omp declare target
void dist_COM_init(dist_COM_data* dist_data,
                   dist_COM_offload_data* offload_data,
                   real* offload_array);
void dist_COM_update_fo(dist_COM_data* dist, B_field_data*Bdata,
                        particle_simd_fo* p_f, particle_simd_fo* p_i, particle_loc* p_loc);
void dist_COM_update_gc(dist_COM_data* dist, B_field_data* Bdata,
                        particle_simd_gc* p_f, particle_simd_gc* p_i);
#pragma omp end declare target

#endif
