/**
 * @file E_3DST.h
 * @brief Header file for E_3DST.c
 */
#ifndef E_3DST_H
#define E_3DST_H
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"

/**
 * @brief 3D electric field parameters on the host
 */
typedef struct {
    int n_r;       /**< number of r grid points in E field data               */
    int n_z;       /**< number of z grid points in E field data               */
    int n_phi;     /**< number of phi grid points in E field data             */
    int n_t;       /**< number of t grid points in E field data               */
    real r_min;    /**< minimum r coordinate in the grid in E field data      */
    real r_max;    /**< maximum r coordinate in the grid in E field data      */
    real r_grid;   /**< r grid interval (r_max-r_min)/(n_r-1) in E field data */
    real z_min;    /**< minimum z coordinate in the grid in E field data      */
    real z_max;    /**< maximum z coordinate in the grid in E field data      */
    real z_grid;   /**< z grid interval (z_max-z_min)/(n_z-1) in E field data */
    real phi_min;  /**< minimum phi coordinate in the grid in E field data    */
    real phi_max;  /**< maximum phi coordinate in the grid in E field data    */
    real phi_grid; /**< phi grid interval 2pi/(n_phi-1) in E field data       */
    real t_min;    /**< minimum t coordinate in the grid in E field data      */
    real t_max;    /**< maximum t coordinate in the grid in E field data      */

    int offload_array_length;   /**< number of elements in offload_array */
} E_3DST_offload_data;

/**
 * @brief 3D electric field parameters on the target
 */
typedef struct {
    interp4D_data E_r;   /**< pointer to E_r interpolation data struct */
    interp4D_data E_phi; /**< pointer to E_phi interpolation data struct */
    interp4D_data E_z;   /**< pointer to E_z interpolation data struct */

} E_3DST_data;

int E_3DST_init_offload(
    E_3DST_offload_data* offload_data, real** offload_array);
void E_3DST_free_offload(
    E_3DST_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void E_3DST_init(E_3DST_data* Edata, E_3DST_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Edata)
a5err E_3DST_eval_E(
    real E[3], real r, real phi, real z, real t, E_3DST_data* Edata);
#pragma omp end declare target
#endif
