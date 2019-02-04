/**
 * @file E_3D.h
 * @brief Header file for E_3D.c
 */
#ifndef E_3D_H
#define E_3D_H
#include "../ascot5.h"
#include "../error.h"
#include "../linint/linint3D.h" /* for 3D interpolation routines */

/**
 * @brief 3D electric field parameters on the host
 */
typedef struct {
    int n_r;              /**< number of r grid points in electric field data */
    int n_z;              /**< number of z grid points in electric field data */
    int n_phi;            /**< number of phi grid points in electric field data */
    real r_min;           /**< minimum r coordinate in the grid in electric field data */
    real r_max;           /**< maximum r coordinate in the grid in electric field data */
    real r_grid;          /**< r grid interval (r_max-r_min)/(n_r-1) in electric field data */
    real z_min;           /**< minimum z coordinate in the grid in electric field data */
    real z_max;           /**< maximum z coordinate in the grid in electric field data */
    real z_grid;          /**< z grid interval (z_max-z_min)/(n_z-1) in electric field data */
    real phi_min;               /**< minimum phi coordinate in the grid in electric field data */
    real phi_max;               /**< maximum phi coordinate in the grid in electric field data */
    real phi_grid;              /**< phi grid interval 2pi/(n_phi-1) in electric field data */
    int offload_array_length;   /**< number of elements in offload_array */
} E_3D_offload_data;

/**
 * @brief 3D electric field parameters on the target
 */
typedef struct {
    linint3D_data E_r;     /**< pointer to start of E_r interpolation data struct */
    linint3D_data E_phi;     /**< pointer to start of E_phi interpolation data struct */
    linint3D_data E_z;     /**< pointer to start of E_z interpolation data struct */

} E_3D_data;

int E_3D_init_offload(E_3D_offload_data* offload_data, real** offload_array);
void E_3D_free_offload(E_3D_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void E_3D_init(E_3D_data* Edata, E_3D_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Edata)
a5err E_3D_eval_E(real E[3], real r, real phi, real z, E_3D_data* Edata);
#pragma omp end declare target   
#endif