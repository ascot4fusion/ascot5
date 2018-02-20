/**
 * @file linint3D.h
 * @brief Header file for linint3D.c
 */
#ifndef LININT3D_H
#define LININT3D_H
#include "../ascot5.h"

/**
 * @brief 3D linear interpolation struct
 */
typedef struct {
    int n_r;                  /**< number of r grid points */
    int n_z;                  /**< number of z grid points */
    int n_phi;                /**< number of phi grid points */
    real r_min;               /**< minimum r coordinate in the grid */
    real r_max;               /**< r grid interval (r_max-r_min)/(n_r-1) */
    real r_grid;              /**< r grid interval (r_max-r_min)/(n_r-1) */
    real z_min;               /**< minimum z coordinate in the grid */
    real z_max;               /**< z grid interval (z_max-z_min)/(n_z-1) */
    real z_grid;              /**< z grid interval (z_max-z_min)/(n_z-1) */
    real phi_min;             /**< minimum phi coordinate in the grid */
    real phi_max;             /**< z coordinate of magnetic axis */
    real phi_grid;            /**< phi grid interval 2pi/(n_phi-1) */
    real* f;                   /**< pointer to array with function values */
} linint3D_data;

#pragma omp declare target
int linint3D_init(linint3D_data* str, real* f, int n_r, int n_phi, int n_z,
		   real r_min, real r_max, real r_grid,
		   real phi_min, real phi_max, real phi_grid,
		   real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
int linint3D_eval(real* f, linint3D_data* str, real r, real phi, real z);
#pragma omp declare simd linear(i) uniform(f, str)
int linint3D_eval_SIMD(int i, real f[NSIMD], linint3D_data* str, real r, real phi, real z);
void linint3D_free(linint3D_data* str);
#pragma omp end declare target
#endif
