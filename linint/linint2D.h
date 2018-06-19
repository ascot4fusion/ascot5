/**
 * @file linint2D.h
 * @brief Header file for linint2D.c
 */
#ifndef LININT2D_H
#define LININT2D_H
#include "../ascot5.h"

/**
 * @brief 2D linear interpolation struct
 */
typedef struct {
    int n_r;                  /**< number of r grid points */
    int n_z;                  /**< number of z grid points */
    real r_min;               /**< minimum r coordinate in the grid */
    real r_max;               /**< r grid interval (r_max-r_min)/(n_r-1) */
    real r_grid;              /**< r grid interval (r_max-r_min)/(n_r-1) */
    real z_min;               /**< minimum z coordinate in the grid */
    real z_max;               /**< z grid interval (z_max-z_min)/(n_z-1) */
    real z_grid;              /**< z grid interval (z_max-z_min)/(n_z-1) */
    real* f;                   /**< pointer to array with function values */
} linint2D_data;

#pragma omp declare target
int linint2D_init(linint2D_data* str, real* f, int n_r, int n_z,
                  real r_min, real r_max, real r_grid,
                  real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
int linint2D_eval(real* f, linint2D_data* str, real r, real z);
#pragma omp declare simd linear(i) uniform(f, str)
int linint2D_eval_SIMD(int i, real f[NSIMD], linint2D_data* str, real r, real z);
void linint2D_free(linint2D_data* str);
#pragma omp end declare target
#endif
