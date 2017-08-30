/**
 * @file interp2D.h
 * @brief Header file for interp2D.c
 */
#ifndef INTERP2D_H
#define INTERP2D_H
#include "../ascot5.h"

/**
 * @brief Interpolation struct bla bla
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
    real* c;                  /**< pointer to array with spline coefficients */
} interp2D_data;

#pragma omp declare target
void interp2D_init(interp2D_data* str, real* f, int n_r, int n_z,
		   real r_min, real r_max, real r_grid,
		   real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
void interp2D_eval_B(real* B, interp2D_data* str, real r, real z);
#pragma omp declare simd uniform(str)
void interp2D_eval_dB(real* B_dB, interp2D_data* str, real r, real z);
void interp2D_free(interp2D_data* str);
#pragma omp end declare target
#endif
