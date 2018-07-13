/**
 * @file linint1D.h
 * @brief Header file for linint1D.c
 */
#ifndef LININT1D_H
#define LININT1D_H
#include "../ascot5.h"

/**
 * @brief 1D linear interpolation struct
 */
typedef struct {
    int n_r;                  /**< number of r grid points */
    real r_min;               /**< minimum r coordinate in the grid */
    real r_max;               /**< r grid interval (r_max-r_min)/(n_r-1) */
    real r_grid;              /**< r grid interval (r_max-r_min)/(n_r-1) */
    real* f;                   /**< pointer to array with function values */
} linint1D_data;

#pragma omp declare target
int linint1D_init(linint1D_data* str, real* f, int n_r,
                  real r_min, real r_max, real r_grid);
#pragma omp declare simd uniform(str)
integer linint1D_eval(real* f, linint1D_data* str, real r);
#pragma omp declare simd linear(i) uniform(f, str)
integer linint1D_eval_SIMD(int i, real f[NSIMD], linint1D_data* str, real r);
void linint1D_free(linint1D_data* str);
#pragma omp end declare target
#endif
