/**
 * @file interp1Dcomp.h
 * @brief Header file for interp1Dcomp.c
 */
#ifndef INTERP1DCOMP_H
#define INTERP1DCOMP_H
#include "../ascot5.h"

/**
 * @brief Interpolation struct
 */
typedef struct {
    int n_r;                  /**< number of r grid points */
    real r_min;               /**< minimum r coordinate in the grid */
    real r_max;               /**< r grid interval (r_max-r_min)/(n_r-1) */
    real r_grid;              /**< r grid interval (r_max-r_min)/(n_r-1) */
    real* c;                  /**< pointer to array with spline coefficients */
} interp1D_data;

#pragma omp declare target
int interp1Dcomp_init(interp1D_data* str, real* f, int n_r,
		      real r_min, real r_max, real r_grid);
#pragma omp declare simd uniform(str) simdlen(8)
int interp1Dcomp_eval_B(real* B, interp1D_data* str, real r);
#pragma omp declare simd linear(i) uniform(B, str) simdlen(8)
int interp1Dcomp_eval_B_SIMD(int i, real B[NSIMD], interp1D_data* str, real r);
#pragma omp declare simd uniform(str) simdlen(8)
int interp1Dcomp_eval_dB(real* B_dB, interp1D_data* str, real r);
#pragma omp declare simd linear(i) uniform(B_dB, str) simdlen(8)
int interp1Dcomp_eval_dB_SIMD(int i, real B_dB[6][NSIMD], interp1D_data* str, real r);
void interp1Dcomp_free(interp1D_data* str);
#pragma omp end declare target
#endif
