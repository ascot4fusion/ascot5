/**
 * @file interp1Dcomp.h
 * @brief Header file for interp1Dcomp.c
 */
#ifndef INTERP1DCOMP_H
#define INTERP1DCOMP_H
#include "../ascot5.h"
#include "interp1D.h"

#pragma omp declare target
int interp1Dcomp_init(interp1D_data* str, real* f, int n_r,
		      real r_min, real r_max, real r_grid);
#pragma omp declare simd uniform(str)
integer interp1Dcomp_eval_B(real* B, interp1D_data* str, real r);
#pragma omp declare simd linear(i) uniform(B, str)
integer interp1Dcomp_eval_B_SIMD(int i, real B[NSIMD], interp1D_data* str, real r);
#pragma omp declare simd uniform(str)
integer interp1Dcomp_eval_dB(real* B_dB, interp1D_data* str, real r);
#pragma omp declare simd linear(i) uniform(B_dB, str)
integer interp1Dcomp_eval_dB_SIMD(int i, real B_dB[6][NSIMD], interp1D_data* str, real r);
void interp1Dcomp_free(interp1D_data* str);
#pragma omp end declare target

#endif
