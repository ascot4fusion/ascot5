/**
 * @file interp2Dcomp.h
 * @brief Header file for interp2Dcomp.c
 */
#ifndef INTERP2DCOMP_H
#define INTERP2DCOMP_H
#include "../ascot5.h"
#include "interp2D.h"

#pragma omp declare target
int interp2Dcomp_init(interp2D_data* str, real* f, int n_r, int n_z,
		      real r_min, real r_max, real r_grid,
		      real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
integer interp2Dcomp_eval_B(real* B, interp2D_data* str, real r, real z);
#pragma omp declare simd linear(i) uniform(B, str)
integer interp2Dcomp_eval_B_SIMD(int i, real B[NSIMD], interp2D_data* str, real r, real z);
#pragma omp declare simd uniform(str)
integer interp2Dcomp_eval_dB(real* B_dB, interp2D_data* str, real r, real z);
#pragma omp declare simd linear(i) uniform(B_dB, str)
integer interp2Dcomp_eval_dB_SIMD(int i, real B_dB[6][NSIMD], interp2D_data* str, real r, real z);
void interp2Dcomp_free(interp2D_data* str);
#pragma omp end declare target
#endif
