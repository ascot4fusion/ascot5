/**
 * @file interp3Dcomp.h
 * @brief Header file for interp3Dcomp.c
 */
#ifndef INTERP3DCOMP_H
#define INTERP3DCOMP_H
#include "../ascot5.h"
#include "interp3D.h"

#pragma omp declare target
int interp3Dcomp_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
		   real r_min, real r_max, real r_grid,
		   real phi_min, real phi_max, real phi_grid,
		   real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
integer interp3Dcomp_eval_B(real* B, interp3D_data* str, real r, real phi, real z);
#pragma omp declare simd linear(i) uniform(B, str)
integer interp3Dcomp_eval_B_SIMD(int i, real B[NSIMD], interp3D_data* str, real r, real phi, real z);
#pragma omp declare simd uniform(str)
integer interp3Dcomp_eval_dB(real* B_dB, interp3D_data* str, real r, real phi, real z);
#pragma omp declare simd linear(i) uniform(B_dB, str)
integer interp3Dcomp_eval_dB_SIMD(int i, real B_dB[10][NSIMD], interp3D_data* str, real r, real phi, real z);
void interp3Dcomp_free(interp3D_data* str);
#pragma omp end declare target
#endif
