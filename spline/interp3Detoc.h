/**
 * @file interp3Detoc.h
 * @brief Header file for interp3Detoc.c
 */
#ifndef INTERP3DETOC_H
#define INTERP3DETOC_H
#include "../ascot5.h"
#include "interp3D.h"

#pragma omp declare target
void interp3Detoc_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
		   real r_min, real r_max, real r_grid,
		   real phi_min, real phi_max, real phi_grid,
		   real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
void interp3Detoc_eval_B(real* B, interp3D_data* str, real r, real phi, real z);
#pragma omp declare simd uniform(str)
void interp3Detoc_eval_dB(real* B_dB, interp3D_data* str, real r, real phi, real z);
void interp3Detoc_free(interp3D_data* str);
#pragma omp end declare target
#endif
