/**
 * @file interp2Detoc.h
 * @brief Header file for interp2Detoc.c
 */
#ifndef INTERP2DETOC_H
#define INTERP2DETOC_H
#include "../ascot5.h"
#include "interp2D.h"

#pragma omp declare target
void interp2Detoc_init(interp2D_data* str, real* f, int n_r, int n_z,
                       real r_min, real r_max, real r_grid,
                       real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
void interp2Detoc_eval_B(real* B, interp2D_data* str, real r, real z);
#pragma omp declare simd uniform(str)
void interp2Detoc_eval_dB(real* B_dB, interp2D_data* str, real r, real z);
void interp2Detoc_free(interp2D_data* str);
#pragma omp end declare target
#endif
