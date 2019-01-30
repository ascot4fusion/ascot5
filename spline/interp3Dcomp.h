/**
 * @file interp3Dcomp.h
 * @brief Header file for interp3Dcomp.c
 */
#ifndef INTERP3DCOMP_H
#define INTERP3DCOMP_H
#include "../ascot5.h"
#include "interp3D.h"

#pragma omp declare target
int interp3Dcomp_init(interp3D_data* str, real* f, int n_x, int n_y, int n_z,
                      real x_min, real x_max, real x_grid,
                      real y_min, real y_max, real y_grid,
                      real z_min, real z_max, real z_grid);
#pragma omp declare simd uniform(str)
integer interp3Dcomp_eval_f(real* f, interp3D_data* str, real x, real y, real z);
#pragma omp declare simd uniform(str)
integer interp3Dcomp_eval_df(real* f_df, interp3D_data* str, real x, real y, real z);
void interp3Dcomp_free(interp3D_data* str);
#pragma omp end declare target
#endif
