/**
 * @file flr_losses.h
 * @brief Header file for flr_losses.c. Contains the routines to 
 * determine the FLR losses of GC markers.
 */
#ifndef FLR_LOSSES_H
#define FLR_LOSSES_H

#include "../ascot5.h"
#include "../offload.h"
#include "../B_field.h"
#include "../wall.h"
#include "../random.h"
#include "../error.h"


GPU_DECLARE_TARGET_SIMD_UNIFORM(B, wall, rnd)
int flr_losses_eval(real r, real phi, real z, real ppar, real mu,
					real mass, real charge, real time,
					B_field_data* B, wall_data* wall, random_data* rnd,
					int* walltile_out, int* err_out);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

#endif 