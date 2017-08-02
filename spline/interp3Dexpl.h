/**
 * @file interp3Dexpl.h
 * @brief Header file for interp3D.c
 */
#include "../ascot5.h"
#include "interp3D.h"

void interp3Dexpl_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
		       real r_min, real r_max, real r_grid,
		       real phi_min, real phi_max, real phi_grid,
		       real z_min, real z_max, real z_grid);
void interp3Dexpl_eval_B(real* B, interp3D_data* str, real r, real phi, real z);
void interp3Dexpl_eval_dB(real* B_dB, interp3D_data* str, real r, real phi, real z);
void interp3Dexpl_free(interp3D_data* str);