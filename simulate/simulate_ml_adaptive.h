/**
 * @file simulate_ml_adaptive.h
 * @brief Header file for simulate_ml_adaptive.c
 */
#ifndef SIMULATE_ML_ADAPTIVE_H
#define SIMULATE_ML_ADAPTIVE_H

#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"

#pragma omp declare target
void simulate_ml_adaptive(int id, int n_particles, particle* particles,
			  sim_offload_data sim,
			  real* B_offload_array,
			  real* E_offload_array,
			  real* plasma_offload_array,
			  real* wall_offload_array,
			  real* dist_offload_array);
#pragma omp end declare target

#endif
