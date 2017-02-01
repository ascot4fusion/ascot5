/**
 * @file simulate_gc_rk4.h
 * @brief Header file for simulate_gc_rk4.c
 */
#ifndef SIMULATE_GC_RK4_H
#define SIMULATE_GC_RK4_H

#include "ascot5.h"
#include "simulate.h"
#include "particle.h"

#pragma omp declare target
void simulate_gc_rk4(int id, int n_particles, particle* particles,
                    sim_offload_data sim,
                    real* B_offload_array,
                    real* plasma_offload_array,
                    real* wall_offload_array,
                    real* dist_offload_array);
#pragma omp end declare target

#endif
