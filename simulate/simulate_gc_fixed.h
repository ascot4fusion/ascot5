/**
 * @file simulate_gc_fixed.h
 * @brief Header file for simulate_gc_fixed.c
 */
#ifndef SIMULATE_GC_FIXED_H
#define SIMULATE_GC_FIXED_H

#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"

#pragma omp declare target
void simulate_gc_fixed(particle_queue_gc* pq, sim_data* sim);
#pragma omp end declare target

#endif
