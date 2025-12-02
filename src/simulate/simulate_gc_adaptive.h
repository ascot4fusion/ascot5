/**
 * @file simulate_gc_adaptive.h
 * @brief Header file for simulate_gc_adaptive.c
 */
#ifndef SIMULATE_GC_ADAPTIVE_H
#define SIMULATE_GC_ADAPTIVE_H

#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"

void simulate_gc_adaptive(particle_queue* pq, sim_data* sim, int mrk_array_size);

#endif
