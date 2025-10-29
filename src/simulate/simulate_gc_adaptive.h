/**
 * @file simulate_gc_adaptive.h
 * @brief Header file for simulate_gc_adaptive.c
 */
#ifndef SIMULATE_GC_ADAPTIVE_H
#define SIMULATE_GC_ADAPTIVE_H

#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"

typedef struct {
    unsigned int crossed_once : 1;
    unsigned int crossed_twice : 1;
    unsigned int first_ppar : 1;
} Crossing;

typedef struct {
    /** Acceleration factor. */
    real acc[NSIMD];

    /** Orbit time [s]. */
    real orbittime[NSIMD];

    /** Collision frequency [1/s]. */
    real collfreq[NSIMD];

    /** Storage for OMP crossing data. */
    Crossing cross[NSIMD];
} Acceleration;

void recalculate_acceleration(
    Acceleration* acc, sim_data* sim, particle_simd_gc* p, particle_simd_gc* p0);

void simulate_gc_adaptive(particle_queue* pq, sim_data* sim);

#endif
