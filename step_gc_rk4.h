/**
 * @file step_gc_rk4.h
 * @brief Header file for step_gc_rk4.c
 */
#ifndef STEP_GC_RK4_H
#define STEP_GC_RK4_H

#include "B_field.h"
#include "particle.h"

#pragma omp declare target
#pragma omp declare simd 
void ydot_gc(real ydot[], real t, real y[], real mass[], real charge[],
             real B_dB[]);
void step_gc_rk4(particle_simd_gc* p, real t, real h, B_field_data* Bdata);
#pragma omp end declare target

#endif
