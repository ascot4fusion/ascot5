/**
 * @file step_fo_lf.h
 * @brief Header file for step_fo_lf.c
 */
#ifndef STEP_FO_LF_H
#define STEP_FO_LF_H
#include "B_field.h"
#include "particle.h"

#pragma omp declare target
void step_fo_lf(particle_simd_fo* p, real t, real h, B_field_data* Bdata);
#pragma omp end declare target

#endif
