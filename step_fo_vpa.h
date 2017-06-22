/**
 * @file step_fo_vpa.h
 * @brief Header file for step_fo_vpa.c
 */
#ifndef STEP_FO_VPA_H
#define STEP_FO_VPA_H
#include "B_field.h"
#include "E_field.h"
#include "particle.h"

#pragma omp declare target
void step_fo_vpa(particle_simd_fo* p, real* h, B_field_data* Bdata, E_field_data* Edata);
#pragma omp end declare target

#endif
