/**
 * @file step_gc_rk4.h
 * @brief Header file for step_gc_rk4.c
 */
#ifndef STEP_GC_RK4_H
#define STEP_GC_RK4_H

#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../boozer.h"
#include "../../mhd.h"
#include "../../particle.h"

#pragma omp declare target
void step_gc_rk4(particle_simd_gc* p, real* h, B_field_data* Bdata,
                 E_field_data* Edata);
void step_gc_rk4_mhd(particle_simd_gc* p, real* h, B_field_data* Bdata,
                     E_field_data* Edata, boozer_data* boozer,
                     mhd_data* mhd);
#pragma omp end declare target

#endif
