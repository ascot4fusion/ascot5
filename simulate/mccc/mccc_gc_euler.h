/**
 * @file mccc_gc_euler.h
 * @brief Header file for mccc_gc_euler.c
 */
#ifndef MCCC_GC_EULER_H
#define MCCC_GC_EULER_H

#include "../../ascot5.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"

#pragma omp declare target
void mccc_gc_euler(particle_simd_gc* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, real* coldata);

#pragma omp end declare target

#endif
