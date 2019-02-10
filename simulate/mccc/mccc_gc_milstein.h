/**
 * @file mccc_gc_milstein.h
 * @brief Header file for mccc_gc_milstein.c
 */
#ifndef MCCC_GC_MILSTEIN_H
#define MCCC_GC_MILSTEIN_H

#include "../../ascot5.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"
#include "mccc_wiener.h"

#pragma omp declare target
void mccc_gc_milstein(particle_simd_gc* p, real* hin, real* hout, real tol,
                      mccc_wienarr** wienarr, B_field_data* Bdata,
                      plasma_data* pdata, random_data* rdata, real* coldata);

#pragma omp end declare target

#endif
