/**
 * @file mccc_fo_euler.h
 * @brief Header file for mccc_fo_euler.c
 */
#ifndef MCCC_FO_EULER_H
#define MCCC_FO_EULER_H

#include "../../ascot5.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"

#pragma omp declare target
void mccc_fo_euler(particle_simd_fo* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, real* coldata);

#pragma omp end declare target

#endif
