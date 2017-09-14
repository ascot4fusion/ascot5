/**
 * @brief Header file for step_gceom.c
 */
#ifndef STEP_GCEOM_H
#define STEP_GCEOM_H

#include "../../ascot5.h"

#pragma omp declare target
#pragma omp simd uniform(y, B_dB, E)
void step_gceom(real* ydot, real* y, real mass, real charge, real* B_dB, real* E);
#pragma omp end declare target

#endif
