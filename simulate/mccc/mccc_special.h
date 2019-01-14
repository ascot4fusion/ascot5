/**
 * @file mccc_special.h
 * @brief Header file for mccc_special.c
*/
#ifndef MCCC_SPECIAL_H
#define MCCC_SPECIAL_H
#include "../../ascot5.h"

#pragma omp declare target
#pragma omp declare simd
void mccc_special_G(real x, real* G, int exact);

#pragma omp declare simd
void mccc_special_GdG(real x, real* GdG, int exact);

#pragma omp declare simd
void mccc_special_fo(real x, real* fdf, int exact, real* coldata);

#pragma omp declare simd
void mccc_special_mu(real u, real th, real* mu, int exact);

#pragma omp declare simd
void mccc_special_mudmu(real u, real th, real* mudmu, int exact);

#pragma omp end declare target

#endif
