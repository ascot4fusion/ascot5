/**
 * @file gctransform.h
 * @brief Header file for gctransform.c
 */
#ifndef GCTRANSFORM_H
#define GCTRANSFORM_H

#include "ascot5.h"


void gctransform_setorder(int order);

#pragma omp declare simd
void gctransform_particle2guidingcenter(
    real mass, real charge, real* B_dB,
    real r, real phi, real z, real pr, real pphi, real pz,
    real* R, real* Phi, real* Z, real* ppar, real* mu, real* zeta);

#pragma omp declare simd
void gctransform_guidingcenter2particle(
    real mass, real charge, real* B_dB,
    real R, real Phi, real Z, real ppar, real mu, real zeta,
    real* r, real* phi, real* z, real* pparprt, real* muprt, real* zetaprt);

#pragma omp declare simd
void gctransform_pparmuzeta2prpphipz(real mass, real charge, real* B_dB,
                                     real phi, real ppar, real mu, real zeta,
                                     real* pr, real* pphi, real* pz);

#endif
