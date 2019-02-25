/**
 * @file gctransform.h
 * @brief Header file for gctransform.c
 */
#ifndef GCTRANSFORM_H
#define GCTRANSFORM_H

#include "ascot5.h"

#pragma omp declare target

void gctransform_setorder(int order);

#pragma omp declare simd
void gctransform_particle2guidingcenter(
    real mass, real charge, real* B_dB,
    real r, real phi, real z, real vr, real vphi, real vz,
    real* R, real* Phi, real* Z, real* vpar, real* mu, real* zeta);

#pragma omp declare simd
void gctransform_guidingcenter2particle(
    real mass, real charge, real* B_dB,
    real R, real Phi, real Z, real vpar, real mu, real zeta,
    real* r, real* phi, real* z, real* vparprt, real* muprt, real* zetaprt);

#pragma omp declare simd
void gctransform_vparmuzeta2vRvphivz(real mass, real charge, real* B_dB,
                                     real phi, real vpar, real mu, real zeta,
                                     real* vr, real* vphi, real* vz);
#pragma omp end declare target

#endif
