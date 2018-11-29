/**
 * @file gctransform.h
 * @brief Header file for gctransform.c
 */
#ifndef GCTRANSFORM_H
#define GCTRANSFORM_H

#include "ascot5.h"

#pragma omp declare target

#pragma omp declare simd
void gctransform_particle2guidingcenter(
    real mass, real charge, real* B_dB,
    real r, real phi, real z, real vr, real vphi, real vz,
    real* R, real* Phi, real* Z, real* vpar, real* mu, real* theta);

#pragma omp declare simd
void gctransform_guidingcenter2particle(
    real mass, real charge, real* B_dB,
    real R, real Phi, real Z, real vpar, real mu, real theta,
    real* r, real* phi, real* z, real* vparprt, real* muprt, real* thetaprt);

#pragma omp declare simd
void gctransform_vparmutheta2vRvphivz(real mass, real charge, real* B_dB,
                                      real phi, real vpar, real mu, real theta,
                                      real* vr, real* vphi, real* vz);
#pragma omp end declare target

#endif
