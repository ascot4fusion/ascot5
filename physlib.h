/**
 * @brief Header file for physlib.c
 */
#ifndef PHYSLIB_H
#define PHYSLIB_H

#include "ascot5.h"

#pragma omp declare target

#pragma omp declare simd
real physlib_relfactorv_fo(real vnorm);

#pragma omp declare simd
real physlib_relfactorp_fo(real mass, real pnorm);

#pragma omp declare simd
real physlib_relfactorv_gc(real mass, real mu, real vpar, real Bnorm);

#pragma omp declare simd
real physlib_relfactorp_gc(real mass, real mu, real ppar, real Bnorm);

#pragma omp declare simd
real physlib_Ekin_fo(real mass, real vtot);

#pragma omp declare simd
real physlib_Ekin_gc(real mass, real mu, real vpar, real Bnorm);

#pragma omp declare simd
void physlib_gc_vxi2muvpar(real mass, real Bnorm, real v, real xi, real* mu, real* vpar);

#pragma omp declare simd
void physlib_gc_muvpar2vxi(real mass, real Bnorm, real mu, real vpar, real* v, real* xi);

#pragma omp declare simd
void physlib_fo2gc(real mass, real charge, real* B_dB,
		   real Rprt, real phiprt, real zprt, real pR, real pphi, real pz,
		   real* R, real* phi, real* z, real* mu, real* ppar, real* theta);

#pragma omp declare simd
void physlib_gc2fo(real mass, real charge, real* B_dB,
		   real R, real phi, real z, real mu, real ppar, real theta,
		   real* Rprt, real* phiprt, real* zprt, real* pR, real* pphi, real* pz);
#pragma omp end declare target

#endif
