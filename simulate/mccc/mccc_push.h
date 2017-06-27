/**
 * @file mccc_push.h
 * @brief Header file for mccc_push.c
*/
#ifndef MCCC_PUSH_H
#define MCCC_PUSH_H

#include <math.h>
#include "../../ascot5.h"
#include "../../consts.h"

#define MCCC_PUSH_ISNAN 20
#define MCCC_PUSH_NOTPHYSICAL 21

#pragma omp declare target
#pragma omp declare simd
void mccc_push_foEM(real F, real Dpara, real Dperp, real dt, real* rnd, real* vin, real* vout, int* err);

#pragma omp declare simd
void mccc_push_gcEM(real K, real mu, real Dpara, real DX, real* B, real dt, real* rnd, 
		    real vin, real* vout, real xiin, real* xiout, real* Xin, real* Xout, real cutoff, int* err);

#pragma omp declare simd
void mccc_push_gcMI(real K, real nu, real Dpara, real DX, real* B, real dt, real* dW, real dQ, real dDpara, real vin, real* vout, 
		    real xiin, real* xiout, real* Xin, real* Xout, real cutoff, real tol, real* kappa_k, real* kappa_d, int* err);

#pragma omp end declare target

#endif
