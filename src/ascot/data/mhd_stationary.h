/**
 * @file mhd_stationary.h
 * Header file for mhd_stat.c
 */
#ifndef MHD_STATIONARY_H
#define MHD_STATIONARY_H

#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd.h"
#include "parallel.h"

int MhdStationary_init(
    MhdStationary *mhd, size_t n, size_t nrho, int nmode[n], int mmode[n],
    real rholim[2], real amplitude[n], real omega[n], real phase[n],
    real alpha[n * nrho], real phi[n * nrho]);

void MhdStationary_free(MhdStationary *mhd);

void MhdStationary_offload(MhdStationary *mhd);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, includemode)
err_t MhdStationary_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, size_t includemode,
    Boozer *boozer, MhdStationary *mhd, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, pertonly, includemode)
err_t MhdStationary_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    size_t includemode, Boozer *boozer, MhdStationary *mhd, Bfield *bfield);

#endif
