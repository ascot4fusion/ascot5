/**
 * @file mhd_dynamic.h
 * Mhd perturbation where eigenmodes evolve in time.
 */
#ifndef MHD_DYNAMIC_H
#define MHD_DYNAMIC_H

#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd.h"
#include "parallel.h"

int MhdDynamic_init(
    MhdDynamic *mhd, size_t n, size_t nrho, size_t ntime, int nmode[n],
    int mmode[n], real rholim[2], real tlim[2], real amplitude[n],
    real omega[n], real phase[n], real alpha[n * nrho * ntime],
    real phi[n * nrho * ntime]);

void MhdDynamic_free(MhdDynamic *mhd);

void MhdDynamic_offload(MhdDynamic *mhd);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, includemode)
err_t MhdDynamic_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, size_t includemode,
    Boozer *boozer, MhdDynamic *mhd, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, pertonly, includemode)
err_t MhdDynamic_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    size_t includemode, Boozer *boozer, MhdDynamic *mhd, Bfield *bfield);

#endif
