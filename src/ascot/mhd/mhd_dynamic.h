/**
 * @file mhd_dynamic.h
 * Header file for mhd_nonstat.c
 */
#ifndef MHD_DYNAMIC_H
#define MHD_DYNAMIC_H

#include "bfield.h"
#include "defines.h"
#include "boozer.h"
#include "errors.h"
#include "mhd.h"
#include "parallel.h"

int mhd_nonstat_init(
    mhd_nonstat_data *data, int nmode, int nrho, int ntime, real rhomin,
    real rhomax, real tmin, real tmax, int *moden, int *modem,
    real *amplitude_nm, real *omega_nm, real *phase_nm, real *alpha, real *phi);

void mhd_nonstat_free(mhd_nonstat_data *data);

void mhd_nonstat_offload(mhd_nonstat_data *data);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, includemode)
a5err mhd_nonstat_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, int includemode,
    Boozer *boozer, mhd_nonstat_data *mhd, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, pertonly, includemode)
a5err mhd_nonstat_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    int includemode, Boozer *boozer, mhd_nonstat_data *mhd,
    Bfield *bfield);

#endif
