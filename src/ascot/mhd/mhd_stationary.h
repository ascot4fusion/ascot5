/**
 * @file mhd_stationary.h
 * Header file for mhd_stat.c
 */
#ifndef MHD_STAT_H
#define MHD_STAT_H

#include "bfield.h"
#include "defines.h"
#include "boozer.h"
#include "errors.h"
#include "mhd.h"
#include "parallel.h"

int mhd_stat_init(
    mhd_stat_data *data, int nmode, int nrho, real rhomin, real rhomax,
    int *moden, int *modem, real *amplitude_nm, real *omega_nm, real *phase_nm,
    real *alpha, real *phi);

void mhd_stat_free(mhd_stat_data *data);

void mhd_stat_offload(mhd_stat_data *data);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, includemode)
a5err mhd_stat_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, int includemode,
    Boozer *boozer, mhd_stat_data *mhd, Bfield *bfield);

DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, pertonly, includemode)
a5err mhd_stat_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    int includemode, Boozer *boozer, mhd_stat_data *mhd,
    Bfield *bfield);

#endif
