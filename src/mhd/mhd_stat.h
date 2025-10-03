/**
 * @file mhd_stat.h
 * @brief Header file for mhd_stat.c
 */
#ifndef MHD_STAT_H
#define MHD_STAT_H

#include "B_field.h"
#include "ascot5.h"
#include "boozer.h"
#include "error.h"
#include "mhd.h"
#include "offload.h"

int mhd_stat_init(
    mhd_stat_data *data, int nmode, int nrho, real rhomin, real rhomax,
    int *moden, int *modem, real *amplitude_nm, real *omega_nm, real *phase_nm,
    real *alpha, real *phi);

void mhd_stat_free(mhd_stat_data *data);

void mhd_stat_offload(mhd_stat_data *data);

DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, includemode)
a5err mhd_stat_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, int includemode,
    boozer_data *boozerdata, mhd_stat_data *mhddata, B_field_data *Bdata);

DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, pertonly, includemode)
a5err mhd_stat_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    int includemode, boozer_data *boozerdata, mhd_stat_data *mhddata,
    B_field_data *Bdata);

#endif
