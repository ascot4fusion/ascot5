/**
 * @file mhd_nonstat.h
 * @brief Header file for mhd_nonstat.c
 */
#ifndef MHD_NONSTAT_H
#define MHD_NONSTAT_H

#include "B_field.h"
#include "ascot5.h"
#include "boozer.h"
#include "error.h"
#include "mhd.h"
#include "offload.h"

int mhd_nonstat_init(
    mhd_nonstat_data *data, int nmode, int nrho, int ntime, real rhomin,
    real rhomax, real tmin, real tmax, int *moden, int *modem,
    real *amplitude_nm, real *omega_nm, real *phase_nm, real *alpha, real *phi);

void mhd_nonstat_free(mhd_nonstat_data *data);

void mhd_nonstat_offload(mhd_nonstat_data *data);

DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, includemode)
a5err mhd_nonstat_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, int includemode,
    boozer_data *boozerdata, mhd_nonstat_data *mhddata, B_field_data *Bdata);

DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, pertonly, includemode)
a5err mhd_nonstat_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    int includemode, boozer_data *boozerdata, mhd_nonstat_data *mhddata,
    B_field_data *Bdata);

#endif
