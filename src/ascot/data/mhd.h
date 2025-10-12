/**
 * @file mhd.h
 * @brief Header file for mhd.c
 *
 * Contains a list declaring all mhd_types, and declaration of
 * mhd_offload_data and mhd_data structs.
 */
#ifndef MHD_H
#define MHD_H

#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd_data.h"

/** @brief includemode parameter to include all modes (default) */
#define MHD_INCLUDE_ALL 0

void Mhd_free(Mhd *data);
void Mhd_offload(Mhd *data);
DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, includemode)
err_t Mhd_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, size_t includemode,
    Boozer *boozer, Mhd *mhd, Bfield *bfield);
DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, pertonly, includemode)
err_t Mhd_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    size_t includemode, Boozer *boozer, Mhd *mhd, Bfield *bfield);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
size_t Mhd_get_n_modes(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const int *Mhd_get_nmode(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(Mhd)
const int *Mhd_get_mmode(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const real *Mhd_get_amplitude(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const real *Mhd_get_frequency(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const real *Mhd_get_phase(Mhd *mhd);
#endif
