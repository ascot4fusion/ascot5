/**
 * @file mhd.h
 * @brief Header file for mhd.c
 *
 * Contains a list declaring all mhd_types, and declaration of
 * mhd_offload_data and mhd_data structs.
 */
#ifndef MHD_H
#define MHD_H

#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "boozer.h"
#include "mhd/mhd_stat.h"
#include "mhd/mhd_nonstat.h"

/** @brief includemode parameter to include all modes (default) */
#define MHD_INCLUDE_ALL -1

/**
 * @brief MHD input types
 */
typedef enum mhd_type {
    mhd_type_stat,   /**< MHD where mode amplitude does not depend on time */
    mhd_type_nonstat /**< MHD where mode amplitude depends on time         */
} mhd_type;

/**
 * @brief MHD simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct {
    mhd_type type;            /**< MHD type wrapped by this struct     */
    mhd_stat_data stat;       /**< Stat field or NULL if not active    */
    mhd_nonstat_data nonstat; /**< Nonstat field or NULL if not active */
} mhd_data;

void mhd_free(mhd_data* data);
void mhd_offload(mhd_data* data);
DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, includemode)
a5err mhd_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
               int includemode, boozer_data* boozerdata, mhd_data* mhddata,
               B_field_data* Bdata);
DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, pertonly,\
                                 includemode)
a5err mhd_perturbations(real pert_field[7], real r, real phi, real z,
                        real t, int pertonly, int includemode,
                        boozer_data* boozerdata, mhd_data* mhddata,
                        B_field_data* Bdata);
DECLARE_TARGET_SIMD_UNIFORM(mhddata)
int mhd_get_n_modes(mhd_data* mhddata);
DECLARE_TARGET_SIMD_UNIFORM(mhddata)
const int* mhd_get_nmode(mhd_data* mhddata);
DECLARE_TARGET_SIMD_UNIFORM(mhddata)
const int* mhd_get_mmode(mhd_data* mhddata);
DECLARE_TARGET_SIMD_UNIFORM(mhddata)
const real* mhd_get_amplitude(mhd_data* mhddata);
DECLARE_TARGET_SIMD_UNIFORM(mhddata)
const real* mhd_get_frequency(mhd_data* mhddata);
DECLARE_TARGET_SIMD_UNIFORM(mhddata)
const real*  mhd_get_phase(mhd_data* mhddata);
#endif
