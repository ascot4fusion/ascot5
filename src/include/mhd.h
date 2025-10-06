/**
 * @file mhd.h
 * @brief Header file for mhd.c
 *
 * Contains a list declaring all mhd_types, and declaration of
 * mhd_offload_data and mhd_data structs.
 */
#ifndef MHD_H
#define MHD_H

#include "defines.h"
#include "bfield.h"
#include "boozer.h"
#include "errors.h"
#include "interp.h"

/** @brief includemode parameter to include all modes (default) */
#define MHD_INCLUDE_ALL -1

/**
 * @brief MHD input types
 */
typedef enum mhd_type
{
    mhd_type_stat,   /**< MHD where mode amplitude does not depend on time */
    mhd_type_nonstat /**< MHD where mode amplitude depends on time         */
} mhd_type;

/**
 * @brief MHD stat parameters on the target
 */
typedef struct
{
    int n_modes;        /**< Number of modes            */
    real rho_min;       /**< rho grid minimum value     */
    real rho_max;       /**< rho grid maximum value     */
    int *nmode;         /**< Toroidal mode numbers      */
    int *mmode;         /**< Poloidal mode numbers      */
    real *amplitude_nm; /**< Amplitude of each mode     */
    real *omega_nm;     /**< Toroidal frequency [rad/s] */
    real *phase_nm;     /**< Phase of each mode [rad]   */

    /**
     * @brief 1D splines (rho) for each mode's magnetic eigenfunction
     */
    interp1D_data *alpha_nm;

    /**
     * @brief 1D splines (rho) for each mode's electric eigenfunction
     */
    interp1D_data *phi_nm;
} mhd_stat_data;

/**
 * @brief MHD parameters
 */
typedef struct
{
    int n_modes;        /**< Number of modes            */
    real rho_min;       /**< rho grid minimum value     */
    real rho_max;       /**< rho grid maximum value     */
    int *nmode;         /**< Toroidal mode numbers      */
    int *mmode;         /**< Poloidal mode numbers      */
    real *amplitude_nm; /**< Amplitude of each mode     */
    real *omega_nm;     /**< Toroidal frequency [rad/s] */
    real *phase_nm;     /**< Phase of each mode [rad]   */

    /**
     * @brief 2D splines (rho,time) for each mode's magnetic eigenfunction
     */
    interp2D_data *alpha_nm;
    /**
     * @brief 2D splines (rho,time) for each mode's electric eigenfunction
     */
    interp2D_data *phi_nm;
} mhd_nonstat_data;

/**
 * @brief MHD simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    mhd_stat_data *stat;       /**< Stat field or NULL if not active    */
    mhd_nonstat_data *nonstat; /**< Nonstat field or NULL if not active */
    mhd_type type;             /**< MHD type wrapped by this struct     */
} Mhd;

void mhd_free(Mhd *data);
void mhd_offload(Mhd *data);
DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, includemode)
a5err mhd_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, int includemode,
    Boozer *boozer, Mhd *mhd, Bfield *bfield);
DECLARE_TARGET_SIMD_UNIFORM(boozer, mhd, bfield, pertonly, includemode)
a5err mhd_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    int includemode, Boozer *boozer, Mhd *mhd, Bfield *bfield);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
int mhd_get_n_modes(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const int *mhd_get_nmode(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(Mhd)
const int *mhd_get_mmode(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const real *mhd_get_amplitude(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const real *mhd_get_frequency(Mhd *mhd);
DECLARE_TARGET_SIMD_UNIFORM(mhd)
const real *mhd_get_phase(Mhd *mhd);
#endif
