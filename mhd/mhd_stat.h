/**
 * @file mhd_stat.h
 * @brief Header file for mhd_stat.c
 */
#ifndef MHD_STAT_H
#define MHD_STAT_H

#include "../ascot5.h"
#include "../error.h"
#include "../boozer.h"
#include "../spline/interp.h"
#include "../B_field.h"

/**
 * @brief MHD stat parameters that will be offloaded to target
 */
typedef struct {

    int n_modes;                          /**< Number of modes                */
    int nrho;                             /**< Number of rho grid points      */
    real rho_min;                         /**< rho grid minimum value         */
    real rho_max;                         /**< rho grid maximum value         */
    int nmode[MHD_MODES_MAX_NUM];         /**< Toroidal mode numbers          */
    int mmode[MHD_MODES_MAX_NUM];         /**< Poloidal mode numbers          */
    real amplitude_nm[MHD_MODES_MAX_NUM]; /**< Amplitude of each mode         */
    real omega_nm[MHD_MODES_MAX_NUM];     /**< Toroidal rotation frequency of
                                               each mode [rad/s]              */
    real phase_nm[MHD_MODES_MAX_NUM];     /**< Phase of each mode [rad]       */

    int offload_array_length; /**< Number of elements in offload_array        */
} mhd_stat_offload_data;

/**
 * @brief MHD stat parameters on the target
 */
typedef struct {
    int n_modes;                          /**< Number of modes                */
    real rho_min;                         /**< rho grid minimum value         */
    real rho_max;                         /**< rho grid maximum value         */
    int nmode[MHD_MODES_MAX_NUM];         /**< Toroidal mode numbers          */
    int mmode[MHD_MODES_MAX_NUM];         /**< Poloidal mode numbers          */
    real amplitude_nm[MHD_MODES_MAX_NUM]; /**< Amplitude of each mode         */
    real omega_nm[MHD_MODES_MAX_NUM];     /**< Toroidal rotation frequency of
                                               each mode [rad/s]              */
    real phase_nm[MHD_MODES_MAX_NUM];     /**< Phase of each mode [rad]       */

    /**< 1D splines (rho) for each mode's magnetic eigenfunction */
    interp1D_data alpha_nm[MHD_MODES_MAX_NUM];
    /**< 1D splines (rho) for each mode's electric eigenfunction */
    interp1D_data phi_nm[MHD_MODES_MAX_NUM];
} mhd_stat_data;

int mhd_stat_init_offload(mhd_stat_offload_data* offload_data,
                          real** offload_array);

void mhd_stat_free_offload(mhd_stat_offload_data* offload_data,
                           real** offload_array);

#pragma omp declare target
void mhd_stat_init(mhd_stat_data* mhdata, mhd_stat_offload_data* offload_data,
                   real* offload_array);
#pragma omp declare simd uniform(boozerdata, mhddata)
a5err mhd_stat_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
                    boozer_data* boozerdata, mhd_stat_data* mhddata,
                    B_field_data* Bdata);
#pragma omp declare simd uniform(boozerdata, mhddata, Bdata)
a5err mhd_stat_perturbations(real pert_field[7], real r, real phi, real z,
                             real t, int pertonly, boozer_data* boozerdata,
                             mhd_stat_data* mhddata, B_field_data* Bdata);

#pragma omp end declare target

#endif
