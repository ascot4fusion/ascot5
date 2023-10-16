/**
 * @file mhd_nonstat.h
 * @brief Header file for mhd_nonstat.c
 */
#ifndef MHD_NONSTAT_H
#define MHD_NONSTAT_H

#include "../ascot5.h"
#include "../error.h"
#include "../boozer.h"
#include "../spline/interp.h"
#include "../B_field.h"

/**
 * @brief MHD parameters that will be offloaded to target
 */
typedef struct {

    int n_modes;                          /**< Number of modes                */
    int nrho;                             /**< Number of rho grid points      */
    real rho_min;                         /**< rho grid minimum value         */
    real rho_max;                         /**< rho grid maximum value         */
    int ntime;                            /**< Number of time grid points     */
    real t_min;                           /**< time grid minimum value        */
    real t_max;                           /**< time grid maximum value        */
    int nmode[MHD_MODES_MAX_NUM];         /**< Toroidal mode numbers          */
    int mmode[MHD_MODES_MAX_NUM];         /**< Poloidal mode numbers          */
    real amplitude_nm[MHD_MODES_MAX_NUM]; /**< Amplitude of each mode         */
    real omega_nm[MHD_MODES_MAX_NUM];     /**< Toroidal rotation frequency of
                                               each mode [rad/s]              */
    real phase_nm[MHD_MODES_MAX_NUM];     /**< Phase of each mode [rad]       */

    int offload_array_length; /**< Number of elements in offload_array        */
} mhd_nonstat_offload_data;

/**
 * @brief MHD parameters on the target
 */
typedef struct {
    int n_modes;                          /**< Number of modes                */
    int nmode[MHD_MODES_MAX_NUM];         /**< Toroidal mode numbers          */
    int mmode[MHD_MODES_MAX_NUM];         /**< Poloidal mode numbers          */
    real amplitude_nm[MHD_MODES_MAX_NUM]; /**< Amplitude of each mode         */
    real omega_nm[MHD_MODES_MAX_NUM];     /**< Toroidal rotation frequency of
                                               each mode [rad/s]              */
    real phase_nm[MHD_MODES_MAX_NUM];     /**< Phase of each mode [rad]       */

    /**
     * @brief 2D splines (rho,time) for each mode's magnetic eigenfunction
     */
    interp2D_data alpha_nm[MHD_MODES_MAX_NUM];
    /**
     * @brief 2D splines (rho,time) for each mode's electric eigenfunction
     */
    interp2D_data phi_nm[MHD_MODES_MAX_NUM];
} mhd_nonstat_data;

int mhd_nonstat_init_offload(mhd_nonstat_offload_data* offload_data,
                             real** offload_array);

void mhd_nonstat_free_offload(mhd_nonstat_offload_data* offload_data,
                              real** offload_array);

void mhd_nonstat_init(mhd_nonstat_data* mhddata,
                      mhd_nonstat_offload_data* offload_data,
                      real* offload_array);
#pragma omp declare simd uniform(boozerdata, mhddata, Bdata, includemode)
a5err mhd_nonstat_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
                       int includemode, boozer_data* boozerdata,
                       mhd_nonstat_data* mhddata, B_field_data* Bdata);
#pragma omp declare simd uniform(boozerdata, mhddata, Bdata, pertonly,\
                                 includemode)
a5err mhd_nonstat_perturbations(real pert_field[7], real r, real phi, real z,
                                real t, int pertonly, int includemode,
                                boozer_data* boozerdata,
                                mhd_nonstat_data* mhddata, B_field_data* Bdata);


#endif
