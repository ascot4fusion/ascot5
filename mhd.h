/**
 * @file mhd.h
 * @brief Header file for mhd.c
 */
#ifndef MHD_H
#define MHD_H

#include "ascot5.h"
#include "error.h"
#include "boozer.h"
#include "spline/interp.h"

#define MHD_MODES_MAX_NUM 20

/**
 * @brief MHD parameters that will be offloaded to target
 */
typedef struct {

    int n_modes;   /**< Number of modes                                       */
    int npsi;      /**< Number of psi grid points                             */
    real psi_min;  /**< psi grid minimum value                                */
    real psi_max;  /**< psi grid maximum value                                */
    int ntime;     /**< Number of time grid points                            */
    real t_min;    /**< time grid minimum value                               */
    real t_max;    /**< time grid maximum value                               */
    int nmode[MHD_MODES_MAX_NUM];         /**< Toroidal mode numbers          */
    int mmode[MHD_MODES_MAX_NUM];         /**< Poloidal mode numbers          */
    real amplitude_nm[MHD_MODES_MAX_NUM]; /**< Amplitude of each mode         */
    real omega_nm[MHD_MODES_MAX_NUM];     /**< Toroidal rotation frequency of
                                               each mode [rad/s]              */

    int offload_array_length; /**< Number of elements in offload_array        */
} mhd_offload_data;

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

    /**< 2D splines (psi,time) for each mode's alpha_nm */
    interp2D_data alpha_nm[MHD_MODES_MAX_NUM];
    /**< 2D splines (psi,time) for each mode's phi_nm */
    interp2D_data phi_nm[MHD_MODES_MAX_NUM];
} mhd_data;

int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array);

void mhd_free_offload(mhd_offload_data* offload_data,
                      real** offload_array);

#pragma omp declare target
void mhd_init(mhd_data* MHDdata, mhd_offload_data* offload_data,
              real* offload_array);

#pragma omp declare simd uniform(boozerdata, mhddata)
a5err mhd_eval(real mhd_dmhd[10], real phase, real r, real phi, real z, real t,
               boozer_data* boozerdata, mhd_data* mhddata);

#pragma omp end declare target

#endif
