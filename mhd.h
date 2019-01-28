/**
 * @file mhd.h
 * @brief Header file for mhd.c
 */
#ifndef MHD_H
#define MHD_H

#include "ascot5.h"
#include "error.h"

#define MHD_MAXIMUM_NUMBER_OF_MODES 20

/**
 * @brief MHD parameters that will be offloaded to target
 */
typedef struct {

    int N_modes; /**< Number of modes                                         */
    int npsi;    /**< Number of psi grid points                               */
    int nmode[MHD_MAXIMUM_NUMBER_OF_MODES]; /**< Toroidal mode numbers        */
    int mmode[MHD_MAXIMUM_NUMBER_OF_MODES]; /**< Poloidal mode numbers        */

    int offload_array_length; /**< Number of elements in offload_array        */
} mhd_offload_data;

/**
 * @brief MHD parameters on the target
 */
typedef struct {
    int N_modes; /**< Number of modes                                         */
    int npsi;    /**< Number of psi grid points                               */
    int nmode[MHD_MAXIMUM_NUMBER_OF_MODES]; /**< Toroidal mode numbers        */
    int mmode[MHD_MAXIMUM_NUMBER_OF_MODES]; /**< Poloidal mode numbers        */

    real* alpha_nm;     /**< Radial profile of alpha_nm                       */
    real* phi_nm;       /**< Radial profile of phi_nm                         */
    real* amplitude_nm; /**< Amplitude of each mode                           */
    real* omega_nm;     /**< Toroidal rotation frequency of each mode [rad/s] */
} mhd_data;

int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array);

#pragma omp declare target
void mhd_init(mhd_data* MHDdata, mhd_offload_data* offload_data,
              real* offload_array);

#pragma omp declare simd uniform(MHDdata)
a5err mhd_eval(mhd_data* MHDdata);

#pragma omp end declare target

#endif
