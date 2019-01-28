/**
 * @file mhd.h
 * @brief Header file for mhd.c
 */
#ifndef MHD_H
#define MHD_H

#include "ascot5.h"
#include "error.h"

/**
 * @brief MHD parameters that will be offloaded to target
 */
typedef struct {
    int offload_array_length; /**< Number of elements in offload_array        */
} mhd_offload_data;

/**
 * @brief MHD parameters on the target
 */
typedef struct {
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
