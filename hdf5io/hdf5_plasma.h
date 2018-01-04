/**
 * @file hdf5_plasma.h
 * @brief Header file for hdf5_plasma.c
 */
#ifndef HDF5_PLASMA_H
#define HDF5_PLASMA_H
#include "../ascot5.h"
#include "../plasma_1d.h"
#include "../plasma_1DS.h"
#include "hdf5.h"

int hdf5_plasma_init_offload(hid_t f, plasma_1d_offload_data* offload_data,
                            real** offload_array);

int hdf5_plasma_init_offload_1DS(hid_t f, plasma_1DS_offload_data* offload_data,
                            real** offload_array);

#endif
