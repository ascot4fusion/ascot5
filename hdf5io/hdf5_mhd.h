/**
 * @file hdf5_mhd.h
 * @brief Header file for hdf5_mhd.c
 */
#ifndef HDF5_MHD_H
#define HDF5_MHD_H
#include "../ascot5.h"
#include "../mhd.h"
#include "hdf5.h"

int hdf5_mhd_init_offload(hid_t f, mhd_offload_data* offload_data,
                          real** offload_array, char* qid);

#endif
