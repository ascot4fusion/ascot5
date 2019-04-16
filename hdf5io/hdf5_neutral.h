/**
 * @file hdf5_neutral.h
 * @brief Header file for hdf5_neutral.c
 */
#ifndef HDF5_NEUTRAL_H
#define HDF5_NEUTRAL_H
#include "../ascot5.h"
#include "../neutral.h"
#include "../neutral/N0_3D.h"
#include "hdf5.h"

int hdf5_neutral_init_offload(hid_t f, neutral_offload_data* offload_data,
                              real** offload_array, char* qid);

#endif
