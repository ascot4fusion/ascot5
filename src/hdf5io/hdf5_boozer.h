/**
 * @file hdf5_boozer.h
 * @brief Header file for hdf5_boozer.c
 */
#ifndef HDF5_BOOZER_H
#define HDF5_BOOZER_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../boozer.h"

int hdf5_boozer_init_offload(hid_t f, boozer_offload_data* offload_data,
                             real** offload_array, char* qid);

#endif
