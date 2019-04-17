/**
 * @file hdf5_options.h
 * @brief Header file for hdf5_options.c
 */
#ifndef HDF5_OPTIONS_H5
#define HDF5_OPTIONS_H5

#include "../simulate.h"
#include <hdf5.h>

int hdf5_options_read(hid_t f, sim_offload_data* sim, char* qid);

#endif
