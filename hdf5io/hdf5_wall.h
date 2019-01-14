/**
 * @file hdf5_wall.h
 * @brief Header file for hdf5_wall.c
 */
#ifndef HDF5_WALL_H
#define HDF5_WALL_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../wall.h"

int hdf5_wall_init_offload(hid_t f, wall_offload_data* offload_data,
                           real** offload_array, char* qid);

#endif
