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
                           real** offload_array, int** int_offload_array,
                           char* qid);

int hdf5_wall_2d_to_offload(
    wall_2d_offload_data *offload_data, real **offload_array,
    int nelements, real *r, real *z);

int hdf5_wall_3d_to_offload(
    wall_3d_offload_data *offload_data, real **offload_array,
    int nelements, real* x1x2x3, real* y1y2y3, real* z1z2z3);
#endif
