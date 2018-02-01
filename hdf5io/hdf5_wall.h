#ifndef HDF5_WALL_H
#define HDF5_WALL_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../wall.h"
#include "../wall/wall_2d.h"
#include "../wall/wall_3d.h"

int hdf5_wall_init_offload(hid_t f, wall_offload_data* offload_data, real** offload_array);

void hdf5_wall_init_offload_2D(hid_t f, wall_2d_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_wall_init_offload_3D(hid_t f, wall_3d_offload_data* offload_data, real** offload_array, char* qid);

#endif
