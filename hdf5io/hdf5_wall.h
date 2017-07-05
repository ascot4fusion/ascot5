#ifndef HDF5_WALL_H
#define HDF5_WALL_H

#include <hdf5.h>
#include "../wall.h"


void hdf5_wall_init_offload(hid_t f, wall_offload_data* offload_data, real** offload_array);

#endif
