#ifndef HDF5_EFIELD_H
#define HDF5_EFIELD_H

#include <hdf5.h>
#include "../E_field.h"


void hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data, real** offload_array);

void hdf5_efield_init_offload_1D(hid_t f, E_1D_offload_data* offload_data, real** offload_array);

#endif
