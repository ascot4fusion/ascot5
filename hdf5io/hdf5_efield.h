/**
 * @file hdf5_efield.h
 * @brief Header file for hdf5_efielc.c
 */
#ifndef HDF5_EFIELD_H
#define HDF5_EFIELD_H

#include <hdf5.h>
#include "../E_field.h"
#include "../Efield/E_TC.h"
#include "../Efield/E_1DS.h"
#include "../Efield/E_3D.h"
#include "../Efield/E_3DS.h"

int hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data,
                             real** offload_array, char* qid);
#endif
