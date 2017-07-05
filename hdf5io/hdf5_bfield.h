/**
 * @file hdf5_bfield.h
 * @brief Header file for hdf5_bfield.c
 */
#ifndef HDF5_BFIELD_H
#define HDF5_BFIELD_H
#include "../ascot5.h"
#include "../B_field.h"
#include "../B_2D.h"
#include "../B_3D.h"
#include "../B_GS.h"
#include "../B_ST.h"
#include "../B_TC.h"
#include "hdf5.h"

int hdf5_bfield_init_offload(hid_t f, B_field_offload_data* offload_data, real** offload_array);

void hdf5_bfield_init_offload_2D(hid_t f, B_2D_offload_data* offload_data, real** offload_array);

void hdf5_bfield_init_offload_3D(hid_t f, B_3D_offload_data* offload_data, real** offload_array);

void hdf5_bfield_init_offload_ST(hid_t f, B_ST_offload_data* offload_data, real** offload_array);

void hdf5_bfield_init_offload_TC(hid_t f, B_TC_offload_data* offload_data, real** offload_array);

#endif
