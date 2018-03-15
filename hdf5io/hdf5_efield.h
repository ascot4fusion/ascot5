/**
 * @file hdf5_efield.h
 * @brief Header file for hdf5_efield.c
 */
#ifndef HDF5_EFIELD_H
#define HDF5_EFIELD_H

#include <hdf5.h>
#include "../E_field.h"
#include "../Efield/E_TC.h"
#include "../Efield/E_1D.h"
#include "../Efield/E_1DS.h"
#include "../Efield/E_3D.h"



int hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data, real** offload_array);

void hdf5_efield_init_offload_1D(hid_t f, E_1D_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_efield_init_offload_1DS(hid_t f, E_1DS_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_efield_init_offload_TC(hid_t f, E_TC_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_efield_init_offload_3D(hid_t f, E_3D_offload_data* offload_data, real** offload_array, char* qid);


#endif
