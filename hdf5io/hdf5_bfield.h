/**
 * @file hdf5_bfield.h
 * @brief Header file for hdf5_bfield.c
 */
#ifndef HDF5_BFIELD_H
#define HDF5_BFIELD_H
#include "../ascot5.h"
#include "../B_field.h"
#include "../Bfield/B_2DS.h"
#include "../Bfield/B_3DS.h"
#include "../Bfield/B_3DS_T.h"
#include "../Bfield/B_STS.h"
#include "../Bfield/B_GS.h"
#include "../Bfield/B_TC.h"
#include "hdf5.h"

int hdf5_bfield_init_offload(hid_t f, B_field_offload_data* offload_data, real** offload_array);

void hdf5_bfield_init_offload_2DS(hid_t f, B_2DS_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_bfield_init_offload_3DS(hid_t f, B_3DS_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_bfield_init_offload_STS(hid_t f, B_STS_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_bfield_init_offload_TC(hid_t f, B_TC_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_bfield_init_offload_GS(hid_t f, B_GS_offload_data* offload_data, real** offload_array, char* qid);

void hdf5_bfield_init_offload_3DS_T(hid_t f, B_3DS_T_offload_data* offload_data, real** offload_array, char* qid); 

#endif
