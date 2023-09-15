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
#include "../Bfield/B_STS.h"
#include "../Bfield/B_GS.h"
#include "../Bfield/B_TC.h"
#include "hdf5.h"

int hdf5_bfield_init_offload(hid_t f, B_field_offload_data* offload_data,
                             real** offload_array, char* qid);

#endif
