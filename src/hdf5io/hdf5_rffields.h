/**
 * @file hdf5_rffields.h
 * @brief Header file for rffields.c
 */
#ifndef HDF5_RF_FIELD_FO_H
#define HDF5_RF_FIELD_FO_H

#include <hdf5.h>
#include "../rf_fields_fo.h"

int hdf5_rffields_init(hid_t f, RF2D_fields* data, char* qid);
#endif
