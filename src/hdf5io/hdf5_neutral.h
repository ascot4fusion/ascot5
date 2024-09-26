/**
 * @file hdf5_neutral.h
 * @brief Header file for hdf5_neutral.c
 */
#ifndef HDF5_NEUTRAL_H
#define HDF5_NEUTRAL_H
#include "../ascot5.h"
#include "../neutral.h"
#include "hdf5.h"

int hdf5_neutral_init(hid_t f, neutral_data* data, char* qid);
#endif
