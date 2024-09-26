/**
 * @file hdf5_boozer.h
 * @brief Header file for hdf5_boozer.c
 */
#ifndef HDF5_BOOZER_H
#define HDF5_BOOZER_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../boozer.h"

int hdf5_boozer_init(hid_t f, boozer_data* data, char* qid);

#endif
