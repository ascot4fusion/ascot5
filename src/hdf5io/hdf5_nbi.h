/**
 * @file hdf5_nbi.h
 * @brief Header file for hdf5_nbi.c
 */
#ifndef HDF5_NBI_H
#define HDF5_NBI_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../nbi.h"

int hdf5_nbi_init(hid_t f, nbi_data* data, char* qid);

#endif
