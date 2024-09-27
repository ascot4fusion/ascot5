/**
 * @file hdf5_asigma.h
 * @brief Header file for hdf5_asigma.c
 */
#ifndef HDF5_ASIGMA_H
#define HDF5_ASIGMA_H

#include <hdf5.h>
#include "../ascot5.h"
#include "../asigma.h"

int hdf5_asigma_init(hid_t f, asigma_data* data, char* qid);
#endif
