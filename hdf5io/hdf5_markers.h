/**
 * @file hdf5_markers.h
 * @brief Header file for hdf5_markers.c
 */
#ifndef HDF5_MARKERS_H
#define HDF5_MARKERS_H

#include <hdf5.h>
#include "../particle.h"

int hdf5_markers_read(hid_t f, int *n, input_particle** p, char* qid);

#endif
