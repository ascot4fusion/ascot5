/**
 * @file hdf5_marker.h
 * @brief Header file for hdf5_marker.c
 */
#ifndef HDF5_MARKER_H
#define HDF5_MARKER_H

#include <hdf5.h>
#include "../particle.h"

int hdf5_marker_read(hid_t f, int *n, input_particle** p, char* qid);
int hdf5_marker_write_particle(hid_t f, int n, input_particle* p, char* qid);
int hdf5_marker_write_particle_shined(hid_t f, int n, input_particle* p, char* qid);

#endif
