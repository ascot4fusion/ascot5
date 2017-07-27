/**
 * @file hdf5_histogram.h
 * @brief Header file for hdf5_histogram.c
 */
#ifndef HDF5_PARTICLESTATE
#define HDF5_PARTICLESTATE

#include <hdf5.h>
#include "../particle.h"

int hdf5_particlestate_write(char* fn, char *state, int n, particle_state* p);

#endif
