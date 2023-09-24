/**
 * @file hdf5_state.h
 * @brief Header file for hdf5_state.c
 */
#ifndef HDF5_STATE
#define HDF5_STATE

#include <hdf5.h>
#include "../particle.h"

int hdf5_state_write(hid_t f, char* run, char *state, integer n,
                     particle_state* p);

#endif
