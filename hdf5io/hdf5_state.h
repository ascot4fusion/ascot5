/**
 * @file hdf5_particlestate.h
 * @brief Header file for hdf5_particlestate.c
 */
#ifndef HDF5_PARTICLESTATE
#define HDF5_PARTICLESTATE

#include <hdf5.h>
#include "../particle.h"

int hdf5_particlestate_write(hid_t f, char* qid, char *state, integer n,
                             particle_state* p);

#endif
