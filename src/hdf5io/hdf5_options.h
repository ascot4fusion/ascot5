/**
 * @file hdf5_options.h
 * @brief Header file for hdf5_options.c
 */
#ifndef HDF5_OPTIONS_H5
#define HDF5_OPTIONS_H5

#include "../simulate.h"
#include <hdf5.h>

int hdf5_options_read(hid_t f, sim_offload_data* sim, char* qid);

#define TOROIDAL_ANGLE_FILLER_VALUE 361 /**< Dummy value in poincare init */
#define POLOIDAL_ANGLE_FILLER_VALUE 361 /**< Dummy value in poincare init */
#define RADIAL_FILLER_VALUE 1000        /**< Dummy value in poincare init */

#endif
