/**
 * @file hdf5_simulate.h
 * @brief Header file for hdf5_simulate.c
 */


#ifndef HDF5_SIMULATE_H5 
#define HDF5_SIMULATE_H5

#include "simulate.h"
#include "hdf5.h"

void hdf5_simulate(hid_t f, sim_offload_data* sim);

#endif
