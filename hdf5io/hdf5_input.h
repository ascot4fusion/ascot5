/**
 * @file hdf5_input.h
 * @brief Header file for hdf5_input.c
 */


#ifndef HDF5_INPUT_H5 
#define HDF5_INPUT_H5

#include "../simulate.h"

int hdf5_input(sim_offload_data* sim, 
	       real** B_offload_array,
	       real** E_offload_array,
	       real** plasma_offload_array,
	       real** wall_offload_array);

#endif
