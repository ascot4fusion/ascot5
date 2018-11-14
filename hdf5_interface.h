/**
 * @file hdf5_interface.h
 * @brief Header file for hdf5_interface.c
 */
#ifndef HDF5_INTERFACE_H5
#define HDF5_INTERFACE_H5

#include "simulate.h"

int hdf5_interface_read_input(sim_offload_data* sim,
                              real** B_offload_array,
                              real** E_offload_array,
                              real** plasma_offload_array,
                              real** neutral_offload_array,
                              real** wall_offload_array,
                              input_particle** p,
                              int* n_markers);

int hdf5_initoutput(sim_offload_data* sim, char* qid);
#endif
