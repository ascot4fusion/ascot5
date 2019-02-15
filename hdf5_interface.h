/**
 * @file hdf5_interface.h
 * @brief Header file for hdf5_interface.c
 */
#ifndef HDF5_INTERFACE_H5
#define HDF5_INTERFACE_H5

#include <hdf5.h>

#include "ascot5.h"
#include "simulate.h"

int hdf5_interface_read_input(sim_offload_data* sim,
                              real** B_offload_array,
                              real** E_offload_array,
                              real** plasma_offload_array,
                              real** neutral_offload_array,
                              real** wall_offload_array,
                              real** boozer_offload_array,
                              real** mhd_offload_array,
                              input_particle** p,
                              int* n_markers);

int hdf5_interface_init_results(sim_offload_data* sim, char* qid);

int hdf5_interface_write_state(char* fn, char* state, integer n,
                               particle_state* p);

int hdf5_interface_write_diagnostics(sim_offload_data* sim,
                                     real* diag_offload_array, char* out);

int hdf5_get_active_qid(hid_t f, const char* group, char qid[11]);
#endif
