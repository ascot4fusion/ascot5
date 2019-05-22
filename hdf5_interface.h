/**
 * @file hdf5_interface.h
 * @brief Header file for hdf5_interface.c
 */
#ifndef HDF5_INTERFACE_H5
#define HDF5_INTERFACE_H5

#include <hdf5.h>

#include "ascot5.h"
#include "simulate.h"

enum {
    hdf5_input_options = 0x1,
    hdf5_input_bfield  = 0x2,
    hdf5_input_efield  = 0x4,
    hdf5_input_plasma  = 0x8,
    hdf5_input_neutral = 0x10,
    hdf5_input_wall    = 0x20,
    hdf5_input_marker  = 0x40,
    hdf5_input_boozer  = 0x80,
    hdf5_input_mhd     = 0x100
};

int hdf5_interface_read_input(sim_offload_data* sim,
                              int input_active,
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
