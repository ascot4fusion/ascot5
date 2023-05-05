/**
 * @file ascot5_main.h
 * @brief Functions to execute main program externally
 *
 * The main progam for ascot5 is in ascot5_main.c. This header file declares
 * the functions needed to carry out the simulation, so that hey can be used by
 * an external program (i.e. via the python interface).
 */
#ifndef ASCOT5_MAIN_H
#define ASCOT5_MAIN_H

#include "ascot5.h"
#include "offload.h"
#include "simulate.h"

int prepare_offload(
    sim_offload_data* sim, int mpi_size, int mpi_rank, int n_tot,
    offload_package* offload_data, real** B_offload_array,
    real** E_offload_array, real** plasma_offload_array,
    real** neutral_offload_array, real** wall_offload_array,
    int** wall_int_offload_array, real** boozer_offload_array,
    real** mhd_offload_array, input_particle** p, particle_state** ps,
    int* nprts, real** offload_array, int** int_offload_array,
    real** diag_offload_array);

int write_rungroup(
    sim_offload_data* sim, int mpi_size, int mpi_rank, int mpi_root,
    int n_tot, particle_state* ps, char* qid);

int offload_and_simulate(
    sim_offload_data* sim, int mpi_rank, int nprts, particle_state* ps,
    offload_package* offload_data, real* offload_array,
    int* int_offload_array, real* diag_offload_array);

int gather_output(
    sim_offload_data* sim, int mpi_size, int mpi_rank, int mpi_root,
    int n_tot, particle_state* ps, particle_state** ps_gathered,
    real* diag_offload_array);

int write_output(
    sim_offload_data* sim, int mpi_rank, int mpi_root,
    particle_state* ps_gathered, int n_tot,
    real* diag_offload_array);

int free_simulation_data(
    sim_offload_data* sim, particle_state** ps, offload_package* offload_data,
    real** offload_array, int** int_offload_array, real** diag_offload_array);

int free_particlestate(particle_state* ps);

void print_marker_summary(particle_state* ps, int n_tot);

#endif
