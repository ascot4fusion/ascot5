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
#include "simulate.h"

int prepare_markers(
    sim_data* sim, int n_tot, input_particle* pin, particle_state** pout,
    int* nprts);

int write_rungroup(sim_data* sim, particle_state* ps, int n_tot, char* qid);

int offload_and_simulate(
    sim_data* sim, int n_tot, int n_proc, particle_state* pin,
    int* n_gather, particle_state** pout);

int write_output(sim_data* sim, particle_state* ps_gathered, int n_tot);

void print_marker_summary(particle_state* ps, int n_tot);

#endif
