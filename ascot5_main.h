#ifndef ASCOT5_MAIN_H
#define ASCOT5_MAIN_H

#include "ascot5.h"
#include "offload.h"
#include "simulate.h"


int offload(
		sim_offload_data *sim,
		real** B_offload_array,
		real** E_offload_array,
		real** plasma_offload_array,
		real** neutral_offload_array,
		real** wall_offload_array,
		int**  wall_int_offload_array,
		real** boozer_offload_array,
		real** mhd_offload_array,
		int n_tot,
		int mpi_rank,
		int mpi_size,
		int mpi_root,
		char *qid,
		int *nprts,
		input_particle **p,
	    int* n_gathered,
	    real **offload_array,
		int  **int_offload_array,
	    offload_package *offload_data,
		particle_state** ps,
	    real** diag_offload_array
);

int free_ps(particle_state *ps);

int cleanup( sim_offload_data sim,    particle_state* ps,     particle_state* ps_gathered,
	    real** diag_offload_array,
	    real* offload_array,
	    offload_package *offload_data
		);

int run(
		int nprts,
		int mpi_rank,
		particle_state *ps,
	    real *offload_array,
	    int* int_offload_array,
	    real *diag_offload_array,
		sim_offload_data *sim,
	    offload_package *offload_data
		);

int gather_output(particle_state *ps, particle_state **ps_gathered,
	    int *n_gathered, int n_tot, int mpi_rank, int mpi_size, int mpi_root,
		sim_offload_data sim, real* diag_offload_array);

int write_output(sim_offload_data sim, int mpi_rank, int mpi_root,
		particle_state *ps_gathered, int n_gathered,
		real* diag_offload_array);


void marker_summary(particle_state* ps, int n);

#endif

