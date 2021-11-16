#ifndef ASCOT5_MAIN_H
#define ASCOT5_MAIN_H

#include "ascot5.h"

int offload(
		sim_offload_data *sim,
		real** B_offload_array,
		real** E_offload_array,
		real** plasma_offload_array,
		real** neutral_offload_array,
		real** wall_offload_array,
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
	    offload_package *offload_data,
		particle_state** ps,
	    real** diag_offload_array_mic0,
	    real** diag_offload_array_mic1,
	    real** diag_offload_array_host
);

#endif
