
#include "../ascot5.h"
#include "../simulate.h"
#include "hdf5_dist.h"

void hdf5_diag_write(sim_offload_data* sim, real* diag_offload_array) {
    
    if(sim->diag_offload_data.dist4D_collect) {
	hdf5_dist_write_rzvv(&sim->diag_offload_data.dist4D, diag_offload_array, sim->hdf5fn);
    }

}

