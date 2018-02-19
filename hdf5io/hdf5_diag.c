#include "../ascot5.h"
#include "../simulate.h"
#include "hdf5_dist.h"

void hdf5_diag_write(sim_offload_data* sim, real* diag_offload_array, char* out, char* qid) {
    if(sim->diag_offload_data.dist5D_collect) {
        hdf5_dist_write_5D(&sim->diag_offload_data.dist5D, diag_offload_array, out, qid);
    }
}

