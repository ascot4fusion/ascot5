#include "../ascot5.h"
#include "../simulate.h"
#include "hdf5_dist.h"

void hdf5_diag_write(sim_offload_data* sim, real* diag_offload_array, char* out, char* qid) {
    if(sim->diag_offload_data.dist5D_collect) {
        hdf5_dist_write_5D(&sim->diag_offload_data.dist5D, &diag_offload_array[sim->diag_offload_data.offload_dist5D_index], out, qid);
    }

    if(sim->diag_offload_data.dist6D_collect) {
        hdf5_dist_write_6D(&sim->diag_offload_data.dist6D, &diag_offload_array[sim->diag_offload_data.offload_dist6D_index], out, qid);
    }
    if(sim->diag_offload_data.distrho5D_collect) {
        hdf5_dist_write_rho5D(&sim->diag_offload_data.distrho5D, &diag_offload_array[sim->diag_offload_data.offload_distrho5D_index], out, qid);
    }

    if(sim->diag_offload_data.distrho6D_collect) {
        hdf5_dist_write_rho6D(&sim->diag_offload_data.distrho6D, &diag_offload_array[sim->diag_offload_data.offload_distrho6D_index], out, qid);
    }
}
