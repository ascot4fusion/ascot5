#include "../diag.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_5D.h"
#include "../simulate.h"
#include "../endcond.h"
#include "../consts.h"
#include "mpi.h"
#include <string.h>

int bmc_walltile_in_target(integer walltile);

void diag_move_distribution(sim_offload_data* sim, diag_data* diag_dest, diag_data* diag_src, int dist_length);

void bmc_update_particles_diag(
    int n_mpi_particles,
    particle_state* ps0,
    particle_state* ps1,
    int* ps_indexes,
    diag_data* diag0,
    diag_data* diag1,
    sim_data* sim,
    int n_montecarlo_steps
);

void bmc_diag_update_gc(
    particle_state* ps0,
    particle_state* ps1,
    int p0_index,
    diag_data* diag_0,
    diag_data* diag_1,
    int n_montecarlo_steps
);
void bmc_diag_update_fo(
    particle_state* ps0,
    particle_state* ps1,
    int p0_index,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps
);

void bmc_diag_5D_update_gc(
    dist_5D_data* dist1,
    dist_5D_data* dist0,
    int p0_index,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);
void bmc_diag_5D_update_fo(
    dist_5D_data* dist1,
    dist_5D_data* dist0,
    int p0_index,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);

void bmc_diag_6D_update_gc(
    dist_6D_data* dist1,
    dist_6D_data* dist0,
    int p0_index,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);
void bmc_diag_6D_update_fo(
    dist_6D_data* dist1,
    dist_6D_data* dist0,
    int p0_index,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);

int bmc_dist5D_gc_index(particle_state* ps, dist_5D_data* dist);
int bmc_dist5D_fo_index(particle_state* ps, dist_5D_data* dist);
int bmc_dist6D_gc_index(particle_state* ps, dist_6D_data* dist);
int bmc_dist6D_fo_index(particle_state* ps, dist_6D_data* dist);