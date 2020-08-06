#include "../diag.h"
#include "../simulate.h"
#include "../endcond.h"
#include "../consts.h"

int bmc_walltile_in_target(integer walltile);

void bmc_diag_update_gc(
    particle_state* ps0,
    particle_state* ps1,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps
);
void bmc_diag_update_fo(
    particle_state* ps0,
    particle_state* ps1,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps
);

void bmc_diag_5D_update_gc(
    dist_5D_data* dist1,
    dist_5D_data* dist0,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);
void bmc_diag_5D_update_fo(
    dist_5D_data* dist1,
    dist_5D_data* dist0,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);

void bmc_diag_6D_update_gc(
    dist_6D_data* dist1,
    dist_6D_data* dist0,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);
void bmc_diag_6D_update_fo(
    dist_6D_data* dist1,
    dist_6D_data* dist0,
    particle_state* ps1,
    particle_state* ps0,
    int n_montecarlo_steps
);

int bmc_dist5D_gc_index(particle_state* ps, dist_5D_data* dist);
int bmc_dist5D_fo_index(particle_state* ps, dist_5D_data* dist);
int bmc_dist6D_gc_index(particle_state* ps, dist_6D_data* dist);
int bmc_dist6D_fo_index(particle_state* ps, dist_6D_data* dist);