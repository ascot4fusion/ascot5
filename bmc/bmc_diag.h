#include "../diag.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_5D.h"
#include "../simulate.h"
#include "../endcond.h"
#include "../consts.h"
#include "mpi.h"
#include "bmc_wall.h"
#include <string.h>

void diag_move_distribution(sim_offload_data* sim, diag_data* diag_dest, diag_data* diag_src, int* updated);

int bmc_update_distr5D(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int* p0_index,
        particle_simd_gc* p1,
        particle_simd_gc* p0,
        int n_simd_particles,
        wall_2d_data* w2d
    );

int fmc_update_distr5D_from_states(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int* p0_indexes,
        particle_state* p1,
        particle_state* p0,
        int n_particles,
        wall_2d_data* w2d
    );


void bmc_dist5D_state_indexes(particle_state* ps0, int* indexes, real* weights, int* target_hit, particle_state* ps, dist_5D_data* dist, wall_2d_data* w2d);
 void bmc_dist5D_gc_indexes(particle_simd_gc* p0, int* indexes, real* weights, int* target_hit, particle_simd_gc* p, int i, dist_5D_data* dist, wall_2d_data* w2d);
int bmc_dist6D_fo_index(particle_state* ps, dist_6D_data* dist);