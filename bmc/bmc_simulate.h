#include "../particle.h"
#include "../simulate.h"
#include "../endcond.h"
#include "../print.h"

void bmc_simulate_timestep_gc(int n_simd_particles, int n_coll_simd_particles, particle_simd_gc* p, particle_simd_gc* p_coll,
        int n_hermite_knots,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        real h, int n_rk4_subcycles
    );

void fmc_simulation(
        particle_state* ps,
        sim_offload_data* sim,
        offload_package* offload_data,
        real* offload_array,
        double* mic1_start, double* mic1_end,
        double* mic0_start, double* mic0_end,
        double* host_start, double* host_end,
        int n_mic,
        int n_host,
        real* diag_offload_array_host,
        real* diag_offload_array_mic0,
        real* diag_offload_array_mic1
    );

void init_particles_coll_simd_hermite(int n_simd_particles, int n_hermite_knots,
        particle_simd_gc* p_coll
    );
void copy_particles_simd_to_coll_simd(int n_simd_particles, int n_hermite_knots,
        particle_simd_gc* p, particle_simd_gc* p_coll
    );