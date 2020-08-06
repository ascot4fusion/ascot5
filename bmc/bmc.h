#include "../diag.h"
#include "../simulate.h"
#include "bmc_diag.h"
#include "mpi.h"

void bmc_init_particles(
        int *n,
        input_particle **p,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        int n_montecarlo_steps
    );

void backward_monte_carlo(
        int n_montecarlo_steps,
        int n_tot_particles,
        int n_mpi_particles,
        input_particle* p_mpi,
        B_field_data* Bdata,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        double* mic1_start, double* mic1_end,
        double* mic0_start, double* mic0_end,
        double* host_start, double* host_end,
        int n_mic,
        int n_host
    );

void bmc_simulate_particles(
        particle_state* ps,
        int n_tot_particles,
        sim_offload_data* sim,
        offload_package* offload_data,
        real* offload_array,
        double* mic1_start, double* mic1_end,
        double* mic0_start, double* mic0_end,
        double* host_start, double* host_end,
        int n_mic,
        int n_host);

void mpi_reduce_distribution(sim_offload_data* sim, diag_data* diag0, diag_data* diag1, int dist_length);