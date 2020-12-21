#include "../diag.h"
#include "../simulate.h"
#include "bmc_diag.h"
#include "mpi.h"
#include "../print.h"
#include "../hdf5_interface.h"
#include "../math.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_5D.h"
#include "bmc_simulate.h"
#include "../simulate/mccc/mccc_coefs.h"
#include "../simulate/mccc/mccc.h"
#include "../physlib.h"

void bmc_setup_endconds(sim_offload_data* sim, real timestep);

int backward_monte_carlo(
        int n_tot_particles,
        int n_mpi_particles,
        int n_hermite_knots,
        particle_state* ps,
        int* ps_indexes,
        B_field_data* Bdata,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        int mpi_rank,
        real t1,
        real t0,
        real h,
        int rk4_subcycles
    );

void backward_monte_carlo_gc(
        particle_state* ps,
        int* p0_indexes,
        int n_mpi_particles,
        int n_hermite_knots,
        sim_offload_data* sim_offload,
        sim_data* sim,
        offload_package* offload_data,
        real* offload_array,
        diag_data* distr0,
        diag_data* distr1,
        B_field_data* Bdata,
        real t1,
        real t0,
        real h,
        int rk4_subcycles
    );

int forward_monte_carlo(
        int n_tot_particles,
        int n_mpi_particles,
        int n_montecarlo_steps,
        particle_state* ps1,
        int* ps1_indexes,
        B_field_data* Bdata,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        double* mic1_start, double* mic1_end,
        double* mic0_start, double* mic0_end,
        double* host_start, double* host_end,
        int n_mic,
        int n_host,
        int mpi_rank,
        bool importance_sampling,
        real t1,
        real t0
    );