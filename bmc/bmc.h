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

void bmc_setup_endconds(sim_offload_data* sim);

int bmc_init_particles(
        int *n,
        particle_state** ps,
        int** ps_indexes,
        int n_per_vertex,
        int use_hermite,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array
    );

int backward_monte_carlo(
        int n_tot_particles,
        int n_mpi_particles,
        particle_state* ps,
        int* ps_indexes,
        B_field_data* Bdata,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        int mpi_rank
    );

void backward_monte_carlo_gc(
        particle_state* ps,
        int* p0_indexes,
        int n_mpi_particles,
        sim_offload_data* sim_offload,
        sim_data* sim,
        offload_package* offload_data,
        real* offload_array,
        diag_data* distr0,
        diag_data* distr1,
        B_field_data* Bdata
    );

void write_probability_distribution(
    sim_offload_data* sim_offload,
    diag_data* distr,
    real* distr_array,
    int mpi_rank
);