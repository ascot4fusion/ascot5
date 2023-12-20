#include "../simulate.h"
#include "./bmc_diag.h"
#include "../print.h"
#include "../physlib.h"

void buildDensityMatrixFromInputParticles(
    real **histogram,
    int dist_length,
    int n_particles,
    particle_state* input_particles,
    dist_5D_offload_data* dist5D,
    wall_data* wdata
);

real ISContinuousDistr(
    real z[5],
    real* probabilityMatrix,
    dist_5D_data* dist5D,
    wall_data* wallData
);

int init_simd_gc_particles(particle_state* ps, int n_ps, particle_simd_gc** p, B_field_data* Bdata);

void buildParticlesWeightsFromProbabilityMatrix(
    real* probabilityMatrix,
    real* weights,
    particle_state* ps,
    int n,
    dist_5D_data* dist,
    wall_data* wdata
);

int fmc_init_importance_sampling_from_source_distribution(
        int *n,
        particle_state** ps,
        diag_data* distr,
        int n_total,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array,
        offload_package* offload_data,
        int importanceSamplingProbability,
        int rk4_subcycles,
        particle_state* input_ps,
        int input_n_ps
    );

void buildISMatrixForMesh(
    int input_n_ps,
    particle_state* input_ps,
    real* ISMatrix,
    int importanceSamplingProbability,
    int dist_length,
    dist_5D_data* dist5D,
    wall_data* wallData
);

void buildISMatrixForParticles(
    diag_data* distr,
    int input_n_ps,
    real* Ekin,
    particle_state* input_ps,
    real* ISMatrix,
    int importanceSamplingProbability,
    int dist_length,
    dist_5D_data* dist5D,
    wall_data* wallData
);

int fmc_init_importance_sampling_mesh(
        int *n,
        particle_state** ps,
        int** ps_indexes,
        int n_total,
        int use_hermite,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array,
        offload_package* offload_data,
        int importanceSamplingProbability,
        int importanceSamplingdensity,
        int importanceSamplingFromParticles,
        real t,
        real m,
        real q,
        int rk4_subcycles,
        particle_state* input_ps,
        int input_n_ps
    );
int fmc_init_importance_sampling(
        int *n,
        particle_state** ps,
        int** ps_indexes,
        int n_total,
        int use_hermite,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array,
        real t,
        real m,
        real q,
        int rk4_subcycles
    );

void fmcInitImportanceSamplingMetropolis(
    int *n,
    particle_state** ps,
    diag_data* distr,
    int n_total,
    sim_offload_data* sim_offload,
    B_field_data* Bdata,
    real* offload_array,
    offload_package* offload_data,
    int importanceSamplingProbability,
    int rk4_subcycles,
    particle_state* input_ps,
    int input_n_ps,
    real t,
    real m,
    real q,
    real displacementPercentage
);

int bmc_init_particles(
        int mpi_rank,
        int mpi_size,
        int *n,
        int* n_tot,
        particle_state** ps,
        int** ps_indexes,
        int n_per_vertex,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array,
        real t,
        real m,
        real q,
        int rk4_subcycles
    );

void buildImportantSamplingHistogram(
    int dist_length,
    real *histogram,
    dist_5D_offload_data* dist5D,
    plasma_data* plasma_data,
    B_field_data* Bdata,
    int importanceSamplingProbability,
    int importanceSamplingdensity,
    int importanceSamplingFromInputParticles,
    real t,
    particle_state* ps,
    int n_ps,
    wall_data* wdata
);

void bmc_init_fo_particle(
        input_particle* p,
        real r, real phi, real z,
        real p_r, real p_phi, real p_z,
        real t,
        real m,
        real q,
        int id
    );
double r2();