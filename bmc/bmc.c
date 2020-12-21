#include "bmc.h"

#define TIMESTEP 1E-7 // TODO use input HDF
#define T0 9E-7
#define T1 1E-6
#define MASS 3.3452438E-27
#define CHARGE 1.60217662E-19
#define RK4_SUBCYCLES 5
#define PI2E0_5 2.50662827463

double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}

// Assumes 0 <= max <= RAND_MAX
// Returns in the closed interval [0, max]
long random_at_most(long max) {
  unsigned long
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    num_bins = (unsigned long) max + 1,
    num_rand = (unsigned long) RAND_MAX + 1,
    bin_size = num_rand / num_bins,
    defect   = num_rand % num_bins;

  long x;
  do {
   x = random();
  }
  // This is carefully written not to overflow
  while (num_rand - defect <= (unsigned long)x);

  // Truncated division is intentional
  return x/bin_size;
}

void bmc_setup_endconds(sim_offload_data* sim) {
    sim->endcond_active = endcond_tmax | endcond_wall;
    sim->fix_usrdef_val = TIMESTEP;
    sim->fix_usrdef_use = 1;

    // TODO: check if this is not necessary anymore
    // force n_time and n_q to be 1 in sim struct
    diag_offload_data* diag_data = &sim->diag_offload_data;
    if (diag_data->dist6D_collect) {
        diag_data->dist6D.n_time = 1;
        diag_data->dist6D.n_q = 1;
    }
    if (diag_data->dist5D_collect) {
        diag_data->dist5D.n_time = 1;
        diag_data->dist5D.n_q = 1;
    }
}

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
        int mpi_rank
    ) {

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nStarting Backward Monte Carlo. N particles: %d.\n", n_mpi_particles);

    /* Allow threads to spawn threads */
    omp_set_nested(1);

    // debug print of first particles
    for(int i = 0; i < 10; i++) {
        printf("Particle %d %f %f %f %f %f %f %f\n", i, ps[i].r, ps[i].phi, ps[i].z, ps[i].ppar, ps[i].rho, ps[i].rprt, ps[i].p_r);
    }

    // initialize distributions
    diag_data distr0, distr1;
    real* distr0_array, *distr1_array;
    diag_init_offload(&sim_offload->diag_offload_data, &distr0_array, 1);
    diag_init(&distr0, &sim_offload->diag_offload_data, distr0_array);
    diag_init_offload(&sim_offload->diag_offload_data, &distr1_array, 1);
    diag_init(&distr1, &sim_offload->diag_offload_data, distr1_array);

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);
    real* ptr = offload_array + 
            sim_offload->B_offload_data.offload_array_length +
            sim_offload->E_offload_data.offload_array_length +
            sim_offload->plasma_offload_data.offload_array_length +
            sim_offload->neutral_offload_data.offload_array_length;
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    if (distr0.dist5D_collect) {
        backward_monte_carlo_gc(ps, ps_indexes, n_mpi_particles, n_hermite_knots, sim_offload, &sim, offload_data, offload_array,
                                &distr0, &distr1, Bdata);
    } else {
        // TODO: FULL ORBIT
    }

    write_probability_distribution(sim_offload, &distr0, distr0_array, mpi_rank, true);

    // Free sitribution data
    diag_free_offload(&sim_offload->diag_offload_data, &distr1_array);
    diag_free_offload(&sim_offload->diag_offload_data, &distr0_array);

    return 0;

}

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
        B_field_data* Bdata
    ) {

    // int n_simd_particles = n_mpi_particles / NSIMD;
    particle_simd_gc *p0, *p1, *pcoll1, *pcoll0;
    int n_simd_particles = init_simd_gc_particles(ps, n_mpi_particles, &p0, Bdata);
    init_simd_gc_particles(ps, n_mpi_particles, &p1, Bdata);

    // compute the number of collisional particles needed.
    // = n_particles * hermite_knots
    int n_coll_simd = (n_mpi_particles * n_hermite_knots + NSIMD - 1) / NSIMD;
    pcoll1 = (particle_simd_gc *)malloc(n_coll_simd * sizeof(particle_simd_gc));  
    pcoll0 = (particle_simd_gc *)malloc(n_coll_simd * sizeof(particle_simd_gc));  
    int *p0_indexes_coll = (int *)malloc(n_mpi_particles * n_hermite_knots * sizeof(int));  
    for (int i=0; i < n_mpi_particles; i++) {
        for (int j=0; j < n_hermite_knots; j++) {
            p0_indexes_coll[i * n_hermite_knots + j] = p0_indexes[i];
        }
    }

    // init Hermite knots and weights for the collisional SIMD array
    init_particles_coll_simd_hermite(n_simd_particles, n_hermite_knots, pcoll1);
    init_particles_coll_simd_hermite(n_simd_particles, n_hermite_knots, pcoll0);
    copy_particles_simd_to_coll_simd(n_simd_particles, n_hermite_knots, p0, pcoll0);

    for (double t=T1; t >= T0; t -= TIMESTEP) {

        // reset particle initial states to ps1 and pcoll1
        memcpy(p1, p0, n_simd_particles * sizeof(particle_simd_gc));

        // set time in particle states
        for(int i = 0; i < n_simd_particles; i++) {
            for(int j = 0; j < NSIMD; j++) {
                if (p1[i].id[j]) {
                    p1[i].time[j] = t;
                }
            }
        }

        #ifdef TARGET
            int n_mic = n_simd_particles / TARGET;
            int n_mic_coll = n_coll_simd / TARGET;
        #else
            int n_mic = 0;
            int n_mic_coll = 0;
        #endif

        // simulate one step of all needed particles.
        // Split particles between offloading cores
        #pragma omp parallel sections num_threads(3)
        {
            #if TARGET >= 1
                #pragma omp section
                {
                    #pragma omp target device(0) map( \
                        p1[0:n_mic], \
                        pcoll1[0:n_mic_coll], \
                        offload_array[0:offload_data.offload_array_length] \
                    )
                    bmc_simulate_timestep_gc(n_mic, n_mic_coll, p1, pcoll1, n_hermite_knots, sim_offload, offload_data, offload_array, TIMESTEP, RK4_SUBCYCLES);
                }
            #endif
            #ifdef TARGET >= 2
                #pragma omp section
                {
                    #pragma omp target device(1) map( \
                        p1[n_mic:2*n_mic], \
                        pcoll1[n_mic_coll:2*n_mic_coll], \
                        offload_array[0:offload_data.offload_array_length] \
                    )
                    bmc_simulate_timestep_gc(n_mic, n_mic_coll, p1 + n_mic, pcoll1 + n_mic_coll, n_hermite_knots, sim_offload, offload_data, offload_array, TIMESTEP, RK4_SUBCYCLES);
                }
            #endif
            #ifndef TARGET
                #pragma omp section
                {
                    bmc_simulate_timestep_gc(n_simd_particles, n_coll_simd, p1, pcoll1, n_hermite_knots, sim_offload, offload_data, offload_array, TIMESTEP, RK4_SUBCYCLES);
                }
            #endif
        }

        // // Update the probability distribution
        int n_updated = bmc_update_distr5D(&distr1->dist5D, &distr0->dist5D, p0_indexes_coll, pcoll1, pcoll0, n_coll_simd, &(sim->wall_data.w2d));
        

        // // shift distributions
        diag_move_distribution(sim_offload, distr0, distr1, &n_updated);
        printf("Time %f Updated %d\n", t, n_updated);
        printf("Distribution moved\n");

        // reset particle initial states to ps1
        memcpy(p1, p0, n_simd_particles * sizeof(particle_simd_gc));
    } 
}

int init_simd_gc_particles(particle_state* ps, int n_ps, particle_simd_gc** p, B_field_data* Bdata) {
    int n_simd = (n_ps + NSIMD - 1) / NSIMD;
    *p = (particle_simd_gc *)malloc(n_simd * sizeof(particle_simd_gc)); 
    int i_ps = 0, i_simd = 0;
    while (i_ps < n_ps) {
        for (int j = 0; j < NSIMD; j++) {
            if (i_ps < n_ps) {
                particle_state_to_gc(ps + i_ps, i_ps, *p + i_simd, j, Bdata);
                i_ps++;
            } else {
                p[0][i_simd].running[j] = 0;
                p[0][i_simd].id[j] = -1;
            }
        }
        i_simd++;
    }
    return n_simd;
}

int bmc_init_particles(
        int *n,
        particle_state** ps,
        int** ps_indexes,
        int n_per_vertex,
        int use_hermite,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array
    ) {
    // vacate the phase space to find the phase-space points in the mesh
    // suitable for simulating the BMC scheme

    // TODO: use_hermite and n_per_vertex are deprecated.
    // Hermite weights and n particles per vertex are due to be computed when creating SIMD particles

    // init hermite params (N=5)
    real hermiteK[5] = {-2.856970, -1.355626, 0.000000, 1.355626, 2.856970};
    real hermiteW[5] = {0.028218, 0.556662, 1.336868, 0.556662, 0.028218};

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);
    
    real* ptr = offload_array + 
            sim_offload->B_offload_data.offload_array_length +
            sim_offload->E_offload_data.offload_array_length +
            sim_offload->plasma_offload_data.offload_array_length +
            sim_offload->neutral_offload_data.offload_array_length;
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    real r;
    real phi;
    real z;
    real ppara;
    real pperp;
    real pr;
    real pphi;
    real pz;

    int n_r, n_phi, n_z, n_ppara, n_pperp, n_pr, n_pphi, n_pz;
    real max_r, max_phi, max_z, max_ppara, max_pperp, max_pr, max_pphi, max_pz;
    real min_r, min_phi, min_z, min_ppara, min_pperp, min_pr, min_pphi, min_pz;
    if (sim_offload->diag_offload_data.dist5D_collect) {
        dist_5D_offload_data dist5D = sim_offload->diag_offload_data.dist5D;
        n_r = dist5D.n_r;
        n_phi = dist5D.n_phi;
        n_z = dist5D.n_z;
        n_ppara = dist5D.n_ppara;
        n_pperp = dist5D.n_pperp;
        max_r = dist5D.max_r;
        max_phi = dist5D.max_phi;
        max_z = dist5D.max_z;
        max_ppara = dist5D.max_ppara;
        max_pperp = dist5D.max_pperp;
        min_r = dist5D.min_r;
        min_phi = dist5D.min_phi;
        min_z = dist5D.min_z;
        min_ppara = dist5D.min_ppara;
        min_pperp = dist5D.min_pperp;
    } else {
        dist_6D_offload_data dist6D = sim_offload->diag_offload_data.dist6D;
        n_r = dist6D.n_r;
        n_phi = dist6D.n_phi;
        n_z = dist6D.n_z;
        n_pr = dist6D.n_pr;
        n_pphi = dist6D.n_pphi;
        n_pz = dist6D.n_pz;
        max_r = dist6D.max_r;
        max_phi = dist6D.max_phi;
        max_z = dist6D.max_z;
        max_pr = dist6D.max_pr;
        max_pphi = dist6D.max_pphi;
        max_pz = dist6D.max_pz;
        min_r = dist6D.min_r;
        min_phi = dist6D.min_phi;
        min_z = dist6D.min_z;
        min_z = dist6D.min_z;
        min_pr = dist6D.min_pr;
        min_pphi = dist6D.min_pphi;
    }

    if (sim_offload->diag_offload_data.dist5D_collect) {
        printf("5D: n_r %d n_phi %d n_z %d n_vpara %d n_vperp %d n_per_vertex %d\n", n_r, n_phi, n_z, n_ppara, n_pperp, n_per_vertex);
        *n = n_r * n_phi * n_z * (n_ppara - 2) * (n_pperp - 2) * n_per_vertex;
    } else {
        printf("6D: n_r %d n_phi %d n_z %d n_vr %d n_vphi %d n_vz %d n_per_vertex %d\n", n_r, n_phi, n_z, n_pr, n_pphi, n_pz, n_per_vertex);
        *n = n_r * n_phi * n_z * n_pr * n_pphi * n_pz * n_per_vertex;
    }

    print_out(VERBOSE_NORMAL, "Mesh size %d.\n", *n / n_per_vertex);

    *ps = (particle_state *)malloc(*n * sizeof(particle_state));
    *ps_indexes = (int *)malloc(*n * sizeof(int));
    input_particle p_tmp; // tmp particle
    particle_state ps_tmp; // tmp particle

    int i = 0;
    for (int i_r = 0; i_r < n_r; ++i_r) {
        r = (max_r - min_r) * i_r / n_r + min_r;
        for (int i_phi = 0; i_phi < n_phi; ++i_phi) {
            phi = (max_phi - min_phi) * i_phi / n_phi + min_phi;
            for (int i_z = 0; i_z < n_z; ++i_z) {
                z = (max_z - min_z) * i_z / n_z + min_z;

                if (!wall_2d_inside(r, z, &sim.wall_data.w2d)) {
                    continue;
                }

                if (sim_offload->diag_offload_data.dist5D_collect) {                
                    for (int i_ppara = 1; i_ppara < n_ppara - 1; ++i_ppara) {
                        ppara = (max_ppara - min_ppara) * i_ppara / n_ppara + min_ppara;
                        for (int i_pperp = 1; i_pperp < n_pperp - 1; ++i_pperp) {
                            pperp = (max_pperp - min_pperp) * i_pperp / n_pperp + min_pperp;
                            bmc_5D_to_particle_state(Bdata, r, phi, z, ppara, pperp, T1, i, &ps_tmp);

                            unsigned long index = dist_5D_index(i_r, i_phi, i_z,
                                    i_ppara, i_pperp,
                                    0, 0,
                                    n_phi, n_z,
                                    n_ppara, n_pperp,
                                    1, 1);

                            if (!ps_tmp.err) {
                                for (int i_mc=0; i_mc<n_per_vertex; ++i_mc) {
                                    ps_tmp.id = i;

                                    if (use_hermite) {
                                        ps_tmp.hermite_weights = hermiteW[i_mc] / PI2E0_5;
                                        ps_tmp.hermite_knots = hermiteK[i_mc];
                                        ps_tmp.use_hermite = 1;
                                    } else {
                                        ps_tmp.hermite_weights = 1. / n_per_vertex;
                                        ps_tmp.use_hermite = 0;
                                    }
                                    memcpy(*ps + i, &ps_tmp, sizeof(particle_state));
                                    (*ps_indexes)[i] = index;
                                    i++;
                                }
                            }
                        }
                    }
                } else {
                    for (int i_pr = 0; i_pr < n_pr; ++i_pr) {
                        pr = (max_pr - min_pr) * i_pr / n_pr + min_pr;
                        for (int i_pphi = 0; i_pphi < n_pphi; ++i_pphi) {
                            pphi = (max_pphi - min_pphi) * i_pphi / n_pphi + min_pphi;
                            for (int i_pz = 0; i_pz < n_pz; ++i_pz) {
                                pz = (max_pz - min_pz) * i_pz / n_pz + min_pz;
                                bmc_init_fo_particle(&p_tmp, r, phi, z, pr, pphi, pz, T1, i+1);
                                particle_input_to_state(&p_tmp, &ps_tmp, Bdata);

                                unsigned long index = dist_6D_index(i_r, i_phi, i_z,
                                    i_pr, i_pphi, i_pz,
                                    0, 0,
                                    n_phi, n_z,
                                    n_pr, n_pphi, n_pz,
                                    1, 1);
                                if (!ps_tmp.err) {
                                    for (int i_mc=0; i_mc<n_per_vertex; ++i_mc) {
                                        ps_tmp.id = i;
                                        memcpy(*ps + i, &ps_tmp, sizeof(particle_state));
                                        i++;
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    *n = i;
    printf("Initialized %d particles\n", i);

    return 0;
}

void buildImportantSamplingHistogram(
    int dist_length,
    real *histogram,
    dist_5D_offload_data* dist5D,
    plasma_data* plasma_data,
    B_field_data* Bdata,
    int importanceSamplingProbability,
    int importanceSamplingdensity
) {

    real *histogram_probability;
    FILE *p_probability;
    if (importanceSamplingProbability) {
        p_probability = fopen("distr_prob", "rb");
        histogram_probability = (real*)malloc(dist_length * sizeof(real));
        fread(histogram_probability, sizeof(real), dist_length, p_probability);
    }

    real r,phi,z,ppara,pperp, dens[10];
    for (int i=0; i<dist_length; i++) {

        histogram[i] = 1.;
        int i_x[5];
        compute_5d_coordinates_from_hist_index(i, i_x, &r, &phi, &z, &ppara, &pperp, dist5D);

        // avoid to fill velocity boundaries
        if ((i_x[4] == 0) || (i_x[4] == dist5D->n_pperp - 1) || (i_x[2] == 0) || (i_x[2] == dist5D->n_ppara - 1)) {
            histogram[i] = 0;
            continue;
        }

        if (importanceSamplingdensity) {

            real psi[1], rho[1], B_dB[15];
            B_field_eval_B_dB(B_dB, r, phi, z, T1, Bdata);
            B_field_eval_psi(psi, r, phi, z, T1, Bdata);
            B_field_eval_rho(rho, psi[0], Bdata);

            plasma_eval_dens(dens, rho[0], r, phi, z, 0, 1, plasma_data);

            // add space distribution from background plasma
            histogram[i] *= dens[0];

            // add velocity distribution
            if ((fabs(ppara) <= 6E-20) && (fabs(ppara) >= 0E-20) && (fabs(pperp) <= 6E-20) && (fabs(pperp) >= 0E-30)) {
                    histogram[i] *= 1;
            } else {
                    histogram[i] = 0;
            }
        } 
        if (importanceSamplingProbability) {
            histogram[i] *= histogram_probability[i];
        } 
    }
}

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
        int importanceSamplingdensity
    ) {
    // vacate the phase space to find the phase-space points in the mesh
    // suitable for simulating the BMC scheme

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);
    real* ptr = offload_unpack(offload_data, offload_array,
            sim_offload->B_offload_data.offload_array_length);
    B_field_init(&sim.B_data, &sim_offload->B_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->E_offload_data.offload_array_length);
    E_field_init(&sim.E_data, &sim_offload->E_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->plasma_offload_data.offload_array_length);
    plasma_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->neutral_offload_data.offload_array_length);
    neutral_init(&sim.neutral_data, &sim_offload->neutral_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->wall_offload_data.offload_array_length);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);


    real* distr0_array;
    diag_init_offload(&sim_offload->diag_offload_data, &distr0_array, 1);

    real r;
    real phi;
    real z;
    real ppara;
    real pperp;
    real pr;
    real pphi;
    real pz;

    int n_r, n_phi, n_z, n_ppara, n_pperp, n_pr, n_pphi, n_pz;
    real max_r, max_phi, max_z, max_ppara, max_pperp, max_pr, max_pphi, max_pz;
    real min_r, min_phi, min_z, min_ppara, min_pperp, min_pr, min_pphi, min_pz;
    dist_5D_offload_data dist5D = sim_offload->diag_offload_data.dist5D;
    n_r = dist5D.n_r;
    n_phi = dist5D.n_phi;
    n_z = dist5D.n_z;
    n_ppara = dist5D.n_ppara;
    n_pperp = dist5D.n_pperp;
    max_r = dist5D.max_r;
    max_phi = dist5D.max_phi;
    max_z = dist5D.max_z;
    max_ppara = dist5D.max_ppara;
    max_pperp = dist5D.max_pperp;
    min_r = dist5D.min_r;
    min_phi = dist5D.min_phi;
    min_z = dist5D.min_z;
    min_ppara = dist5D.min_ppara;
    min_pperp = dist5D.min_pperp;

    int dist_length = sim_offload->diag_offload_data.offload_array_length;
    real *histogram = (real*)malloc(dist_length * sizeof(real));

    buildImportantSamplingHistogram(dist_length, histogram, &dist5D, &sim.plasma_data, Bdata, importanceSamplingProbability, importanceSamplingdensity);

    int *nparticlesHistogram = (int*)malloc(dist_length * sizeof(int));

    real sum = 0;
    for (int i=0; i<dist_length; i++) {
       sum += histogram[i];
    }
    printf("Init initial sum %f\n", sum);
    *n = 0;
    for (int i=0; i<dist_length; i++) {
        nparticlesHistogram[i] = ceil(histogram[i] / sum * n_total);
        // if (histogram[i] / sum * n_total > 0) {
        //     printf("hist %d %f %d\n", i, histogram[i] / sum * n_total, nparticlesHistogram[i]);
        // }
        *n += nparticlesHistogram[i];
    }

    *ps = (particle_state *)malloc(*n * sizeof(particle_state));
    *ps_indexes = (int *)malloc(*n * sizeof(int));
    input_particle p_tmp; // tmp particle
    particle_state ps_tmp; // tmp particle

    int i = 0;
    for (int i_r = 0; i_r < n_r; ++i_r) {
        r = (max_r - min_r) * i_r / n_r + min_r;
        for (int i_phi = 0; i_phi < n_phi; ++i_phi) {
            phi = (max_phi - min_phi) * i_phi / n_phi + min_phi;
            for (int i_z = 0; i_z < n_z; ++i_z) {
                z = (max_z - min_z) * i_z / n_z + min_z;

                if (!wall_2d_inside(r, z, &sim.wall_data.w2d)) {
                    continue;
                }

                for (int i_ppara = 0; i_ppara < n_ppara; ++i_ppara) {
                    ppara = (max_ppara - min_ppara) * i_ppara / n_ppara + min_ppara;
                    for (int i_pperp = 0; i_pperp < n_pperp; ++i_pperp) {
                        pperp = (max_pperp - min_pperp) * i_pperp / n_pperp + min_pperp;
                        bmc_5D_to_particle_state(Bdata, r, phi, z, ppara, pperp, T1, 0, &ps_tmp);

                        unsigned long index = dist_5D_index(i_r, i_phi, i_z,
                                i_ppara, i_pperp,
                                0, 0,
                                n_phi, n_z,
                                n_ppara, n_pperp,
                                1, 1);

                        if (!ps_tmp.err) {
                            for (int i_mc=0; i_mc<nparticlesHistogram[index]; ++i_mc) {
                                ps_tmp.id = i;

                                ps_tmp.use_hermite = 0;
                                ps_tmp.hermite_weights = 1./nparticlesHistogram[index];

                                memcpy(*ps + i, &ps_tmp, sizeof(particle_state));
                                (*ps_indexes)[i] = index;
                                i++;
                            }
                        }
                    }
                }
            }
        }
    }



    printf("Initialized %d %d particles\n", *n, i);

    return 0;
}
int fmc_init_importance_sampling(
        int *n,
        particle_state** ps,
        int** ps_indexes,
        int n_total,
        int use_hermite,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array
    ) {
    // vacate the phase space to find the phase-space points in the mesh
    // suitable for simulating the BMC scheme

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);
    
    real* ptr = offload_array + 
            sim_offload->B_offload_data.offload_array_length +
            sim_offload->E_offload_data.offload_array_length +
            sim_offload->plasma_offload_data.offload_array_length +
            sim_offload->neutral_offload_data.offload_array_length;
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    real* distr0_array;
    diag_init_offload(&sim_offload->diag_offload_data, &distr0_array, 1);

    real r;
    real phi;
    real z;
    real ppara;
    real pperp;
    real pr;
    real pphi;
    real pz;
    real r_new;
    real phi_new;
    real z_new;
    real ppara_new;
    real pperp_new;
    real pr_new;
    real pphi_new;
    real pz_new;

    int n_r, n_phi, n_z, n_ppara, n_pperp, n_pr, n_pphi, n_pz;
    real max_r, max_phi, max_z, max_ppara, max_pperp, max_pr, max_pphi, max_pz;
    real min_r, min_phi, min_z, min_ppara, min_pperp, min_pr, min_pphi, min_pz;
    dist_5D_offload_data dist5D = sim_offload->diag_offload_data.dist5D;
    n_r = dist5D.n_r;
    n_phi = dist5D.n_phi;
    n_z = dist5D.n_z;
    n_ppara = dist5D.n_ppara;
    n_pperp = dist5D.n_pperp;
    max_r = dist5D.max_r;
    max_phi = dist5D.max_phi;
    max_z = dist5D.max_z;
    max_ppara = dist5D.max_ppara;
    max_pperp = dist5D.max_pperp;
    min_r = dist5D.min_r;
    min_phi = dist5D.min_phi;
    min_z = dist5D.min_z;
    min_ppara = dist5D.min_ppara;
    min_pperp = dist5D.min_pperp;

    // printf("5D: n_r %d n_phi %d n_z %d n_vpara %d n_vperp %d n_per_vertex %d\n", n_r, n_phi, n_z, n_ppara, n_pperp, n_per_vertex);
    *n = n_total;
    // print_out(VERBOSE_NORMAL, "Mesh size %d.\n", *n / n_per_vertex);

    *ps = (particle_state *)malloc(*n * sizeof(particle_state));
    *ps_indexes = (int *)malloc(*n * sizeof(int));
    input_particle p_tmp; // tmp particle
    particle_state ps_tmp; // tmp particle

    // open precomputed distr file
    FILE *fp;
    fp = fopen("distr_out", "rb");
    int dist_length = sim_offload->diag_offload_data.offload_array_length;
    printf("%d\n", dist_length);
    real *histogram = (real*)malloc(dist_length * sizeof(real));
    fread(histogram, sizeof(real), dist_length, fp);

    int weights[dist_length];
    for (int i=0; i< dist_length; i++) {
        weights[i] = 0;
    }

    unsigned long index_old, index_new;
    int i_x[5], i_x_new[5];
    i_x[0] = n_r / 2;
    i_x[1] = n_phi / 2;
    i_x[2] = n_z / 2;
    i_x[3] = n_ppara / 2;
    i_x[4] = n_pperp / 2;
    // int i_r = n_r / 2, i_phi = n_phi / 2, i_z = n_z / 2;
    // int i_ppara = n_ppara / 2, i_pperp = n_pperp / 2;
    // int i_r_new, i_phi_new, i_z_new, i_ppara_new, i_pperp_new;
    compute_element_5d_coordinates(i_x, &r, &phi, &z, &ppara, &pperp, &dist5D);

    index_old = dist_5D_index(i_x[0], i_x[1], i_x[2],
        i_x[3], i_x[4],
        0, 0,
        n_phi, n_z,
        n_ppara, n_pperp,
        1, 1);

    int accepted = 0;
    int i = 0;
    while (i < *n) { 

        importance_sampling_random_move(i_x, i_x_new, i);
        // printf("%f %f %f %f %f %d %d\n", r, z, dist5D.min_r, dist5D.min_z, dist5D.max_z, i_x_new[2], i_x[2]);
        compute_element_5d_coordinates(i_x_new, &r, &phi, &z, &ppara, &pperp, &dist5D);

        if (!wall_2d_inside(r, z, &sim.wall_data.w2d)) {
            continue;
        }

        printf("%d %d\n", i, *n);

        // while (!wall_2d_inside(r, z, &sim.wall_data.w2d)) {
        //     importance_sampling_random_move(&i_r, &i_phi, &i_z, &i_ppara, &i_pperp);
        //     compute_element_5d_coordinates(i_r, i_phi, i_z, i_ppara, i_pperp, &r_new, &phi_new, &z_new, &ppara_new, &pperp_new);
        // }

        index_new = dist_5D_index(i_x_new[0], i_x_new[1], i_x_new[2],
            i_x_new[3], i_x_new[4],
            0, 0,
            n_phi, n_z,
            n_ppara, n_pperp,
            1, 1);

        double r_rand = r2(); 
        // printf("%d %d %e %e\n", index_new, index_old, histogram[index_new] / histogram[index_old], histogram[index_old]);
        if ((histogram[index_new] / histogram[index_old]) > r_rand ) {
            // printf("pass\n");
            index_old = index_new;
            r = r_new;
            phi = phi_new;
            z = z_new;
            ppara = ppara_new;
            pperp = pperp_new;
            memcpy(i_x, i_x_new, 5*sizeof(int));
        }

        bmc_5D_to_particle_state(Bdata, r, phi, z, ppara, pperp, T1, i, &ps_tmp);

        if (!ps_tmp.err) {
            weights[index_old]++;
            ps_tmp.id = i;
            ps_tmp.use_hermite = 0;
            // ps_tmp.hermite_weights = 1. / n_per_vertex;
            memcpy(*ps + i, &ps_tmp, sizeof(particle_state));
            (*ps_indexes)[i] = index_old;
            i++;
        } else {
            printf("err\n");
        }
    }

    printf("Computing particles weights\n");

    for (int i = 0; i < *n; i++) {
        (*ps)[i].hermite_weights = 1. / weights[(*ps_indexes)[i]];
    }

    printf("Initialized %d particles\n", n);

    return 0;
}

void importance_sampling_random_move(int* i_x, int* i_x_new, int i) {
    int dim = i % 5;
    memcpy(i_x_new, i_x, 5 * sizeof(int));
    i_x_new[dim] += random_at_most(4) - 2;
}

void compute_5d_coordinates_from_hist_index(int i, int* i_x, real* r, real* phi, real* z, real* ppara, real* pperp, dist_5D_offload_data* dist5D) {
    int i1 = i;
    i_x[4] = i1 % dist5D->n_pperp;
    i1 = i1 / dist5D->n_pperp;
    i_x[3] = i1 % dist5D->n_ppara;
    i1 = i1 / dist5D->n_ppara;
    i_x[2] = i1 % dist5D->n_z;
    i1 = i1 / dist5D->n_z;
    i_x[1] = i1 % dist5D->n_phi;
    i_x[0] = i1 / dist5D->n_phi;
    compute_element_5d_coordinates(i_x, r, phi, z, ppara, pperp, dist5D);
}

void compute_element_5d_coordinates(int* i_x_new, real* r, real* phi, real* z, real* ppara, real* pperp, dist_5D_offload_data* dist) {
    *r = dist->min_r + (dist->max_r - dist->min_r) / dist->n_r * i_x_new[0];
    *phi = dist->min_phi + (dist->max_phi - dist->min_phi) / dist->n_phi * i_x_new[1];
    *z = dist->min_z + (dist->max_z - dist->min_z) / dist->n_z * i_x_new[2];
    *ppara = dist->min_ppara + (dist->max_ppara - dist->min_ppara) / dist->n_ppara * i_x_new[3];
    *pperp = dist->min_pperp + (dist->max_pperp - dist->min_pperp) / dist->n_pperp * i_x_new[4];
}

void bmc_5D_to_particle_state(
        B_field_data* Bdata,
        real r, real phi, real z,
        real ppara, real pperp,
        real t,
        int id,
        particle_state* ps
    ) {

    a5err err = 0;

    real psi[1], rho[1], B_dB[15];
    err = B_field_eval_B_dB(B_dB, r, phi, z, t, Bdata);
    if (!err) {
        err = B_field_eval_psi(psi, r, phi, z,
                                t, Bdata);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], Bdata);
    }

    if(!err) {
        ps->rho        = rho[0];
        ps->B_r        = B_dB[0];
        ps->B_phi      = B_dB[4];
        ps->B_z        = B_dB[8];
        ps->B_r_dr     = B_dB[1];
        ps->B_phi_dr   = B_dB[5];
        ps->B_z_dr     = B_dB[9];
        ps->B_r_dphi   = B_dB[2];
        ps->B_phi_dphi = B_dB[6];
        ps->B_z_dphi   = B_dB[10];
        ps->B_r_dz     = B_dB[3];
        ps->B_phi_dz   = B_dB[7];
        ps->B_z_dz     = B_dB[11];
    }

    /* Input is in (Ekin,xi) coordinates but state needs (mu,vpar) so we need to do that
        * transformation first. */
    real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
    real mu = 0.5 * pperp * pperp / MASS / Bnorm;
    if(!err && mu < 0)          {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        ps->n_t_subcycles = RK4_SUBCYCLES;
        ps->r        = r;
        ps->phi      = phi;
        ps->z        = z;
        ps->mu       = mu;
        ps->ppar     = ppara;
        ps->zeta     = 0;
        ps->mass     = MASS;
        ps->charge   = CHARGE;
        ps->anum     = 0;
        ps->znum     = 1;
        ps->weight   = 1;
        ps->time     = t;
        ps->theta    = atan2(ps->z-B_field_get_axis_z(Bdata, ps->phi),
                                ps->r-B_field_get_axis_r(Bdata, ps->phi));
        ps->id       = id;
        ps->endcond  = 0;
        ps->walltile = 0;
        ps->cputime  = 0;
    }

    /* Guiding center transformation to get particle coordinates */
    real rprt, phiprt, zprt, pr, pphi, pz;
    if(!err) {
        real pparprt, muprt, zetaprt;
        gctransform_guidingcenter2particle(
            ps->mass, ps->charge, B_dB,
            ps->r, ps->phi, ps->z, ps->ppar, ps->mu, ps->zeta,
            &rprt, &phiprt, &zprt, &pparprt, &muprt, &zetaprt);

        B_field_eval_B_dB(B_dB, rprt, phiprt, zprt, ps->time, Bdata);

        gctransform_pparmuzeta2prpphipz(
            ps->mass, ps->charge, B_dB,
            phiprt, pparprt, muprt, zetaprt,
            &pr, &pphi, &pz);
    }
    if(!err && rprt <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        ps->rprt   = rprt;
        ps->phiprt = phiprt;
        ps->zprt   = zprt;
        ps->p_r    = pr;
        ps->p_phi  = pphi;
        ps->p_z    = pz;

        ps->err = 0;
    } else {
        ps->err = err;
    }

}

void bmc_init_fo_particle(
        input_particle* p,
        real r, real phi, real z,
        real p_r, real p_phi, real p_z,
        real t,
        int id
    ) {
    p->p.r = r;
    p->p.phi = phi;
    p->p.z = z;
    p->p.p_r = p_r;
    p->p.p_phi = p_phi;
    p->p.p_z = p_z;
    p->p.mass = MASS;
    p->p.charge = CHARGE;
    p->p.anum = 0;
    p->p.znum = 1;
    p->p.weight = 1;
    p->p.time = t;
    p->p.id = id;
    p->type = input_particle_type_p;
}


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
        bool importance_sampling
    ) {

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nStarting Backward Monte Carlo. N particles: %d.\n", n_mpi_particles);

    /* Allow threads to spawn threads */
    omp_set_nested(1);

    // ps0 holds the initial particle states (constant space in vertices, and changing time)
    // ps1 are simulated and will hold the final state at eatch time step
    particle_state* ps0 = (particle_state*) malloc(n_mpi_particles * sizeof(particle_state));
    memcpy(ps0, ps1, n_mpi_particles * sizeof(particle_state));

    for(int i = 0; i < 50; i++) {
        printf("Particle %d %f %f %f %f %f %f %f\n", i, ps1[i].r, ps1[i].phi, ps1[i].z, ps1[i].ppar, ps1[i].rho, ps1[i].rprt, ps1[i].p_r);
    }

    // initialize distributions
    diag_data distr0, distr1;
    real* distr0_array, *distr1_array;
    diag_init_offload(&sim_offload->diag_offload_data, &distr0_array, 1);
    diag_init(&distr0, &sim_offload->diag_offload_data, distr0_array);
    diag_init_offload(&sim_offload->diag_offload_data, &distr1_array, 1);
    diag_init(&distr1, &sim_offload->diag_offload_data, distr1_array);
    int dist_length = sim_offload->diag_offload_data.offload_array_length;

    /* Initialize diagnostics offload data.
     * Separate arrays for host and target */
    // diagnostic might be useless for BMC. Remove this?
    real* diag_offload_array_mic0, *diag_offload_array_mic1, *diag_offload_array_host;
    #ifdef TARGET
        diag_init_offload(&sim_offload->diag_offload_data, &diag_offload_array_mic0, n_tot_particles);
        diag_init_offload(&sim_offload->diag_offload_data, &diag_offload_array_mic1, n_tot_particles);
    #else
        diag_init_offload(&sim_offload->diag_offload_data, &diag_offload_array_host, n_tot_particles);
    #endif


    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);
    real* ptr = offload_unpack(offload_data, offload_array,
            sim_offload->B_offload_data.offload_array_length);
    B_field_init(&sim.B_data, &sim_offload->B_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->E_offload_data.offload_array_length);
    E_field_init(&sim.E_data, &sim_offload->E_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->plasma_offload_data.offload_array_length);
    plasma_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->neutral_offload_data.offload_array_length);
    neutral_init(&sim.neutral_data, &sim_offload->neutral_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->wall_offload_data.offload_array_length);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    // setup time end conditions.
    // By setting the end time to be the initial time,
    // the simulation is forced to end after 1 timestep
    sim_offload->endcond_max_simtime = T1;

    // set time in particle states
    for(int i = 0; i < n_mpi_particles; i++) {
        ps1[i].time = T0;
        ps0[i].time = T0;
    }

    // simulate one step of all needed particles
    fmc_simulation(ps1, sim_offload, offload_data, offload_array,
                        mic1_start, mic1_end, mic0_start, mic0_end, host_start, host_end, n_mic, n_host,
                        diag_offload_array_host, diag_offload_array_mic0, diag_offload_array_mic1);
    

    // // // Update the probability distribution
    int n_updated;
    if (distr0.dist5D_collect) {
        n_updated = fmc_update_distr5D_from_states(&distr1.dist5D, &distr0.dist5D, ps1_indexes, ps1, ps0, n_mpi_particles, &(sim.wall_data.w2d));
        
    } else {
        // TODO: Full orbit
    }

    // // shift distributions. Required since distr1 is partitioned through all the MPI nodes,
    // // and can't be written directly to disk
    diag_move_distribution(sim_offload, &distr0, &distr1, &n_updated);
    printf("Updated %d\n", n_updated);

    real sum = 0, dens[5];
    for (int i=0; i < dist_length; i++) {
        sum += distr0.dist5D.histogram[i];
    }
    printf("Value of integrated signal:%f\n", sum);

    write_probability_distribution(sim_offload, &distr0, distr0_array, mpi_rank, true);

    // // Free diagnostic data
    #ifdef TARGET
        diag_free_offload(&sim_offload->diag_offload_data, &diag_offload_array_mic0);
        diag_free_offload(&sim_offload->diag_offload_data, &diag_offload_array_mic1);
    #else
        diag_free_offload(&sim_offload->diag_offload_data, &diag_offload_array_host);
    #endif

    return 0;

}

/**
 * Write the probability distribution to the output HDF file
 * 
 * @param sim_offload pointer to simulation offload data
 * @param distr pointer to distribution struct
 * @param distr_array pointer to offload array of the distribution
 * @param dist_length length of the distribution data
 * @param mpi_rank MPI node rank id
 * 
 **/
void write_probability_distribution(
    sim_offload_data* sim_offload,
    diag_data* distr,
    real* distr_array,
    int mpi_rank,
    bool write_for_importance_sampling
) {

    FILE *fp;

    if (write_for_importance_sampling)
        fp = fopen("distr_prob", "wb");

    int dist_length = sim_offload->diag_offload_data.offload_array_length;

    if (write_for_importance_sampling)
        fwrite(distr->dist5D.histogram, sizeof(distr->dist5D.histogram), dist_length, fp);

    for (int i = 0; i<dist_length; i++) {
        // printf("%f\n", distr->dist5D.histogram[i]);
        if (sim_offload->diag_offload_data.dist5D_collect) {
            if (distr->dist5D.histogram[i] > 1.0001) {
                // printf("Warning: unpysical probability: %f\n", distr->dist5D.histogram[i]);
            }
        }
        if (sim_offload->diag_offload_data.dist6D_collect) {
            if (distr->dist6D.histogram[i] > 1.0001) {
                printf("Warning: unpysical probability: %f\n", distr->dist6D.histogram[i]);
            }
        }
    }
    if (write_for_importance_sampling)
        fclose(fp);

    /* Combine distribution and write it to HDF5 file */
    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nWriting BMC probability distribution.\n");
    if (mpi_rank == 0) {
        int err_writediag = 0;
        err_writediag = hdf5_interface_write_diagnostics(sim_offload, distr_array, sim_offload->hdf5_out);
        if(err_writediag) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                    "\nWriting diagnostics failed.\n"
                    "See stderr for details.\n");
        }
        else {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                    "BMC distributions written.\n");
        }

    }
}
