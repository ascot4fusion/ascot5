#include "bmc.h"

void bmc_setup_endconds(sim_offload_data* sim, real timestep) {
    sim->endcond_active = endcond_tmax | endcond_wall;
    sim->fix_usrdef_val = timestep;
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
        int mpi_rank,
        real t1,
        real t0,
        real h,
        int rk4_subcycles
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
                                &distr0, &distr1, Bdata, t1, t0, h, rk4_subcycles);
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
        B_field_data* Bdata,
        real t1,
        real t0,
        real h,
        int rk4_subcycles
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

    for (double t=t1; t >= t0; t -= h) {

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
                    bmc_simulate_timestep_gc(n_simd_particles, n_coll_simd, p1, pcoll1, n_hermite_knots, sim_offload, offload_data, offload_array, h, rk4_subcycles);
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
    sim_offload->endcond_max_simtime = t1;

    // set time in particle states
    for(int i = 0; i < n_mpi_particles; i++) {
        ps1[i].time = t0;
        ps0[i].time = t0;
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
