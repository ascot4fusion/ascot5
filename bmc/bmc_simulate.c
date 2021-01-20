#include "bmc_simulate.h"

/**
 * @brief Execute marker simulation
 *
 * This simulates markers using given inputs and options. All different types of
 * simulations are initialized and run via this function.
 *
 * @param id target id where this function is executed, zero if on host
 * @param n_particles total number of markers to be simulated
 * @param p pointer to array storing all marker states to be simulated
 * @param sim_offload pointer to simulation offload data
 * @param offload_data pointer to the rest of the offload data
 * @param offload_array pointer to input data offload array
 * @param diag_offload_array pointer to diagnostics offload array
 *
 * @todo Reorganize this function so that it conforms to documentation.
 */

#define PI2E0_5 2.50662827463

void bmc_simulate_timestep_gc(int n_simd_particles, int n_coll_simd_particles, particle_simd_gc* p, particle_simd_gc* p_coll,
        int n_hermite_knots,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        real h, int n_rk4_subcycles
    ) {

    /**************************************************************************/
    /* 1. Input offload data is unpacked and initialized by calling           */
    /*    respective init functions.                                          */
    /*                                                                        */
    /**************************************************************************/
    sim_data sim;
    sim_init(&sim, sim_offload);

    real* ptr;
    offload_data->unpack_pos = 0;
    ptr = offload_unpack(offload_data, offload_array,
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

    random_init(&sim.random_data, 0);

    /******************************************************************/
    /* 2. Prepare time steps */
    /*                                                                */
    /******************************************************************/
    real h_rk4[NSIMD] __memalign__;
    real h_coll[NSIMD] __memalign__;
    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        h_rk4[i] = h / n_rk4_subcycles;
        h_coll[i] = h;
    }

    /* MAIN SIMULATION LOOP
     * - Store current state
     * - Integrate motion due to background EM-field (orbit-following)
     * - Integrate scattering due to Coulomb collisions
     * - Advance time
     * - Check for end condition(s)
     * - Update diagnostics
     */
    particle_simd_gc p0;
    for (int i_simd = 0; i_simd < n_simd_particles; i_simd++) {

        /* RK4 method for orbit-following */
        for (int nt = 0; nt < n_rk4_subcycles; ++nt) {
            memcpy(&p0, p + i_simd, sizeof(particle_simd_gc));
            step_gc_rk4(p + i_simd, h_rk4, &sim.B_data, &sim.E_data);

            #pragma omp simd
            for(int i = 0; i < NSIMD; i++) {
                if (p[i_simd].running[i]) {
                    real w_coll = 0;
                    int tile = wall_hit_wall(p0.r[i], p0.phi[i], p0.z[i],
                                            p[i_simd].r[i], p[i_simd].phi[i], p[i_simd].z[i],
                                            &(sim.wall_data), &w_coll);
                    if(tile > 0) {
                        p[i_simd].walltile[i] = tile;
                        p[i_simd].endcond[i] |= endcond_wall;
                        p[i_simd].running[i] = 0;
                    }
                }
            }
        }
    }

    // copy the result of RK4 to N_HERMITE different copy of the particle.
    // compute the hermite weights for each copy of the particle
    copy_particles_simd_to_coll_simd(n_simd_particles, n_hermite_knots, p, p_coll);
    /* Euler-Maruyama method for collisions */
    if(sim.enable_clmbcol) {
        for (int i_simd = 0; i_simd < n_coll_simd_particles; i_simd++) {
            mccc_gc_euler(p_coll + i_simd, h_coll, &sim.B_data, &sim.plasma_data,
                            &sim.random_data, &sim.mccc_data);
        }
    }
}

void init_particles_coll_simd_hermite(int n_simd_particles, int n_hermite_knots,
        particle_simd_gc* p_coll
    ) {
    real hermiteK[5] = {-2.856970, -1.355626, 0.000000, 1.355626, 2.856970};
    real hermiteW[5] = {0.028218, 0.556662, 1.336868, 0.556662, 0.028218};
    int i_coll = 0;
    for (int i=0; i < n_simd_particles; i++) {
        for (int j=0; j < NSIMD; j++) {
            for (int k = 0; k < n_hermite_knots; k++) {
                p_coll[i_coll / NSIMD].hermite_knots[i_coll % NSIMD] = hermiteK[k];
                p_coll[i_coll / NSIMD].hermite_weights[i_coll % NSIMD] = hermiteW[k] / PI2E0_5;
                p_coll[i_coll / NSIMD].use_hermite[i_coll % NSIMD] = 1;
                i_coll++;
            }
        }
    }
}

void copy_particles_simd_to_coll_simd(int n_simd_particles, int n_hermite_knots,
        particle_simd_gc* p, particle_simd_gc* p_coll
    ) {
    int i_coll = 0;
    for (int i=0; i < n_simd_particles; i++) {
        for (int j=0; j < NSIMD; j++) {
            for (int k = 0; k < n_hermite_knots; k++) {
                particle_copy_gc(&p[i], j, &p_coll[i_coll / NSIMD], i_coll % NSIMD);

                i_coll++;
            }
        }
    }
}

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
    ) {

    offload_data->unpack_pos = 0;

    /* Actual marker simulation happens here. Threads are spawned which
    * distribute the execution between target(s) and host. Both input and
    * diagnostic offload arrays are mapped to target. Simulation is initialized
    * at the target and completed within the simulate() function.*/
    #pragma omp parallel sections num_threads(3)
    {
        /* Run simulation on first target */
        #if TARGET >= 1
            #pragma omp section
            {
                *mic0_start = omp_get_wtime();

                #pragma omp target device(0) map( \
                    ps[0:n_mic], \
                    offload_array[0:offload_data.offload_array_length], \
                    diag_offload_array_mic0[0:sim.diag_offload_data.offload_array_length] \
                )
                simulate(1, n_mic, ps, sim, offload_data, offload_array,
                    diag_offload_array_mic0);

                *mic0_end = omp_get_wtime();
            }
        #endif

        /* Run simulation on second target */
        #ifdef TARGET >= 2
            #pragma omp section
            {
                *mic1_start = omp_get_wtime();

                #pragma omp target device(1) map( \
                    ps[n_mic:2*n_mic], \
                    offload_array[0:offload_data.offload_array_length], \
                    diag_offload_array_mic1[0:sim.diag_offload_data.offload_array_length] \
                )
                simulate(2, n_mic, ps+n_mic, sim, offload_data, offload_array,
                    diag_offload_array_mic1);

                *mic1_end = omp_get_wtime();
            }
        #endif
            /* No target, marker simulation happens where the code execution began.
            * Offloading is only emulated. */
        #ifndef TARGET
            #pragma omp section
            {
                *host_start = omp_get_wtime();
                simulate(0, n_host, ps+2*n_mic, sim, offload_data,
                    offload_array, diag_offload_array_host);
                *host_end = omp_get_wtime();
            }
        #endif
    }
}