#include "bmc.h"

void bmc_setup_endconds(sim_offload_data* sim) {
    // TODO: error assertion check
    sim->endcond_active = endcond_tmax | endcond_wall;
    sim->endcond_max_simtime = sim->fix_usrdef_val;

    // force n_time and n_q to be 1 in sim struct
    diag_offload_data diag_data = sim->diag_offload_data;
    if (diag_data.dist6D_collect) {
        diag_data.dist6D.n_time = 1;
        diag_data.dist6D.n_q = 1;
    }
    if (diag_data.dist5D_collect) {
        diag_data.dist5D.n_time = 1;
        diag_data.dist5D.n_q = 1;
    }
}

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
    ) {

    // hardcoded times
    // TODO: use input hdf
    // TODO: make sure distributions have nt = 1
    double ti = 0;
    double tf = 10;
    double dt = 0.1;

    /* Allow threads to spawn threads */
    omp_set_nested(1);

    // create 2 copies of the particles.
    // ps0 holds the initial particle states (constant space in vertices, and changing time)
    // ps1 are simulated and will hold the final state at eatch time step
    particle_state* ps1 = (particle_state*) malloc(n_mpi_particles * sizeof(particle_state));
    particle_state* ps0 = (particle_state*) malloc(n_mpi_particles * sizeof(particle_state));
    for(int i = 0; i < n_mpi_particles; i++) {
        particle_input_to_state(&p_mpi[i], &ps1[i], Bdata);
        particle_input_to_state(&p_mpi[i], &ps0[i], Bdata);
    }

    // initialize distributions
    diag_data diag0, diag1;
    real* diag_offload_array;
    diag_init_offload(&sim_offload->diag_offload_data, &diag_offload_array, n_tot_particles);
    diag_init(&diag0, &sim_offload->diag_offload_data, diag_offload_array);
    diag_init(&diag1, &sim_offload->diag_offload_data, diag_offload_array);

    int dist_length = sim_offload->diag_offload_data.offload_array_length;

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);

    for (double t=tf; t >= ti; t -= dt) {

        // set time in particle states
        for(int i = 0; i < n_mpi_particles; i++) {
            ps1[i].time = t;
            ps0[i].time = t;
        }

        // simulate one step of all needed particles
        bmc_simulate_particles(ps1, n_tot_particles, sim_offload, offload_data, offload_array, mic1_start, mic1_end, mic0_start, mic0_end, host_start, host_end, n_mic, n_host);

        // Update the probability distribution
        bmc_update_particles_diag(n_mpi_particles, ps0, ps1, diag0, diag1, sim, n_montecarlo_steps);

        // reduce and shift distributions
        #ifdef MPI
            // MPI_Allreduce(diag1.histogram, diag0.histogram, dist_length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            mpi_reduce_distribution(sim_offload, &diag0, &diag1, dist_length);
        #else
            // memcpy(diag0.histogram, diag1.histogram, dist_length*sizeof(real));
            diag_copy_distribution(&diag0, &diag1);
        #endif
    } 
}

void mpi_reduce_distribution(sim_offload_data* sim, diag_data* diag0, diag_data* diag1, int dist_length) {
    if (sim->diag_offload_data.dist5D_collect) {
        MPI_Allreduce(diag1->dist5D.histogram, diag0->dist5D.histogram, dist_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    if (sim->diag_offload_data.dist6D_collect) {
        MPI_Allreduce(diag1->dist6D.histogram, diag0->dist6D.histogram, dist_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}

void bmc_init_particles(
        int *n,
        input_particle **p,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        int n_montecarlo_steps
    ) {
    // vacate the phase space to find the phase-space points in the mesh
    // suitable for simulating the BMC scheme

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);

    // TODO: use input hdf for time
    double tf = 10;
    *n = 0;

    real r;
    real phi;
    real z;
    real vpara;
    real vperp;
    real vr;
    real vphi;
    real vz;
    real B_dB[15];

    int n_r, n_phi, n_z, n_vpara, n_vperp, n_vr, n_vphi, n_vz;
    int max_r, max_phi, max_z, max_vpara, max_vperp, max_vr, max_vphi, max_vz;
    int min_r, min_phi, min_z, min_vpara, min_vperp, min_vr, min_vphi, min_vz;
    if (sim_offload->diag_offload_data.dist5D_collect) {
        dist_5D_offload_data dist5D = sim_offload->diag_offload_data.dist5D;
        n_r = dist5D.n_r;
        n_phi = dist5D.n_phi;
        n_z = dist5D.n_z;
        n_vpara = dist5D.n_vpara;
        n_vperp = dist5D.n_vperp;
        max_r = dist5D.max_r;
        max_phi = dist5D.max_phi;
        max_z = dist5D.max_z;
        max_vpara = dist5D.max_vpara;
        max_vperp = dist5D.max_vperp;
        min_r = dist5D.min_r;
        min_phi = dist5D.min_phi;
        min_z = dist5D.min_z;
        min_vpara = dist5D.min_vpara;
        min_vperp = dist5D.min_vperp;
    } else {
        dist_6D_offload_data dist6D = sim_offload->diag_offload_data.dist6D;
        n_r = dist6D.n_r;
        n_phi = dist6D.n_phi;
        n_z = dist6D.n_z;
        n_vr = dist6D.n_vr;
        n_vphi = dist6D.n_vphi;
        n_vz = dist6D.n_vz;
        max_r = dist6D.max_r;
        max_phi = dist6D.max_phi;
        max_z = dist6D.max_z;
        max_vr = dist6D.max_vr;
        max_vphi = dist6D.max_vphi;
        max_vz = dist6D.max_vz;
        min_r = dist6D.min_r;
        min_phi = dist6D.min_phi;
        min_z = dist6D.min_z;
        min_z = dist6D.min_z;
        min_vr = dist6D.min_vr;
        min_vphi = dist6D.min_vphi;
    }


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
                    // compute magnetic field
                    B_field_eval_B_dB(B_dB, r, phi, z, tf, Bdata);

                    for (int i_vpara = 0; i_vpara < n_vpara; ++i_vpara) {
                        vpara = (max_vpara - min_vpara) * i_vpara / n_vpara + min_vpara;
                        for (int i_vperp = 0; i_vperp < n_vperp; ++i_vperp) {
                            vperp = (max_vperp - min_vperp) * i_vperp / n_vperp + min_vperp;
                            bmc_5D_to_fo(B_dB, r, phi, z, vpara, vperp, &vr, &vphi, &vz);
                        }
                    }
                } else {
                    for (int i_vr = 0; i_vr < n_vr; ++i_vr) {
                        vr = (max_vr - min_vr) * i_vr / n_vr + min_vr;
                        for (int i_vphi = 0; i_vphi < n_vphi; ++i_vphi) {
                            vphi = (max_vphi - min_vphi) * i_vphi / n_vphi + min_vphi;
                            for (int i_vz = 0; i_vz < n_vz; ++i_vz) {
                                vz = (max_vz - min_vz) * i_vz / n_vz + min_vz;
                            }
                        }
                    }
                }

                // create as many particles as needed for the integration of this mesh vertex
                for (int i=0; i<n_montecarlo_steps; ++i) {
                    *n = *n + 1;
                    *p = (input_particle *)realloc(*p, *n * sizeof(input_particle));
                    bmc_init_fo_particle(*p+*n, r, phi, z, vr, vphi, vz, tf, *n);
                }
            }
        }
    }
}

void bmc_5D_to_fo(
        real B_dB[15],
        real r, real phi, real z,
        real vpara, real vperp,
        real *vr,
        real *vphi,
        real *vz
    ) {

    // find the parallel unit vector
    real Bpar_norm = sqrt(B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]);
    real e_par_r = B_dB[0] / Bpar_norm;
    real e_par_phi = B_dB[4] / Bpar_norm;
    real e_par_z = B_dB[8] / Bpar_norm;
    
    // find a perpendicular unit vector --> i.e. normalized of (-Bz, 0, Br)
    real Bperp_norm = sqrt(B_dB[0]*B_dB[0] + B_dB[8]*B_dB[8]);
    real e_perp_r = - B_dB[8] / Bperp_norm;
    real e_perp_phi = 0;
    real e_perp_z = B_dB[0] / Bperp_norm;

    *vr = vpara * e_par_r + vperp * e_perp_r;
    *vphi = vpara * e_par_phi + vperp * e_perp_phi;
    *vz = vpara * e_par_z + vperp * e_perp_z;
}

void bmc_init_fo_particle(
        input_particle* p,
        real r, real phi, real z,
        real v_r, real v_phi, real v_z,
        real t,
        int id
    ) {
    p->p.r = r;
    p->p.phi = phi;
    p->p.z = z;
    p->p.v_r = v_r;
    p->p.v_phi = v_phi;
    p->p.v_z = v_z;
    p->p.mass = 1;
    p->p.charge = 1;
    p->p.anum = 0;
    p->p.znum = 1;
    p->p.weight = 1;
    p->p.time = t;
    p->p.id = id;
    p->type = input_particle_type_p;
}

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
        int n_host
    ) {

    /* Initialize diagnostics offload data.
     * Separate arrays for host and target */
    // diagnostic might be useless for BMC. Remove this?
    #ifdef TARGET
        real* diag_offload_array_mic0;
        real* diag_offload_array_mic1;
        diag_init_offload(&sim->diag_offload_data, &diag_offload_array_mic0, n_tot_particles);
        diag_init_offload(&sim->diag_offload_data, &diag_offload_array_mic1, n_tot_particles);
    #else
        real* diag_offload_array_host;
        diag_init_offload(&sim->diag_offload_data, &diag_offload_array_host, n_tot_particles);
    #endif

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
                simulate(1, n_mic, ps, &sim, &offload_data, offload_array,
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
                simulate(2, n_mic, ps+n_mic, &sim, &offload_data, offload_array,
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
                simulate(0, n_host, ps+2*n_mic, &sim, &offload_data,
                    offload_array, diag_offload_array_host);
                *host_end = omp_get_wtime();
            }
        #endif
    }
}