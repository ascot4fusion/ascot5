#include "bmc.h"

#define TIMESTEP 1E-7 // TODO use input HDF
#define T0 0.
#define T1 4E-7
// #define T1 
#define MASS 9.10938356E-31
#define CHARGE 1.60217662E-19
#define RK4_SUBCYCLES 1
#define PI2E0_5 2.50662827463

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
        printf("Particle %d %f %f %f %f %f %f %f\n", i, ps[i].r, ps[i].phi, ps[i].z, ps[i].vpar, ps[i].rho, ps[i].rprt, ps[i].rdot);
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
        backward_monte_carlo_gc(ps, ps_indexes, n_mpi_particles, sim_offload, &sim, offload_data, offload_array,
                                &distr0, &distr1, Bdata);
    } else {
        // TODO: FULL ORBIT
    }

    write_probability_distribution(sim_offload, &distr0, distr0_array, mpi_rank);

    // Free sitribution data
    diag_free_offload(&sim_offload->diag_offload_data, &distr1_array);
    diag_free_offload(&sim_offload->diag_offload_data, &distr0_array);

    return 0;

}

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
    ) {

    // int n_simd_particles = n_mpi_particles / NSIMD;
    particle_simd_gc *p0, *p1;
    int n_simd_particles = init_simd_gc_particles(ps, n_mpi_particles, &p0, Bdata);
    init_simd_gc_particles(ps, n_mpi_particles, &p1, Bdata);
    
    for (double t=T1; t >= T0; t -= TIMESTEP) {

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
            int n_host = 0;
        #else
            int n_mic = 0;
            int n_host = n_simd_particles;
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
                        offload_array[0:offload_data.offload_array_length] \
                    )
                    bmc_simulate_timestep_gc(n_mic, p1 + n_mic, sim_offload, offload_data, offload_array, TIMESTEP, RK4_SUBCYCLES);
                }
            #endif
            #ifdef TARGET >= 2
                #pragma omp section
                {
                    #pragma omp target device(1) map( \
                        p1[n_mic:2*n_mic], \
                        offload_array[0:offload_data.offload_array_length] \
                    )
                    bmc_simulate_timestep_gc(n_mic, p1 + n_mic, sim_offload, offload_data, offload_array, TIMESTEP, RK4_SUBCYCLES);
                }
            #endif
            #ifndef TARGET
                #pragma omp section
                {
                    bmc_simulate_timestep_gc(n_simd_particles, p1, sim_offload, offload_data, offload_array, TIMESTEP, RK4_SUBCYCLES);
                }
            #endif
        }

        // // Update the probability distribution
        // bmc_update_distr_gc(n_mpi_particles, ps0, ps1, ps1_indexes, &distr0, &distr1, &sim);
        int n_updated = bmc_update_distr5D(&distr1->dist5D, &distr0->dist5D, p0_indexes, p1, p0, n_simd_particles, &(sim->wall_data.w2d));
        
        printf("Time %f Updated %d\n", t, n_updated);

        // // shift distributions
        diag_move_distribution(sim_offload, distr0, distr1);

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
    real vpara;
    real vperp;
    real vr;
    real vphi;
    real vz;

    int n_r, n_phi, n_z, n_vpara, n_vperp, n_vr, n_vphi, n_vz;
    real max_r, max_phi, max_z, max_vpara, max_vperp, max_vr, max_vphi, max_vz;
    real min_r, min_phi, min_z, min_vpara, min_vperp, min_vr, min_vphi, min_vz;
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

    if (sim_offload->diag_offload_data.dist5D_collect) {
        printf("5D: n_r %d n_phi %d n_z %d n_vpara %d n_vperp %d n_per_vertex %d\n", n_r, n_phi, n_z, n_vpara, n_vperp, n_per_vertex);
        *n = n_r * n_phi * n_z * n_vpara * n_vperp * n_per_vertex;
    } else {
        printf("6D: n_r %d n_phi %d n_z %d n_vr %d n_vphi %d n_vz %d n_per_vertex %d\n", n_r, n_phi, n_z, n_vr, n_vphi, n_vz, n_per_vertex);
        *n = n_r * n_phi * n_z * n_vr * n_vphi * n_vz * n_per_vertex;
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
                    for (int i_vpara = 0; i_vpara < n_vpara; ++i_vpara) {
                        vpara = (max_vpara - min_vpara) * i_vpara / n_vpara + min_vpara;
                        for (int i_vperp = 0; i_vperp < n_vperp; ++i_vperp) {
                            vperp = (max_vperp - min_vperp) * i_vperp / n_vperp + min_vperp;
                            bmc_5D_to_particle_state(Bdata, r, phi, z, vpara, vperp, T1, i, &ps_tmp);

                            unsigned long index = dist_5D_index(i_r, i_phi, i_z,
                                    i_vpara, i_vperp,
                                    0, 0,
                                    n_phi, n_z,
                                    n_vpara, n_vperp,
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
                    for (int i_vr = 0; i_vr < n_vr; ++i_vr) {
                        vr = (max_vr - min_vr) * i_vr / n_vr + min_vr;
                        for (int i_vphi = 0; i_vphi < n_vphi; ++i_vphi) {
                            vphi = (max_vphi - min_vphi) * i_vphi / n_vphi + min_vphi;
                            for (int i_vz = 0; i_vz < n_vz; ++i_vz) {
                                vz = (max_vz - min_vz) * i_vz / n_vz + min_vz;
                                bmc_init_fo_particle(&p_tmp, r, phi, z, vr, vphi, vz, T1, i+1);
                                particle_input_to_state(&p_tmp, &ps_tmp, Bdata);

                                unsigned long index = dist_6D_index(i_r, i_phi, i_z,
                                    i_vr, i_vphi, i_vz,
                                    0, 0,
                                    n_phi, n_z,
                                    n_vr, n_vphi, n_vz,
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

void bmc_5D_to_particle_state(
        B_field_data* Bdata,
        real r, real phi, real z,
        real vpara, real vperp,
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
    real mu = 0.5 * vperp * vperp * MASS / Bnorm;
    if(!err && mu < 0)          {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
    if(!err && vpara >= CONST_C) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        ps->n_t_subcycles = RK4_SUBCYCLES;
        ps->r        = r;
        ps->phi      = phi;
        ps->z        = z;
        ps->mu       = mu;
        ps->vpar     = vpara;
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
    real rprt, phiprt, zprt, vR, vphi, vz;
    if(!err) {
        real vparprt, muprt, zetaprt;
        gctransform_guidingcenter2particle(
            ps->mass, ps->charge, B_dB,
            ps->r, ps->phi, ps->z, ps->vpar, ps->mu, ps->zeta,
            &rprt, &phiprt, &zprt, &vparprt, &muprt, &zetaprt);

        B_field_eval_B_dB(B_dB, rprt, phiprt, zprt, ps->time, Bdata);
        if(!err && vparprt >= CONST_C) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && -vparprt >= CONST_C) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

        gctransform_vparmuzeta2vRvphivz(
            ps->mass, ps->charge, B_dB,
            phiprt, vparprt, muprt, zetaprt,
            &vR, &vphi, &vz);
    }
    if(!err && rprt <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        ps->rprt   = rprt;
        ps->phiprt = phiprt;
        ps->zprt   = zprt;
        ps->rdot   = vR;
        ps->phidot = vphi / ps->rprt;
        ps->zdot   = vz;

        ps->err = 0;
    } else {
        ps->err = err;
    }

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
        int mpi_rank
    ) {

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nStarting Backward Monte Carlo. N particles: %d.\n", n_mpi_particles);

    /* Allow threads to spawn threads */
    omp_set_nested(1);

    // ps0 holds the initial particle states (constant space in vertices, and changing time)
    // ps1 are simulated and will hold the final state at eatch time step
    particle_state* ps0 = (particle_state*) malloc(n_mpi_particles * sizeof(particle_state));
    memcpy(ps0, ps1, n_mpi_particles * sizeof(particle_state));

    for(int i = 0; i < 50; i++) {
        printf("Particle %d %f %f %f %f %f %f %f\n", i, ps1[i].r, ps1[i].phi, ps1[i].z, ps1[i].vpar, ps1[i].rho, ps1[i].rprt, ps1[i].rdot);
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
    real* ptr = offload_array + 
        sim_offload->B_offload_data.offload_array_length +
        sim_offload->E_offload_data.offload_array_length +
        sim_offload->plasma_offload_data.offload_array_length +
        sim_offload->neutral_offload_data.offload_array_length;
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
    fmc_simulation(ps1, n_tot_particles, sim_offload, offload_data, offload_array,
                        mic1_start, mic1_end, mic0_start, mic0_end, host_start, host_end, n_mic, n_host,
                        diag_offload_array_host, diag_offload_array_mic0, diag_offload_array_mic1);

    // // Update the probability distribution
    if (distr0.dist5D_collect) {
        // int n_updated = bmc_update_distr5D_from_states(n_mpi_particles, ps0, ps1, ps1_indexes, &distr0.dist5D, &distr1.dist5D, &(sim.wall_data.w2d));
        int n_updated = bmc_update_distr5D_from_states(&distr1.dist5D, &distr0.dist5D, ps1_indexes, ps1, ps0, n_mpi_particles, &(sim.wall_data.w2d));
        
        printf("Updated %d\n", n_updated);
    } else {
        // TODO: Full orbit
    }

    // shift distributions. Required since distr1 is partitioned through all the MPI nodes,
    // and can't be written directly to disk
    diag_move_distribution(sim_offload, &distr0, &distr1);

    write_probability_distribution(sim_offload, &distr0, distr0_array, mpi_rank);

    // Free diagnostic data
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
    int mpi_rank
) {
        int dist_length = sim_offload->diag_offload_data.offload_array_length;
        for (int i = 0; i<dist_length; i++) {
        if (sim_offload->diag_offload_data.dist5D_collect) {
            if (distr->dist5D.histogram[i] > 1.0001) {
                printf("Warning: unpysical probability: %f\n", distr->dist5D.histogram[i]);
            }
        }
        if (sim_offload->diag_offload_data.dist6D_collect) {
            if (distr->dist6D.histogram[i] > 1.0001) {
                printf("Warning: unpysical probability: %f\n", distr->dist6D.histogram[i]);
            }
        }
    }

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