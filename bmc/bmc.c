#include "bmc.h"

#define TIMESTEP 1E-7 // TODO use input HDF
#define T0 0.999999
#define T1 1.
#define N_MC_STEPS 1
#define MASS 9.10938356E-31
#define CHARGE 1.60217662E-19

void bmc_setup_endconds(sim_offload_data* sim) {
    sim->endcond_active = endcond_tmax | endcond_wall;
    sim->fix_usrdef_val = TIMESTEP;
    sim->fix_usrdef_use = 1;

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

int backward_monte_carlo(
        int n_tot_particles,
        int n_mpi_particles,
        particle_state* ps1,
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

    // for(int i = 0; i < 50; i++) {
    //     printf("Particle %d %f %f %f %f %f %f %f\n", i, ps1[i].r, ps1[i].phi, ps1[i].z, ps1[i].vpar, ps1[i].rho, ps1[i].rprt, ps1[i].rdot);
    // }

    // initialize distributions
    diag_data diag0, diag1;
    real* diag0_array, *diag1_array;
    diag_init_offload(&sim_offload->diag_offload_data, &diag0_array, n_tot_particles);
    diag_init_offload(&sim_offload->diag_offload_data, &diag1_array, n_tot_particles);
    diag_init(&diag0, &sim_offload->diag_offload_data, diag0_array);
    diag_init(&diag1, &sim_offload->diag_offload_data, diag1_array);

    int dist_length = sim_offload->diag_offload_data.offload_array_length;

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);

    for (double t=T1; t >= T0; t -= TIMESTEP) {

        // setup time end conditions.
        // By setting the end time to be the initial time,
        // the simulation is forced to end after 1 timestep
        sim_offload->endcond_max_simtime = t;

        // set time in particle states
        for(int i = 0; i < n_mpi_particles; i++) {
            ps1[i].time = t;
            ps0[i].time = t;
        }

        // simulate one step of all needed particles
        bmc_simulate_particles(ps1, n_tot_particles, sim_offload, offload_data, offload_array, mic1_start, mic1_end, mic0_start, mic0_end, host_start, host_end, n_mic, n_host);

        // Update the probability distribution
        bmc_update_particles_diag(n_mpi_particles, ps0, ps1, &diag0, &diag1, &sim, N_MC_STEPS);

        // shift distributions
        diag_copy_distribution(sim_offload, &diag0, &diag1, dist_length);

        // reset particle initial states to ps1
        memcpy(ps1, ps0, n_mpi_particles * sizeof(particle_state));
    } 

    /* Combine distribution and write it to HDF5 file */
    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nWriting BMC probability distribution.\n");
    if (mpi_rank == 0) {
        int err_writediag = 0;
        err_writediag = hdf5_interface_write_diagnostics(sim_offload, diag0_array, sim_offload->hdf5_out);
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

    return 0;

}

int bmc_init_particles(
        int *n,
        particle_state** ps,
        sim_offload_data* sim_offload,
        B_field_data* Bdata,
        real* offload_array
    ) {
    // vacate the phase space to find the phase-space points in the mesh
    // suitable for simulating the BMC scheme

    // init sim data
    sim_data sim;
    sim_init(&sim, sim_offload);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, offload_array);

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
        *n = n_r * n_phi * n_z * n_vpara * n_vperp * N_MC_STEPS;
    } else {
        *n = n_r * n_phi * n_z * n_vr * n_vphi * n_vz * N_MC_STEPS;
    }

    *ps = (particle_state *)malloc(*n * sizeof(particle_state));
    input_particle p_tmp; // tmp particle
    particle_state ps_tmp; // tmp particle

    int i = 0;
    for (int i_r = 0; i_r < n_r; ++i_r) {
        r = (max_r - min_r) * i_r / n_r + min_r;
        for (int i_phi = 0; i_phi < n_phi; ++i_phi) {
            phi = (max_phi - min_phi) * i_phi / n_phi + min_phi;
            for (int i_z = 0; i_z < n_z; ++i_z) {
                z = (max_z - min_z) * i_z / n_z + min_z;

                // if (!wall_2d_inside(r, z, &sim.wall_data.w2d)) {
                //     continue;
                // }

                if (sim_offload->diag_offload_data.dist5D_collect) {                
                    // compute magnetic field
                    for (int i_vpara = 0; i_vpara < n_vpara; ++i_vpara) {
                        vpara = (max_vpara - min_vpara) * i_vpara / n_vpara + min_vpara;
                        for (int i_vperp = 0; i_vperp < n_vperp; ++i_vperp) {
                            vperp = (max_vperp - min_vperp) * i_vperp / n_vperp + min_vperp;
                            // bmc_5D_to_fo(B_dB, r, phi, z, vpara, vperp, &vr, &vphi, &vz);
                            bmc_5D_to_particle_state(Bdata, r, phi, z, vpara, vperp, T1, i, &ps_tmp);

                            if (!ps_tmp.err) {
                                for (int i_mc=0; i_mc<N_MC_STEPS; ++i_mc) {
                                    ps_tmp.id = i;
                                    memcpy(*ps + i, &ps_tmp, sizeof(particle_state));
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
                                bmc_init_fo_particle(&p_tmp, r, phi, z, vr, vphi, vz, T1, i);
                                particle_input_to_state(&p_tmp, &ps_tmp, Bdata);
                            }

                            if (!ps_tmp.err) {
                                for (int i_mc=0; i_mc<N_MC_STEPS; ++i_mc) {
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

    offload_data->unpack_pos = 0;

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