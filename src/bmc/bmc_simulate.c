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
        real* offload_array, int* int_offload_array,
        real h, int n_rk4_subcycles
    ) {

    /**************************************************************************/
    /* 1. Input offload data is unpacked and initialized by calling           */
    /*    respective init functions.                                          */
    /*                                                                        */
    /**************************************************************************/
    sim_data sim;
    sim_init(&sim, sim_offload);

    real* ptr; int* ptrint;
    offload_unpack(offload_data, offload_array,
                   sim_offload->B_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    B_field_init(&sim.B_data, &sim_offload->B_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->E_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    E_field_init(&sim.E_data, &sim_offload->E_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->plasma_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    plasma_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->neutral_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    neutral_init(&sim.neutral_data, &sim_offload->neutral_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->wall_offload_data.offload_array_length,
                   int_offload_array,
                   sim_offload->wall_offload_data.int_offload_array_length,
                   &ptr, &ptrint);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr, ptrint);

    random_init(&sim.random_data, time(NULL));

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
            bmc_step_deterministic(p + i_simd, h_rk4, &sim.B_data, &sim.E_data, &sim.plasma_data,
                                sim.enable_clmbcol, sim.enable_mhd,
                                &sim.random_data, &sim.mccc_data, &sim.boozer_data, &sim.mhd_data);
            bmc_check_simd_particle_wallhit(p+i_simd, &p0, &(sim.wall_data));
        }
    }

    // copy the result of RK4 to N_HERMITE different copy of the particle.
    // compute the hermite weights for each copy of the particle
    copy_particles_simd_to_coll_simd(n_simd_particles, n_hermite_knots, p, p_coll);
    /* Euler-Maruyama method for collisions */
    if(sim.enable_clmbcol) {
        for (int i_simd = 0; i_simd < n_coll_simd_particles; i_simd++) {
            memcpy(&p0, p_coll + i_simd, sizeof(particle_simd_gc));
            bmc_step_stochastic(p_coll + i_simd, h_coll, &sim.B_data, &sim.E_data, &sim.plasma_data, &sim.random_data, &sim.mccc_data);
            bmc_check_simd_particle_wallhit(p_coll+i_simd, &p0, &(sim.wall_data));
        }
    }
}

void init_particles_coll_simd_hermite(int n_simd_particles, int n_hermite_knots,
        particle_simd_gc* p_coll
    ) {
    real hermiteK[10];
    real hermiteW[10];

    if (n_hermite_knots == 5) {
        real ktmp[5] = {-2.856970, -1.355626, 0.000000, 1.355626, 2.856970};
        real wtmp[5] = {0.028218, 0.556662, 1.336868, 0.556662, 0.028218};
        memcpy(hermiteK, ktmp, 5*sizeof(int));
        memcpy(hermiteW, wtmp, 5*sizeof(int));
    } else if (n_hermite_knots == 10) {
        real ktmp[10] = {-4.859463, -3.581823, -2.484326, -1.465989, -0.484936, 0.484936, 1.465989, 2.484326, 3.581823, 4.859463};
        real wtmp[10] = {0.000011, 0.001900, 0.047906, 0.339607, 0.863890, 0.863890, 0.339607, 0.047906, 0.001900, 0.000011};
        memcpy(hermiteK, ktmp, 10*sizeof(int));
        memcpy(hermiteW, wtmp, 10*sizeof(int));
    }
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

    // TODO add hermite knots for n != {1,5}
    real hermiteK[5] = {-2.856970, -1.355626, 0.000000, 1.355626, 2.856970};
    real hermiteW[5] = {0.028218, 0.556662, 1.336868, 0.556662, 0.028218};

    for (int i = 0; i < 5; i++) {
        hermiteW[i] /= PI2E0_5;
    }

    if (n_hermite_knots == 1)
    {
        hermiteK[0] = 0;
        hermiteW[0] = 1;
    }

    int i_coll = 0;
    for (int i=0; i < n_simd_particles; i++) {
        for (int j=0; j < NSIMD; j++) {
            for (int k = 0; k < n_hermite_knots; k++) {
                particle_copy_gc(&p[i], j, &p_coll[i_coll / NSIMD], i_coll % NSIMD);
                p_coll[i_coll / NSIMD].err[i_coll % NSIMD] = p[i].err[j];

                // update hermite knots and weights
                p_coll[i_coll / NSIMD].hermite_knots[i_coll % NSIMD] = hermiteK[k];
                p_coll[i_coll / NSIMD].hermite_weights[i_coll % NSIMD] = hermiteW[k];

                i_coll++;
            }
        }
    }
    // pad dummy particles in the last SIMD element
    if ((i_coll % NSIMD) != 0) {
        for (int i = (i_coll % NSIMD); i < NSIMD; i++) {
            p_coll[i_coll / NSIMD].id[i] = -1;
        }
    }
}

void fmc_simulation(
        particle_state* ps,
        sim_offload_data* sim,
        offload_package* offload_data,
        real* offload_array,
        int* int_offload_array,
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
                    offload_array, int_offload_array, diag_offload_array_host);
                *host_end = omp_get_wtime();
            }
        #endif
    }
}


/**
 * @brief Integrate a guiding center step for a struct of markers with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD
 * markers with RK4 simultaneously using SIMD instructions. All arrays in the
 * function are of NSIMD length so vectorization can be performed directly
 * without gather and scatter operations.
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void bmc_step_deterministic(particle_simd_gc *p, real *h, B_field_data *Bdata,
                 E_field_data *Edata, plasma_data* pdata, int enable_clmbcol, int enable_mhd,
                 random_data* rdata, mccc_data* mdata, boozer_data* boozer, mhd_data* mhd)
{

    real rnd[5 * NSIMD];
    random_normal_simd(rdata, 5 * NSIMD, rnd);
    int n_species = plasma_get_n_species(pdata);
    const real *qb = plasma_get_species_charge(pdata);
    const real *mb = plasma_get_species_mass(pdata);

    if (enable_mhd) {
        step_gc_rk4_mhd(p, h, Bdata, Edata, boozer, mhd);
    }
    else {
        step_gc_rk4(p, h, Bdata, Edata);
    }

    if (!enable_clmbcol)
        return;

    /// DETERMINISTIC PART OF COLLISION OPERATOR
    int i;
    #pragma omp simd aligned(h : 64)
    for (i = 0; i < NSIMD; i++)
    {
        if (p->running[i]) {
            a5err errflag = 0;

            /* Initial (R,z) position and magnetic field are needed for later */
            real Brpz[3] = {p->B_r[i], p->B_phi[i], p->B_z[i]};
            real Bnorm = math_norm(Brpz);
            real Bxyz[3];
            math_vec_rpz2xyz(Brpz, Bxyz, p->phi[i]);
            real R0 = p->r[i];
            real z0 = p->z[i];

            /* Move guiding center to (x, y, z, vnorm, xi) coordinates */
            real vin, pin, xiin, Xin_xyz[3];
            pin = physlib_gc_p(p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            vin = physlib_vnorm_pnorm(p->mass[i], pin);
            Xin_xyz[0] = p->r[i] * cos(p->phi[i]);
            Xin_xyz[1] = p->r[i] * sin(p->phi[i]);
            Xin_xyz[2] = p->z[i];

            /* Evaluate plasma density and temperature */
            real nb[MAX_SPECIES], Tb[MAX_SPECIES];
            if (!errflag)
            {
                errflag = plasma_eval_densandtemp(nb, Tb, p->rho[i],
                                                  p->r[i], p->phi[i], p->z[i],
                                                  p->time[i], pdata);
            }

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, p->mass[i], p->charge[i], vin,
                            n_species, mb, qb, nb, Tb);

            /* Evaluate collision coefficients and sum them for each *
             * species                                               */
            real gyrofreq = phys_gyrofreq_pnorm(p->mass[i], p->charge[i],
                                                pin, Bnorm);
            real K = 0, Dpara = 0, nu = 0, DX = 0;
            for (int j = 0; j < n_species; j++)
            {
                real vb = sqrt(2 * Tb[j] / mb[j]);
                real x = vin / vb;
                real mufun[3];
                mccc_coefs_mufun(mufun, x, mdata); // eq. 2.83 PhD Hirvijoki

                real Qb = mccc_coefs_Q(p->mass[i], p->charge[i], mb[j],
                                       qb[j], nb[j], vb, clogab[j],
                                       mufun[0]);
                real Dparab = mccc_coefs_Dpara(p->mass[i], p->charge[i], vin,
                                               qb[j], nb[j], vb, clogab[j],
                                               mufun[0]);
                real Dperpb = mccc_coefs_Dperp(p->mass[i], p->charge[i], vin,
                                               qb[j], nb[j], vb, clogab[j],
                                               mufun[1]);
                real dDparab = mccc_coefs_dDpara(p->mass[i], p->charge[i], vin,
                                                 qb[j], nb[j], vb, clogab[j],
                                                 mufun[0], mufun[2]);

                K += mccc_coefs_K(vin, Dparab, dDparab, Qb);
                Dpara += Dparab;
                nu += mccc_coefs_nu(vin, Dperpb); // eq.41
                DX += mccc_coefs_DX(xiin, Dparab, Dperpb, gyrofreq);
            }

            /* Evaluate collisions */
            real sdt = sqrt(h[i]);
            real dW[5];
            dW[0] = 0;
            dW[1] = 0;
            dW[2] = 0;
            dW[3] = 0;
            dW[4] = 0;

            real bhat[3];
            math_unit(Bxyz, bhat);

            real kc1 = sqrt(2 * DX);
            real kc2 = math_dot(bhat, dW);

            real vout, xiout, Xout_xyz[3];
            Xout_xyz[0] = Xin_xyz[0] + kc1 * (dW[0] - kc2 * bhat[0]);
            Xout_xyz[1] = Xin_xyz[1] + kc1 * (dW[1] - kc2 * bhat[1]);
            Xout_xyz[2] = Xin_xyz[2] + kc1 * (dW[2] - kc2 * bhat[2]);
            vout = vin + K * h[i] + sqrt(2 * Dpara) * dW[3];
            xiout = xiin - xiin * nu * h[i] + sqrt((1 - xiin * xiin) * nu) * dW[4];

            /* Enforce boundary conditions */
            real cutoff = MCCC_CUTOFF * sqrt(Tb[0] / p->mass[i]);
            if (vout < cutoff)
            {
                vout = 2 * cutoff - vout;
            }

            if (fabs(xiout) > 1)
            {
                xiout = ((xiout > 0) - (xiout < 0)) * (2 - fabs(xiout));
            }

            /* Remove energy or pitch change or spatial diffusion from the    *
             * results if that is requested                                   */
            if (!mdata->include_energy)
            {
                vout = vin;
            }
            if (!mdata->include_pitch)
            {
                xiout = xiin;
            }
            if (!mdata->include_gcdiff)
            {
                Xout_xyz[0] = Xin_xyz[0];
                Xout_xyz[1] = Xin_xyz[1];
                Xout_xyz[2] = Xin_xyz[2];
            }
            real pout = physlib_pnorm_vnorm(p->mass[i], vout);

            /* Back to cylindrical coordinates */
            real Xout_rpz[3];
            math_xyz2rpz(Xout_xyz, Xout_rpz);

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real B_dB[15], psi[1], rho[1];
            if (!errflag)
            {
                errflag = B_field_eval_B_dB(B_dB, Xout_rpz[0], Xout_rpz[1],
                                            Xout_rpz[2], p->time[i] + h[i],
                                            Bdata);
            }
            if (!errflag)
            {
                errflag = B_field_eval_psi(psi, Xout_rpz[0], Xout_rpz[1],
                                           Xout_rpz[2], p->time[i] + h[i],
                                           Bdata);
            }
            if (!errflag)
            {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
            }

            if (!errflag)
            {
                /* Update marker coordinates at the new position */
                p->B_r[i] = B_dB[0];
                p->B_r_dr[i] = B_dB[1];
                p->B_r_dphi[i] = B_dB[2];
                p->B_r_dz[i] = B_dB[3];

                p->B_phi[i] = B_dB[4];
                p->B_phi_dr[i] = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i] = B_dB[7];

                p->B_z[i] = B_dB[8];
                p->B_z_dr[i] = B_dB[9];
                p->B_z_dphi[i] = B_dB[10];
                p->B_z_dz[i] = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);

                p->r[i] = Xout_rpz[0];
                p->z[i] = Xout_rpz[2];
                p->ppar[i] = physlib_gc_ppar(pout, xiout);
                p->mu[i] = physlib_gc_mu(p->mass[i], pout, xiout, Bnorm);

                /* Evaluate phi and theta angles so that they are cumulative */
                real rz[2];
                B_field_get_axis_rz(rz, Bdata, p->phi[i]);
                p->theta[i] += atan2((R0 - rz[0]) * (p->z[i] - rz[1]) - (z0 - rz[1]) * (p->r[i] - rz[0]),
                                     (R0 - rz[0]) * (p->r[i] - rz[0]) + (z0 - rz[1]) * (p->z[i] - rz[1]));
                p->phi[i] += atan2(Xin_xyz[0] * Xout_xyz[1] - Xin_xyz[1] * Xout_xyz[0],
                                   Xin_xyz[0] * Xout_xyz[0] + Xin_xyz[1] * Xout_xyz[1]);
            }

            /* Error handling */
            if (errflag)
            {
                p->err[i] = errflag;
                p->running[i] = 0;
            }
        }
    }
}/**
 * @brief Integrate a guiding center step for a struct of markers with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD
 * markers with RK4 simultaneously using SIMD instructions. All arrays in the
 * function are of NSIMD length so vectorization can be performed directly
 * without gather and scatter operations.
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void bmc_step_stochastic(particle_simd_gc *p, real *h, B_field_data *Bdata,
                 E_field_data *Edata, plasma_data* pdata, random_data* rdata, mccc_data* mdata)
{

    real rnd[5 * NSIMD];
    random_normal_simd(rdata, 5 * NSIMD, rnd);
    int n_species = plasma_get_n_species(pdata);
    const real *qb = plasma_get_species_charge(pdata);
    const real *mb = plasma_get_species_mass(pdata);

    int i;
/* Following loop will be executed simultaneously for all i */
#pragma omp simd aligned(h : 64)
    for (i = 0; i < NSIMD; i++)
    {
        if (p->running[i])
        {
            a5err errflag = 0;

            /// STOCHASTIC PART OF COLLISION OPERATOR
            /* Initial (R,z) position and magnetic field are needed for later */
            real Brpz[3] = {p->B_r[i], p->B_phi[i], p->B_z[i]};
            real Bnorm = math_norm(Brpz);
            real Bxyz[3];
            math_vec_rpz2xyz(Brpz, Bxyz, p->phi[i]);
            real R0 = p->r[i];
            real z0 = p->z[i];

            /* Move guiding center to (x, y, z, vnorm, xi) coordinates */
            real vin, pin, xiin, Xin_xyz[3];
            pin = physlib_gc_p(p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            vin = physlib_vnorm_pnorm(p->mass[i], pin);
            Xin_xyz[0] = p->r[i] * cos(p->phi[i]);
            Xin_xyz[1] = p->r[i] * sin(p->phi[i]);
            Xin_xyz[2] = p->z[i];

            /* Evaluate plasma density and temperature */
            real nb[MAX_SPECIES], Tb[MAX_SPECIES];
            if (!errflag)
            {
                errflag = plasma_eval_densandtemp(nb, Tb, p->rho[i],
                                                  p->r[i], p->phi[i], p->z[i],
                                                  p->time[i], pdata);
            }

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, p->mass[i], p->charge[i], vin,
                            n_species, mb, qb, nb, Tb);

            /* Evaluate collision coefficients and sum them for each *
             * species                                               */
            real gyrofreq = phys_gyrofreq_pnorm(p->mass[i], p->charge[i],
                                                pin, Bnorm);
            real K = 0, Dpara = 0, nu = 0, DX = 0;
            for (int j = 0; j < n_species; j++)
            {
                real vb = sqrt(2 * Tb[j] / mb[j]);
                real x = vin / vb;
                real mufun[3];
                mccc_coefs_mufun(mufun, x, mdata); // eq. 2.83 PhD Hirvijoki

                real Qb = mccc_coefs_Q(p->mass[i], p->charge[i], mb[j],
                                       qb[j], nb[j], vb, clogab[j],
                                       mufun[0]);
                real Dparab = mccc_coefs_Dpara(p->mass[i], p->charge[i], vin,
                                               qb[j], nb[j], vb, clogab[j],
                                               mufun[0]);
                real Dperpb = mccc_coefs_Dperp(p->mass[i], p->charge[i], vin,
                                               qb[j], nb[j], vb, clogab[j],
                                               mufun[1]);
                real dDparab = mccc_coefs_dDpara(p->mass[i], p->charge[i], vin,
                                                 qb[j], nb[j], vb, clogab[j],
                                                 mufun[0], mufun[2]);

                K += mccc_coefs_K(vin, Dparab, dDparab, Qb);
                Dpara += Dparab;
                nu += mccc_coefs_nu(vin, Dperpb); // eq.41
                DX += mccc_coefs_DX(xiin, Dparab, Dperpb, gyrofreq);
            }

            /* Evaluate collisions */
            real sdt = sqrt(h[i]);
            real dW[5];
            dW[0] = 0;
            dW[1] = 0;
            dW[2] = 0;
            dW[3] = 0;
            dW[4] = 0;
            
            if (p->hermite_weights[i] > 0)
            {
                dW[4] = sdt * p->hermite_knots[i];
            }
            else
            {
                printf("warning: no hermite knots for the particle\n");
                dW[4] = sdt * rnd[4 * NSIMD + i]; // For v
            }

            real bhat[3];
            math_unit(Bxyz, bhat);

            real kc1 = sqrt(2 * DX);
            real kc2 = math_dot(bhat, dW);

            real vout, xiout, Xout_xyz[3];
            Xout_xyz[0] = Xin_xyz[0] + kc1 * dW[0];
            Xout_xyz[1] = Xin_xyz[1] + kc1 * dW[1];
            Xout_xyz[2] = Xin_xyz[2] + kc1 * dW[2];
            vout = vin + sqrt(2 * Dpara) * dW[3];
            xiout = xiin + sqrt((1 - xiin * xiin) * nu) * dW[4];

            /* Enforce boundary conditions */
            real cutoff = MCCC_CUTOFF * sqrt(Tb[0] / p->mass[i]);
            if (vout < cutoff)
            {
                vout = 2 * cutoff - vout;
            }

            if (fabs(xiout) > 1)
            {
                xiout = ((xiout > 0) - (xiout < 0)) * (2 - fabs(xiout));
            }

            /* Remove energy or pitch change or spatial diffusion from the    *
             * results if that is requested                                   */
            if (!mdata->include_energy)
            {
                vout = vin;
            }
            if (!mdata->include_pitch)
            {
                xiout = xiin;
            }
            if (!mdata->include_gcdiff)
            {
                Xout_xyz[0] = Xin_xyz[0];
                Xout_xyz[1] = Xin_xyz[1];
                Xout_xyz[2] = Xin_xyz[2];
            }
            real pout = physlib_pnorm_vnorm(p->mass[i], vout);

            /* Back to cylindrical coordinates */
            real Xout_rpz[3];
            math_xyz2rpz(Xout_xyz, Xout_rpz);

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real B_dB[15], psi[1], rho[1];
            if (!errflag)
            {
                errflag = B_field_eval_B_dB(B_dB, Xout_rpz[0], Xout_rpz[1],
                                            Xout_rpz[2], p->time[i] + h[i],
                                            Bdata);
            }
            if (!errflag)
            {
                errflag = B_field_eval_psi(psi, Xout_rpz[0], Xout_rpz[1],
                                           Xout_rpz[2], p->time[i] + h[i],
                                           Bdata);
            }
            if (!errflag)
            {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
            }

            if (!errflag)
            {
                /* Update marker coordinates at the new position */
                p->B_r[i] = B_dB[0];
                p->B_r_dr[i] = B_dB[1];
                p->B_r_dphi[i] = B_dB[2];
                p->B_r_dz[i] = B_dB[3];

                p->B_phi[i] = B_dB[4];
                p->B_phi_dr[i] = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i] = B_dB[7];

                p->B_z[i] = B_dB[8];
                p->B_z_dr[i] = B_dB[9];
                p->B_z_dphi[i] = B_dB[10];
                p->B_z_dz[i] = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);

                p->r[i] = Xout_rpz[0];
                p->z[i] = Xout_rpz[2];
                p->ppar[i] = physlib_gc_ppar(pout, xiout);
                p->mu[i] = physlib_gc_mu(p->mass[i], pout, xiout, Bnorm);

                /* Evaluate phi and theta angles so that they are cumulative */
                real rz[2];
                B_field_get_axis_rz(rz, Bdata, p->phi[i]);
                p->theta[i] += atan2((R0 - rz[0]) * (p->z[i] - rz[1]) - (z0 - rz[1]) * (p->r[i] - rz[0]),
                                     (R0 - rz[0]) * (p->r[i] - rz[0]) + (z0 - rz[1]) * (p->z[i] - rz[1]));
                p->phi[i] += atan2(Xin_xyz[0] * Xout_xyz[1] - Xin_xyz[1] * Xout_xyz[0],
                                   Xin_xyz[0] * Xout_xyz[0] + Xin_xyz[1] * Xout_xyz[1]);
            }

            /* Error handling */
            if (errflag)
            {
                p->err[i] = errflag;
                p->running[i] = 0;
            }
        }
    }
}