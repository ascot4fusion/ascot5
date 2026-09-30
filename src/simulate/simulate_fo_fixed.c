/**
 * @file simulate_fo_fixed.c
 * @brief Simulate particles using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "../ascot5.h"
#include "../physlib.h"
#include "../simulate.h"
#include "../particle.h"
#include "../wall.h"
#include "../diag.h"
#include "../B_field.h"
#include "../E_field.h"
#include "../plasma.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"
#include "simulate_fo_fixed.h"
#include "step/step_fo_vpa.h"
#include "mccc/mccc.h"
#include "atomic.h"
#ifdef USE_TORCHFORT_BFIELD
#include "../bfield_nn.h"
#endif

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i);

/**
 * @brief Simulates particles using fixed time-step
 *
 * The simulation includes:
 * - orbit-following with Volume-Preserving Algorithm
 * - Coulomb collisions with Euler-Maruyama method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined: either a directly given fixed value
 * or a given fraction of gyrotime.
 *
 * @param pq particles to be simulated
 * @param sim simulation data struct
 * @param mrk_array_size size of particle arrays
 */
void simulate_fo_fixed(particle_queue* pq, sim_data* sim, int mrk_array_size) {
    // Indicates whether a new marker was initialized
    int* cycle = (int*) malloc(mrk_array_size*sizeof(int));
    // Time-step
    real* hin = (real*) malloc(mrk_array_size*sizeof(real));

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states
    particle_allocate_fo(&p, mrk_array_size);
    particle_allocate_fo(&p0, mrk_array_size);

    /* Init dummy markers */
    for(int i = 0; i < mrk_array_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < mrk_array_size; i++) {
        if(cycle[i] > 0) {
            hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
        }
    }

    cputime_last = A5_WTIME;

    /* MAIN SIMULATION LOOP
     * - Store current state
     * - Integrate motion due to background EM-field (orbit-following)
     * - Integrate scattering due to Coulomb collisions
     * - Advance time
     * - Check for end condition(s)
     * - Update diagnostics
     */
    particle_offload_fo(&p);
    particle_offload_fo(&p0);
#ifdef USE_TORCHFORT_BFIELD
    bfield_nn_init("config_bfield_B_dB.yaml", 0);
#endif
    real* rnd = (real*) malloc(3*mrk_array_size*sizeof(real));
#ifdef USE_TORCHFORT_BFIELD
    const int nn_batch_capacity = (mrk_array_size < 2048) ? mrk_array_size : 2048;
    float* nn_input = (float*) malloc(3 * nn_batch_capacity * sizeof(float));
    float* nn_label = (float*) malloc(13 * nn_batch_capacity * sizeof(float));
    float* nn_output = (float*) malloc(13 * nn_batch_capacity * sizeof(float));
    real nn_rmin = p.r[0], nn_rmax = p.r[0];
    real nn_zmin = p.z[0], nn_zmax = p.z[0];
    for(int i = 1; i < p.n_mrk; i++) {
        if(p.r[i] < nn_rmin) nn_rmin = p.r[i];
        if(p.r[i] > nn_rmax) nn_rmax = p.r[i];
        if(p.z[i] < nn_zmin) nn_zmin = p.z[i];
        if(p.z[i] > nn_zmax) nn_zmax = p.z[i];
    }
    printf("[bfield_nn] replay domain R=[%g,%g] Z=[%g,%g]\n", nn_rmin, nn_rmax, nn_zmin, nn_zmax);
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size],rnd[0:3*mrk_array_size],nn_input[0:3*nn_batch_capacity],nn_label[0:13*nn_batch_capacity],nn_output[0:13*nn_batch_capacity])
#else
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size],rnd[0:3*mrk_array_size])
#endif
    while(n_running > 0) {
        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            particle_copy_fo(&p, i, &p0, i);
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_fo_vpa_mhd(
                    &p, hin, &sim->B_data, &sim->E_data, &sim->boozer_data,
                    &sim->mhd_data, sim->enable_aldforce);
            }
            else {
                step_fo_vpa(&p, hin, &sim->B_data, &sim->E_data,
                            sim->enable_aldforce);
#ifdef USE_TORCHFORT_BFIELD
                int nn_count = (p.n_mrk < nn_batch_capacity) ? p.n_mrk : nn_batch_capacity;
                if(nn_count > 0) {
                    const int step = bfield_nn_get_step();
                    const float bscale[3] = {1.953079f, 17.53601f, 3.749450f};
                    const float gscale[3][3] = {
                        {3.55015f, 0.52316f, 4.50280f},
                        {0.53320f, 0.06598f, 0.00752f},
                        {2.94440f, 0.12928f, 2.03438f}
                    };
                    const int nn_marker_count = nn_count / 2;
                    GPU_PARALLEL_LOOP_ALL_LEVELS
                    for(int i = 0; i < nn_count; i++) {
                        real r, phi, z;
                        real B_dB[15];
                        if(i < nn_marker_count) {
                            int idx = (int)(((long long)step * nn_marker_count + i) % p.n_mrk);
                            r = p.r[idx];
                            phi = fmod(p.phi[idx], 2.0 * M_PI);
                            if(phi < 0.0) phi += 2.0 * M_PI;
                            z = p.z[idx];
                            B_dB[0] = p.B_r[idx];
                            B_dB[1] = p.B_r_dr[idx];
                            B_dB[2] = p.B_r_dphi[idx];
                            B_dB[3] = p.B_r_dz[idx];
                            B_dB[4] = p.B_phi[idx];
                            B_dB[5] = p.B_phi_dr[idx];
                            B_dB[6] = p.B_phi_dphi[idx];
                            B_dB[7] = p.B_phi_dz[idx];
                            B_dB[8] = p.B_z[idx];
                            B_dB[9] = p.B_z_dr[idx];
                            B_dB[10] = p.B_z_dphi[idx];
                            B_dB[11] = p.B_z_dz[idx];
                        } else {
                            int j = i - nn_marker_count;
                            real ur = fmod((j + 1) * 0.618033989 + (step + 1) * 0.754877666, 1.0);
                            real uz = fmod((j + 1) * 0.414213562 + (step + 1) * 0.569840291, 1.0);
                            real up = fmod((j + 1) * 0.732050808 + (step + 1) * 0.438579021, 1.0);
                            r = nn_rmin + ur * (nn_rmax - nn_rmin);
                            z = nn_zmin + uz * (nn_zmax - nn_zmin);
                            phi = 2.0 * M_PI * up;
                            (void)B_field_eval_B_dB(B_dB, r, phi, z, 0.0, &sim->B_data);
                        }

                        nn_input[i*3 + 0] = (float)((r - 6.0) / 2.0);
                        nn_input[i*3 + 1] = (float)((phi - M_PI) / M_PI);
                        nn_input[i*3 + 2] = (float)(z / 2.0);

                        nn_label[i*13 + 0] = (float)(B_dB[0] / bscale[0]);
                        nn_label[i*13 + 1] = (float)(2.0 * B_dB[1] / (bscale[0] * gscale[0][0]));
                        nn_label[i*13 + 2] = (float)(M_PI * B_dB[2] / (bscale[0] * gscale[0][1]));
                        nn_label[i*13 + 3] = (float)(2.0 * B_dB[3] / (bscale[0] * gscale[0][2]));
                        nn_label[i*13 + 4] = (float)(B_dB[4] / bscale[1]);
                        nn_label[i*13 + 5] = (float)(2.0 * B_dB[5] / (bscale[1] * gscale[1][0]));
                        nn_label[i*13 + 6] = (float)(M_PI * B_dB[6] / (bscale[1] * gscale[1][1]));
                        nn_label[i*13 + 7] = (float)(2.0 * B_dB[7] / (bscale[1] * gscale[1][2]));
                        nn_label[i*13 + 8] = (float)(B_dB[8] / bscale[2]);
                        nn_label[i*13 + 9] = (float)(2.0 * B_dB[9] / (bscale[2] * gscale[2][0]));
                        nn_label[i*13 + 10] = (float)(M_PI * B_dB[10] / (bscale[2] * gscale[2][1]));
                        nn_label[i*13 + 11] = (float)(2.0 * B_dB[11] / (bscale[2] * gscale[2][2]));
                        real divB = B_dB[1] + B_dB[0] / r + B_dB[6] / r + B_dB[11];
                        nn_label[i*13 + 12] = (float)(0.5477225575 * divB / 3.03743887);
                    }


                    #pragma acc host_data use_device(nn_input, nn_label, nn_output)
                    {
                        bfield_nn_train_gpu(nn_input, nn_label, nn_output, nn_count);
                    }

                    if(bfield_nn_validation_ready()) {
                        float b_rel_sum = 0.0f, b_rel_max = 0.0f;
                        float grad_err_sum = 0.0f, grad_err_max = 0.0f;
                        float div_err_sq_sum = 0.0f;
                        #pragma acc parallel loop present(nn_label[0:13*nn_count], nn_output[0:13*nn_count]) \
                            reduction(+:b_rel_sum,grad_err_sum,div_err_sq_sum) \
                            reduction(max:b_rel_max,grad_err_max)
                        for(int i = 0; i < nn_count; i++) {
                            float db0 = nn_output[i*13 + 0] - nn_label[i*13 + 0];
                            float db1 = nn_output[i*13 + 4] - nn_label[i*13 + 4];
                            float db2 = nn_output[i*13 + 8] - nn_label[i*13 + 8];
                            float bnorm = sqrtf(db0*db0 + db1*db1 + db2*db2);
                            float btrue = sqrtf(nn_label[i*13 + 0]*nn_label[i*13 + 0]
                                              + nn_label[i*13 + 4]*nn_label[i*13 + 4]
                                              + nn_label[i*13 + 8]*nn_label[i*13 + 8]);
                            float brel = bnorm / (btrue + 1.0e-8f);
                            if(brel > b_rel_max) b_rel_max = brel;
                            b_rel_sum += brel;

                            float gerr = 0.0f;
                            for(int k = 0; k < 12; k++) {
                                if(k == 0 || k == 4 || k == 8) continue;
                                float e = nn_output[i*13 + k] - nn_label[i*13 + k];
                                gerr += e*e;
                            }
                            gerr = sqrtf(gerr / 9.0f);
                            if(gerr > grad_err_max) grad_err_max = gerr;
                            grad_err_sum += gerr;
                            float de = (nn_output[i*13 + 12] - nn_label[i*13 + 12]) / 0.5477225575f;
                            div_err_sq_sum += de*de;
                        }
                        printf("[bfield_nn validation] step %d, markers %d\n", bfield_nn_get_step(), nn_count);
                        printf("  |B| relative mean %.6e max %.6e\n", b_rel_sum/nn_count, b_rel_max);
                        printf("  scaled dB gradient error mean %.6e max %.6e\n", grad_err_sum/nn_count, grad_err_max);
                        printf("  divB error RMS / term scale %.6e\n", sqrtf(div_err_sq_sum/nn_count));
                    }
                }
#endif
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if(sim->enable_clmbcol) {
            random_normal_simd(&sim->random_data, 3*p.n_mrk, rnd);
            mccc_fo_euler(&p, hin, &sim->plasma_data, &sim->mccc_data, rnd);
        }
        /* Atomic reactions */
        if(sim->enable_atomic) {
            atomic_fo(&p, hin, &sim->plasma_data, &sim->neutral_data,
                      &sim->random_data, &sim->asigma_data);
        }
        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(p.running[i]){
                p.time[i]    += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) ) * hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_fo(&p, &p0, sim);

        /* Update diagnostics */
        if(!(sim->record_mode)) {
            /* Record particle coordinates */
            diag_update_fo(&sim->diag_data, &sim->B_data, &p, &p0);
        }
        else {
            /* Instead of particle coordinates we record guiding center */

            // Dummy guiding centers
            particle_simd_gc gc_f;
            particle_simd_gc gc_i;

            /* Particle to guiding center transformation */
            #pragma omp simd
            for(int i=0; i<p.n_mrk; i++) {
                if(p.running[i]) {
                    particle_fo_to_gc(&p, i, &gc_f, &sim->B_data);
                }
                else {
                    gc_f.id[i] = p.id[i];
                    gc_f.running[i] = 0;
                }
                if(p0.running[i]) {
                    particle_fo_to_gc(&p0, i, &gc_i, &sim->B_data);
                }
                else {
                    gc_i.id[i] = p0.id[i];
                    gc_i.running[i] = 0;
                }
            }
            diag_update_gc(&sim->diag_data, &sim->B_data, &gc_f, &gc_i);
        }

        /* Update running particles */
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for(int i = 0; i < p.n_mrk; i++)
        {
            if(p.running[i] > 0) n_running++;
        }
#else
        n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif
#ifndef GPU
        /* Determine simulation time-step for new particles */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
            }
        }
#endif
    }
    /* All markers simulated! */
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim[0:1])
    particle_onload_fo(&p);
    n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif
    free(cycle);
    free(hin);
    free(rnd);
#ifdef USE_TORCHFORT_BFIELD
    bfield_nn_finalize();
#endif
#ifdef USE_TORCHFORT_BFIELD
    free(nn_input);
    free(nn_label);
    free(nn_output);
#endif
}

/**
 * @brief Calculates time step value
 *
 * The time step is calculated as a user-defined fraction of gyro time,
 * whose formula accounts for relativity, or an user defined value
 * is used as is depending on simulation options.
 *
 * @param sim pointer to simulation data struct
 * @param p SIMD array of markers
 * @param i index of marker for which time step is assessed
 *
 * @return Calculated time step
 */
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i) {

    real h;

    /* Value defined directly by user */
    if(sim->fix_usrdef_use) {
        h = sim->fix_usrdef_val;
    }
    else {
      /* Value calculated from gyrotime */
        real Bnorm = math_normc( p->B_r[i], p->B_phi[i], p->B_z[i] );
        real pnorm = math_normc( p->p_r[i], p->p_phi[i], p->p_z[i] );
        real gyrotime = CONST_2PI/
            phys_gyrofreq_pnorm(p->mass[i], p->charge[i], pnorm, Bnorm);
        h = gyrotime/sim->fix_gyrodef_nstep;
    }

    return h;
}


