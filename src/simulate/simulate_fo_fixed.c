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
    bfield_nn_init("config_bfield_B.yaml", 0);
#endif
    real* rnd = (real*) malloc(3*mrk_array_size*sizeof(real));
    float* nn_input     = (float*) malloc(3 * mrk_array_size * sizeof(float));
    
    float* nn_label_B   = (float*) malloc(3 * mrk_array_size * sizeof(float));
    float* nn_output_B  = (float*) malloc(3 * mrk_array_size * sizeof(float));
    
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size],rnd[0:3*mrk_array_size],nn_input[0:3*mrk_array_size],nn_label_B[0:3*mrk_array_size],nn_output_B[0:3*mrk_array_size])
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

GPU_PARALLEL_LOOP_ALL_LEVELS
for(int i = 0; i < p.n_mrk; i++) {

    double phi_mod = fmod(p.phi[i], 2.0 * M_PI);
    if(phi_mod < 0.0) {
        phi_mod += 2.0 * M_PI;
    }

    const float bscale_B[3] = {
        1.953079e+00f,
        1.753601e+01f,
        3.749450e+00f
    };

    nn_input[i*3 + 0] = (float)((p.r[i]  - 6.0) / 2.0);
    nn_input[i*3 + 1] = (float)((phi_mod - M_PI) / M_PI);
    nn_input[i*3 + 2] = (float)(p.z[i] / 2.0);

    nn_label_B[i*3 + 0] = (float)p.B_r[i]   / bscale_B[0];
    nn_label_B[i*3 + 1] = (float)p.B_phi[i] / bscale_B[1];
    nn_label_B[i*3 + 2] = (float)p.B_z[i]   / bscale_B[2];
}
 
float rmin = 1e30f, rmax = -1e30f;
float zmin = 1e30f, zmax = -1e30f;

#pragma acc parallel loop reduction(min:rmin,zmin) \
                          reduction(max:rmax,zmax)
for(int i = 0; i < p.n_mrk; i++) {

    float r = (float)p.r[i];
    float z = (float)p.z[i];

    if(r < rmin) rmin = r;
    if(r > rmax) rmax = r;

    if(z < zmin) zmin = z;
    if(z > zmax) zmax = z;
}

static int printed_bounds = 0;

if(!printed_bounds) {
    printf("\n=== INPUT DOMAIN ===\n");
    printf("R : [%f , %f]\n", rmin, rmax);
    printf("Z : [%f , %f]\n", zmin, zmax);
    printf("====================\n\n");

    printed_bounds = 1;
}

#pragma acc host_data use_device(nn_input, nn_label_B, nn_output_B)
{
    bfield_nn_train_gpu(
        nn_input,
        nn_label_B,
        nn_output_B,
        p.n_mrk
    );
}
 
if (bfield_nn_validation_ready()) {

    float err_Bnorm_sum = 0.0f;
    float err_Bnorm_max = 0.0f;

    float abs_Br_sum   = 0.0f;
    float abs_Bphi_sum = 0.0f;
    float abs_Bz_sum   = 0.0f;

    float abs_Br_max   = 0.0f;
    float abs_Bphi_max = 0.0f;
    float abs_Bz_max   = 0.0f;
    
    int n_eval = 0;

    #pragma acc parallel loop \
    present(nn_label_B[0:3*p.n_mrk], nn_output_B[0:3*p.n_mrk]) \
      reduction(+:err_Bnorm_sum,abs_Br_sum,abs_Bphi_sum,abs_Bz_sum,n_eval) \
      reduction(max:err_Bnorm_max,abs_Br_max,abs_Bphi_max,abs_Bz_max)
    for(int i = 0; i < p.n_mrk; i++) {

        float Br_true   = nn_label_B[i*3 + 0];
        float Bphi_true = nn_label_B[i*3 + 1];
        float Bz_true   = nn_label_B[i*3 + 2];

        float Br_pred   = nn_output_B[i*3 + 0];
        float Bphi_pred = nn_output_B[i*3 + 1];
        float Bz_pred   = nn_output_B[i*3 + 2];

        float dBr   = Br_pred   - Br_true;
        float dBphi = Bphi_pred - Bphi_true;
        float dBz   = Bz_pred   - Bz_true;

        float Berr_norm = sqrtf(dBr*dBr + dBphi*dBphi + dBz*dBz);
        float Btrue_norm = sqrtf(
            Br_true*Br_true +
            Bphi_true*Bphi_true +
            Bz_true*Bz_true
        );

        float eB = Berr_norm / (Btrue_norm + 1.0e-8f);

        float abs_Br   = fabsf(dBr);
        float abs_Bphi = fabsf(dBphi);
        float abs_Bz   = fabsf(dBz);


 err_Bnorm_sum += eB;
        if(eB > err_Bnorm_max) {
            err_Bnorm_max = eB;
        }

        abs_Br_sum   += abs_Br;
        abs_Bphi_sum += abs_Bphi;
        abs_Bz_sum   += abs_Bz;

        if(abs_Br > abs_Br_max) {
            abs_Br_max = abs_Br;
        }

        if(abs_Bphi > abs_Bphi_max) {
            abs_Bphi_max = abs_Bphi;
        }

        if(abs_Bz > abs_Bz_max) {
            abs_Bz_max = abs_Bz;
        }

        n_eval++;
    }

    if(n_eval > 0) {
        printf("[bfield_nn validation] step %d\n", bfield_nn_get_step());
        printf("  |B| relative mean %.6e max %.6e\n",
               err_Bnorm_sum / n_eval, err_Bnorm_max);
 printf("  abs Br   mean %.6e max %.6e\n",
               abs_Br_sum / n_eval, abs_Br_max);
        printf("  abs Bphi mean %.6e max %.6e\n",
               abs_Bphi_sum / n_eval, abs_Bphi_max);
        printf("  abs Bz   mean %.6e max %.6e\n",
               abs_Bz_sum / n_eval, abs_Bz_max);
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
    free(nn_input);
    free(nn_label_B);
    free(nn_output_B);
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


