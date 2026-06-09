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
    bfield_nn_init("config_bfield.yaml", 0);
#endif
 real* rnd = (real*) malloc(3*mrk_array_size*sizeof(real));
    float* nn_input = (float*) malloc(3 * mrk_array_size * sizeof(float));
    float* nn_label = (float*) malloc(12 * mrk_array_size * sizeof(float));
    float* nn_output = (float*) malloc(12 * mrk_array_size * sizeof(float));
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size],rnd[0:3*mrk_array_size],nn_input[0:3*mrk_array_size],nn_label[0:12*mrk_array_size],nn_output[0:12*mrk_array_size])
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
    double phi_mod = fmod(p.phi[i], 2.0*M_PI);
    if(phi_mod < 0.0) phi_mod += 2.0*M_PI;
    
const float bscale[12] = {
    1.953079e+00f,
    3.799192e+00f,
    5.698458e-01f,
    4.110970e+00f,
    1.753601e+01f,
    1.407675e+01f,
    5.582170e-01f,
    3.261646e-02f,
    3.749450e+00f,
    7.147724e+00f,
    7.610648e-02f,
    4.386748e+00f
};

 nn_input[i*3 + 0] = (float)((p.r[i]   - 6.0) / 2.0);
    nn_input[i*3 + 1] = (float)((phi_mod  - M_PI) / M_PI);
    nn_input[i*3 + 2] = (float)(p.z[i] / 2.0);
    
    nn_label[i*12 + 0]  = (float)p.B_r[i] / bscale[0];
    nn_label[i*12 + 1]  = (float)p.B_r_dr[i] / bscale[1];
    nn_label[i*12 + 2]  = (float)p.B_r_dphi[i] / bscale[2];
    nn_label[i*12 + 3]  = (float)p.B_r_dz[i] / bscale[3];
    
    nn_label[i*12 + 4]  = (float)p.B_phi[i] / bscale[4];
    nn_label[i*12 + 5]  = (float)p.B_phi_dr[i] / bscale[5];
    nn_label[i*12 + 6]  = (float)p.B_phi_dphi[i] / bscale[6];
    nn_label[i*12 + 7]  = (float)p.B_phi_dz[i] / bscale[7];
    
    nn_label[i*12 + 8]  = (float)p.B_z[i] / bscale[8];
    nn_label[i*12 + 9]  = (float)p.B_z_dr[i] / bscale[9];
    nn_label[i*12 + 10] = (float)p.B_z_dphi[i] / bscale[10];
    nn_label[i*12 + 11] = (float)p.B_z_dz[i] / bscale[11];
  }
float maxabs0 = 0.0f;
float maxabs1 = 0.0f;
float maxabs2 = 0.0f;
float maxabs3 = 0.0f;
float maxabs4 = 0.0f;
float maxabs5 = 0.0f;
float maxabs6 = 0.0f;
float maxabs7 = 0.0f;
float maxabs8 = 0.0f;
float maxabs9 = 0.0f;
float maxabs10 = 0.0f;
float maxabs11 = 0.0f;

#pragma acc parallel loop present(nn_label[0:12*p.n_mrk]) \
    reduction(max:maxabs0,maxabs1,maxabs2,maxabs3,maxabs4,maxabs5,maxabs6,maxabs7,maxabs8,maxabs9,maxabs10,maxabs11)
for(int i = 0; i < p.n_mrk; i++) {
    float v0  = fabsf(nn_label[i*12 + 0]);
    float v1  = fabsf(nn_label[i*12 + 1]);
    float v2  = fabsf(nn_label[i*12 + 2]);
    float v3  = fabsf(nn_label[i*12 + 3]);
    float v4  = fabsf(nn_label[i*12 + 4]);
    float v5  = fabsf(nn_label[i*12 + 5]);
    float v6  = fabsf(nn_label[i*12 + 6]);
    float v7  = fabsf(nn_label[i*12 + 7]);
    float v8  = fabsf(nn_label[i*12 + 8]);
    float v9  = fabsf(nn_label[i*12 + 9]);
    float v10 = fabsf(nn_label[i*12 + 10]);
    float v11 = fabsf(nn_label[i*12 + 11]);

    if(v0  > maxabs0)  maxabs0  = v0;
    if(v1  > maxabs1)  maxabs1  = v1;
    if(v2  > maxabs2)  maxabs2  = v2;
    if(v3  > maxabs3)  maxabs3  = v3;
    if(v4  > maxabs4)  maxabs4  = v4;
    if(v5  > maxabs5)  maxabs5  = v5;
    if(v6  > maxabs6)  maxabs6  = v6;
    if(v7  > maxabs7)  maxabs7  = v7;
    if(v8  > maxabs8)  maxabs8  = v8;
    if(v9  > maxabs9)  maxabs9  = v9;
    if(v10 > maxabs10) maxabs10 = v10;
    if(v11 > maxabs11) maxabs11 = v11;
}

static int printed_scale = 0;
if(!printed_scale) {
    printf("\n=== BFIELD LABEL SCALE ===\n");
    printf("0  %.6e\n", maxabs0);
    printf("1  %.6e\n", maxabs1);
    printf("2  %.6e\n", maxabs2);
    printf("3  %.6e\n", maxabs3);
    printf("4  %.6e\n", maxabs4);
    printf("5  %.6e\n", maxabs5);
    printf("6  %.6e\n", maxabs6);
    printf("7  %.6e\n", maxabs7);
    printf("8  %.6e\n", maxabs8);
    printf("9  %.6e\n", maxabs9);
    printf("10 %.6e\n", maxabs10);
    printf("11 %.6e\n", maxabs11);
    printf("==========================\n\n");
    printed_scale = 1;
}
#pragma acc host_data use_device(nn_input, nn_label, nn_output)
{
    bfield_nn_train_gpu(nn_input, nn_label, nn_output, p.n_mrk);
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

float abs_d0_sum = 0.0f, abs_d0_max = 0.0f;
float abs_d1_sum = 0.0f, abs_d1_max = 0.0f;
float abs_d2_sum = 0.0f, abs_d2_max = 0.0f;
float abs_d3_sum = 0.0f, abs_d3_max = 0.0f;
float abs_d4_sum = 0.0f, abs_d4_max = 0.0f;
float abs_d5_sum = 0.0f, abs_d5_max = 0.0f;
float abs_d6_sum = 0.0f, abs_d6_max = 0.0f;
float abs_d7_sum = 0.0f, abs_d7_max = 0.0f;
float abs_d8_sum = 0.0f, abs_d8_max = 0.0f;

#pragma acc parallel loop \
present(nn_label[0:12*p.n_mrk], nn_output[0:12*p.n_mrk]) \
reduction(+:err_Bnorm_sum,abs_Br_sum,abs_Bphi_sum,abs_Bz_sum,n_eval, \
             abs_d0_sum,abs_d1_sum,abs_d2_sum,abs_d3_sum,abs_d4_sum, \
             abs_d5_sum,abs_d6_sum,abs_d7_sum,abs_d8_sum) \
reduction(max:err_Bnorm_max,abs_Br_max,abs_Bphi_max,abs_Bz_max, \
              abs_d0_max,abs_d1_max,abs_d2_max,abs_d3_max,abs_d4_max, \
              abs_d5_max,abs_d6_max,abs_d7_max,abs_d8_max)
    for(int i = 0; i < p.n_mrk; i++) {

        float Br_true   = nn_label[i*12 + 0];
        float Bphi_true = nn_label[i*12 + 4];
        float Bz_true   = nn_label[i*12 + 8];

        float Br_pred   = nn_output[i*12 + 0];
        float Bphi_pred = nn_output[i*12 + 4];
        float Bz_pred   = nn_output[i*12 + 8];

        float eps = 1.0e-8f;

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
	
	if(eB > err_Bnorm_max)
	  err_Bnorm_max = eB;
	
	abs_Br_sum   += abs_Br;
	abs_Bphi_sum += abs_Bphi;
	abs_Bz_sum   += abs_Bz;
	
	if(abs_Br > abs_Br_max)
	  abs_Br_max = abs_Br;
	
	if(abs_Bphi > abs_Bphi_max)
	  abs_Bphi_max = abs_Bphi;
	
	if(abs_Bz > abs_Bz_max)
	  abs_Bz_max = abs_Bz;
	
	float ad0 = fabsf(nn_output[i*12 + 1]  - nn_label[i*12 + 1]);
	float ad1 = fabsf(nn_output[i*12 + 2]  - nn_label[i*12 + 2]);
	float ad2 = fabsf(nn_output[i*12 + 3]  - nn_label[i*12 + 3]);
	
	float ad3 = fabsf(nn_output[i*12 + 5]  - nn_label[i*12 + 5]);
	float ad4 = fabsf(nn_output[i*12 + 6]  - nn_label[i*12 + 6]);
	float ad5 = fabsf(nn_output[i*12 + 7]  - nn_label[i*12 + 7]);
	
	float ad6 = fabsf(nn_output[i*12 + 9]  - nn_label[i*12 + 9]);
	float ad7 = fabsf(nn_output[i*12 + 10] - nn_label[i*12 + 10]);
	float ad8 = fabsf(nn_output[i*12 + 11] - nn_label[i*12 + 11]);
	
	abs_d0_sum += ad0; if(ad0 > abs_d0_max) abs_d0_max = ad0;
	abs_d1_sum += ad1; if(ad1 > abs_d1_max) abs_d1_max = ad1;
	abs_d2_sum += ad2; if(ad2 > abs_d2_max) abs_d2_max = ad2;
	
	abs_d3_sum += ad3; if(ad3 > abs_d3_max) abs_d3_max = ad3;
	abs_d4_sum += ad4; if(ad4 > abs_d4_max) abs_d4_max = ad4;
	abs_d5_sum += ad5; if(ad5 > abs_d5_max) abs_d5_max = ad5;
	
	abs_d6_sum += ad6; if(ad6 > abs_d6_max) abs_d6_max = ad6;
	abs_d7_sum += ad7; if(ad7 > abs_d7_max) abs_d7_max = ad7;
	abs_d8_sum += ad8; if(ad8 > abs_d8_max) abs_d8_max = ad8;
	
	n_eval++;
    }

    if(n_eval > 0) {
      printf("[bfield_nn validation] step %d\n", bfield_nn_get_step());
      printf("  |B| relative mean %.6e max %.6e\n",
	     err_Bnorm_sum/n_eval, err_Bnorm_max);
      printf("  abs Br   mean %.6e max %.6e\n",
	     abs_Br_sum/n_eval, abs_Br_max);
      printf("  abs Bphi mean %.6e max %.6e\n",
	     abs_Bphi_sum/n_eval, abs_Bphi_max);
      printf("  abs Bz   mean %.6e max %.6e\n",
	     abs_Bz_sum/n_eval, abs_Bz_max);
      printf("  derivatives abs errors:\n");
      printf("    dBr/dR    mean %.6e max %.6e\n", abs_d0_sum/n_eval, abs_d0_max);
      printf("    dBr/dphi  mean %.6e max %.6e\n", abs_d1_sum/n_eval, abs_d1_max);
      printf("    dBr/dZ    mean %.6e max %.6e\n", abs_d2_sum/n_eval, abs_d2_max);
      
      printf("    dBphi/dR   mean %.6e max %.6e\n", abs_d3_sum/n_eval, abs_d3_max);
      printf("    dBphi/dphi mean %.6e max %.6e\n", abs_d4_sum/n_eval, abs_d4_max);
      printf("    dBphi/dZ   mean %.6e max %.6e\n", abs_d5_sum/n_eval, abs_d5_max);
      
      printf("    dBz/dR    mean %.6e max %.6e\n", abs_d6_sum/n_eval, abs_d6_max);
      printf("    dBz/dphi  mean %.6e max %.6e\n", abs_d7_sum/n_eval, abs_d7_max);
      printf("    dBz/dZ    mean %.6e max %.6e\n", abs_d8_sum/n_eval, abs_d8_max);
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
    free(nn_label);
    free(nn_output);
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


