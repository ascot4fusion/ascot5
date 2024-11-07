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
#include "../copytogpu.h"
#include "simulate_fo_fixed.h"
#include "step/step_fo_vpa.h"
#include "mccc/mccc.h"
#include "atomic.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i);

DECLARE_TARGET_SIMD_UNIFORM(sim)
void sort_by_key_int_wrapper(int *data, int *value, int N);
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
    int cycle[mrk_array_size];// Indicates whether a new marker was initialized
    real hin[mrk_array_size];// Time step
    real hin2[mrk_array_size];// Time step

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p2;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states
    particle_allocate_fo(&p, mrk_array_size);
    particle_allocate_fo(&p2, mrk_array_size);
    particle_allocate_fo(&p0, mrk_array_size);

    /* Init dummy markers */
    for(int i=0; i< mrk_array_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
	p.initialIndex[i] = i;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < mrk_array_size; i++) {
        if(cycle[i] > 0) {
            hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
	    hin2[i] = hin[i];
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
    particle_simd_fo *p_ptr = &p;
    particle_simd_fo *p0_ptr = &p0;
    particle_simd_fo *p2_ptr = &p2;
    real* hin_ptr = hin;
    real* hin2_ptr = hin2;
    real rnd[3*mrk_array_size];
    int sort_index[mrk_array_size];
    int ps[mrk_array_size];
    for (int i=0;i<mrk_array_size;i++) {
      sort_index[i] = i;
      ps[i] = -1*p_ptr->running[i];
    }

    int step = 0;
    p.n_mrk_ref = p.n_mrk;  
    p2.n_mrk_ref = p.n_mrk;  
    int packing = 0;
    int ipack = 0;
    particle_simd_fo *p_tmp_ptr;
    real* hin_tmp_ptr;
    real rpack_coef = 1.;
    int n_running_ref = n_running;
#ifdef PACK_COEF
    rpack_coef = PACK_COEF;
#endif    
    printf("Pack coefficient = %f",rpack_coef);
    particle_offload_fo(p_ptr);
    particle_offload_fo(p2_ptr);
    particle_offload_fo(p0_ptr);
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size],hin2[0:mrk_array_size], rnd[0:3*mrk_array_size], ps[0:mrk_array_size],sort_index[0:mrk_array_size])
    while(n_running > 0) {
      step++;
      if(step % 10000 == 0) { // print every 100th step
        printf("Current time: %f n_running: %d s\n", step*hin[0], n_running);
      }
        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_running_ref; i++) {
	    if(p_ptr->running[i]){
               particle_copy_fo(p_ptr, i, p0_ptr, i);
	    }
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_running_ref; i++) {
            if(sim->reverse_time) {
	       if(p_ptr->running[i]){
                   hin[i] = -hin[i];
	       }
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_fo_vpa_mhd(p_ptr, hin_ptr, &sim->B_data, &sim->E_data,
                                &sim->boozer_data, &sim->mhd_data);
            }
            else {
                step_fo_vpa(p_ptr, hin_ptr, &sim->B_data, &sim->E_data);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_running_ref; i++) {
            if(sim->reverse_time) {
	       if(p_ptr->running[i]){
                  hin[i]  = -hin[i];
	       }
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if(sim->enable_clmbcol) {
            random_normal_simd(&sim->random_data, 3*p_ptr->n_mrk, rnd);
            mccc_fo_euler(p_ptr, hin_ptr, &sim->plasma_data, &sim->mccc_data, rnd);
        }
        /* Atomic reactions */
        if(sim->enable_atomic) {
            atomic_fo(p_ptr, hin_ptr, &sim->plasma_data, &sim->neutral_data,
                      &sim->random_data, &sim->asigma_data);
        }
        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_running_ref; i++) {
            if(p_ptr->running[i]){
                p_ptr->time[i]    += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) ) * hin[i];
                p_ptr->mileage[i] += hin[i];
                p_ptr->cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_fo(p_ptr, p0_ptr, sim);

        /* Update diagnostics */
        if(!(sim->record_mode)) {
            /* Record particle coordinates */
            diag_update_fo(&sim->diag_data, &sim->B_data, p_ptr, p0_ptr);
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
        /* Update, sort and pack running particles */
#ifdef GPU
	if ((n_running_ref - n_running) > rpack_coef * p_ptr->n_mrk  ) {
	  packing = 1;
	  n_running_ref = n_running;
	}
#pragma acc kernels
	{
	  p_ptr->n_mrk_ref = n_running_ref;
	}
	if (packing == 1) {
	  ipack++;
	  printf("PACKING: n_running = %d => packing number %d\n",n_running,ipack);
	  GPU_PARALLEL_LOOP_ALL_LEVELS
	  for (int i=0;i<n_running_ref;i++) {
	    ps[i] = -1*p_ptr->running[i];
	    sort_index[i] = i;
	  }

#pragma acc host_data use_device(ps,sort_index)
	  {
	    sort_by_key_int_wrapper(ps,sort_index,n_running_ref);
	  }

	  GPU_PARALLEL_LOOP_ALL_LEVELS
	  for(int iloc = 0; iloc < NSIMD; iloc++)
	    {
	      int i = sort_index[iloc];
	      particle_copy_fo(p_ptr, i, p2_ptr, iloc);
	      //hinbis_ptr[iloc] = hin_ptr[i];
	      //use this to circumvent nvhpc compiler bug 
	      //hin_copy_fo(hin_ptr, i, hin2_ptr, iloc);
	    }
	}

        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for(int i = 0; i < p.n_mrk; i++)
        {
            if(p_ptr->running[i] > 0) n_running++;
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
    // Swap pointers to particle arrays
#ifdef GPU
	if (packing == 1) {
	  hin_tmp_ptr = hin_ptr;
	  hin_ptr = hin2_ptr;
	  hin2_ptr = hin_tmp_ptr;
	  p_tmp_ptr = p_ptr;
	  p_ptr = p2_ptr;
	  p2_ptr = p_tmp_ptr;

	  packing = 0;
	}
	if (n_running == 0) {
	  GPU_PARALLEL_LOOP_ALL_LEVELS
	  for(int iloc = 0; iloc < p_ptr->n_mrk; iloc++)
	    {
	      sort_index[iloc] = p_ptr->initialIndex[iloc];
	    }

	  GPU_PARALLEL_LOOP_ALL_LEVELS
          for(int iloc = 0; iloc < p_ptr->n_mrk; iloc++)
	    {
	      int i = sort_index[iloc];
	      particle_copy_fo(p_ptr, iloc, p2_ptr, i);
	      //use this to circumvent nvhpc compiler bug
	      //hin_copy_fo(hin_ptr, iloc, hin2_ptr, i);
	    }
	  p_tmp_ptr = p_ptr;
	  p_ptr = p2_ptr;
	  p2_ptr = p_tmp_ptr;
	}
#endif
    }

    /* All markers simulated! */
#ifdef GPU
    simulate_fo_fixed_copy_from_gpu(sim, p_ptr, p2_ptr);
    n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
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


