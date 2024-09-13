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
 * @param n_queue_size size of particle arrays
 */
void simulate_fo_fixed(particle_queue* pq, sim_data* sim, int n_queue_size) {
    int cycle[n_queue_size];//  __memalign__; // Flag indigating whether a new marker was initialized
    real hin[n_queue_size];//  __memalign__;  // Time step

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states
    particle_to_fo_alloc(&p, n_queue_size);
    particle_to_fo_alloc(&p0, n_queue_size);

    /* Init dummy markers */
    for(int i=0; i< n_queue_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle, n_queue_size);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < n_queue_size; i++) {
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
    int diag_data_field_size = sim->diag_data.diagorb.Nmrk*sim->diag_data.diagorb.Npnt;
    particle_simd_fo *p_ptr=&p;
    particle_simd_fo *p0_ptr=&p0;
    B_field_data* Bdata = &sim->B_data;
    E_field_data* Edata = &sim->E_data;
    particle_loc  p_loc;
    particle_loc_alloc(&p_loc, n_queue_size);
    real rnd[3*n_queue_size];

#ifdef GPU
    simulate_fo_fixed_copy_to_gpu(sim, p_ptr, p0_ptr, Bdata, Edata, &p_loc, hin, rnd, n_queue_size);
#endif    
    while(n_running > 0) {
        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_queue_size; i++) {
	  particle_copy_fo(p_ptr, i, p0_ptr, i);
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_queue_size; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
#ifdef GPU
	      printf("step_fo_vpa_mhd NOT YET PORTED TO GPU");
	      exit(1);
#endif
                step_fo_vpa_mhd(&p, hin, &sim->B_data, &sim->E_data,
                                &sim->boozer_data, &sim->mhd_data);
            }
            else {
	      step_fo_vpa(p_ptr, hin, &sim->B_data, &sim->E_data, n_queue_size);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_queue_size; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if(sim->enable_clmbcol) {
#if !defined(GPU) || defined(RANDOM_LCG)
	  mccc_fo_euler(p_ptr, hin, &sim->plasma_data,
#if defined(RANDOM_LCG)
			&sim->random_data,
#else
			sim->random_data,
#endif			
			&sim->mccc_data,
			rnd, n_queue_size);
#else
	  printf("mccc_fo_euler ported on GPU only for RANDOM_LCG");
	  exit(1);	  
#endif	  
        }
        /* Atomic reactions */
        if(sim->enable_atomic) {
#ifdef GPU
	  printf("atomic_fo NOT YET PORTED TO GPU");
	  exit(1);
#endif
            atomic_fo(p_ptr, hin, &sim->plasma_data, &sim->neutral_data,
                      &sim->random_data, &sim->asigma_data,
                      &sim->enable_atomic);
        }
        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_queue_size; i++) {
            if(p.running[i]){
                p.time[i]    += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) ) * hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_fo(p_ptr, p0_ptr, sim, n_queue_size);

        /* Update diagnostics */
        if(!(sim->record_mode)) {
            /* Record particle coordinates */
	  diag_update_fo(&sim->diag_data, &sim->B_data, p_ptr, p0_ptr, &p_loc, n_queue_size);
        }
        else {
#ifdef GPU
	  printf("particle_fo_to_gc NOT YET PORTED TO GPU");
	  exit(1);
#endif	  
	  /* Instead of particle coordinates we record guiding center */

            // Dummy guiding centers
            particle_simd_gc gc_f;
            particle_simd_gc gc_i;

            /* Particle to guiding center transformation */
            #pragma omp simd
            for(int i=0; i<n_queue_size; i++) {
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
	for(int i = 0; i < n_queue_size; i++)
	  {
	    if(p_ptr->running[i] > 0) n_running++;
	  }
#else
	n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle, n_queue_size);
#endif
#ifndef GPU	
        /* Determine simulation time-step for new particles */
	GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < n_queue_size; i++) {
	  if(cycle[i] > 0)
	    {
	      hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
	    }
        }
#endif	
    }
    /* All markers simulated! */

#ifdef GPU
    simulate_fo_fixed_copy_from_gpu(sim, p_ptr, n_queue_size);
    n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle, n_queue_size);
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


