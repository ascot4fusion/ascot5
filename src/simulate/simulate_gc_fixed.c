/**
 * @file simulate_gc_fixed.c
 * @brief Simulate guiding centers using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "../ascot5.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"
#include "../physlib.h"
#include "../simulate.h"
#include "../particle.h"
#include "../wall.h"
#include "../diag.h"
#include "../B_field.h"
#include "../E_field.h"
#include "../plasma.h"
#include "simulate_gc_fixed.h"
#include "step/step_gc_rk4.h"
#include "mccc/mccc.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);

/**
 * @brief Simulates guiding centers using fixed time-step
 *
 * The simulation includes:
 * - orbit-following with RK4 method
 * - Coulomb collisions with Euler-Maruyama method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined either directly or as a fraction of
 * gyrotime.
 *
 * @param pq particles to be simulated
 * @param sim simulation data
 */
void simulate_gc_fixed(particle_queue* pq, sim_data* sim, int mrk_array_size) {
    int* cycle = (int*) malloc(mrk_array_size*sizeof(int)); /* Flag indigating whether a new marker
                                                              was initialized */
    real* hin = (real*) malloc(mrk_array_size*sizeof(real));/* Time step given as an input into the
                                                                integrators. Almost always default.*/
    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states
    particle_allocate_gc(&p, mrk_array_size);
    particle_allocate_gc(&p0, mrk_array_size);
    
    /* Init dummy markers */
    for(int i=0; i< mrk_array_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < mrk_array_size; i++) {
        if(cycle[i] > 0) {
            hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
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
    particle_offload_gc(&p);
    particle_offload_gc(&p0);
    real* rnd = (real*) malloc(5*mrk_array_size*sizeof(real));
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size], rnd[0:5*mrk_array_size])
    while(n_running > 0) {

        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            particle_copy_gc(&p, i, &p0, i);
        }

        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* RK4 method for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_gc_rk4_mhd(&p, hin, &sim->B_data, &sim->E_data,
                                &sim->boozer_data, &sim->mhd_data);
            }
            else {
                step_gc_rk4(&p, hin, &sim->B_data, &sim->E_data);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama method for collisions */
        if(sim->enable_clmbcol) {
            random_normal_simd(&sim->random_data, 5*p.n_mrk, rnd);
            mccc_gc_euler(&p, hin, &sim->B_data, &sim->plasma_data,
                          &sim->mccc_data, rnd);
        }

        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(p.running[i]) {
                p.time[i]    += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) ) * hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(&sim->diag_data, &sim->B_data, &p, &p0);

        /* Update running particles */
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for(int i = 0; i < p.n_mrk; i++)
        {
            if(p.running[i] > 0) n_running++;
        }
#else
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);
#endif
	
#ifndef GPU
        /* Determine simulation time-step */
        #pragma omp simd
        for(int i = 0; i < p.n_mrk; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
            }
        }
#endif

    }

    /* All markers simulated! */
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim[0:1])
    particle_onload_gc(&p);
    n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);
#endif
    free(cycle);
    free(hin);
    free(rnd);

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
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i) {
    real h;

    /* Value defined directly by user */
    if(sim->fix_usrdef_use) {
        h = sim->fix_usrdef_val;
    }
    else {
        /* Value calculated from gyrotime */
        real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
        real gyrotime = CONST_2PI /
            phys_gyrofreq_ppar(p->mass[i], p->charge[i],
                               p->mu[i], p->ppar[i], Bnorm);
        h = gyrotime/sim->fix_gyrodef_nstep;
    }

    return h;
}
