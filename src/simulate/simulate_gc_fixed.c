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
#include "../rfof.h"
#include "../plasma.h"
#include "simulate_gc_fixed.h"
#include "step/step_gc_rk4.h"
#include "mccc/mccc.h"

#pragma omp declare target
#pragma omp declare simd uniform(sim)
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);
#pragma omp end declare target

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
void simulate_gc_fixed(particle_queue* pq, sim_data* sim) {
    int cycle[NSIMD]      __memalign__; /* Flag indigating whether a new marker
                                           was initialized */
    real hin[NSIMD]       __memalign__;  // Fixed time step (used by default)
    real hin_default[NSIMD] __memalign__;
    //real hnext_recom[NSIMD]     __memalign__; /* Next time step (same as hin almost
     //                                      always, except for when RFOF needs
       //                                    smaller)                           */
    real hout_rfof[NSIMD] __memalign__; /* The time step that RFOF recommends.
                                            Small positive means that resonance
                                            is close, small negative means that
                                            the step failed because the marker
                                            overshot the resonance and that the
                                            time step should be retaken with a
                                            smaller time step given by the
                                            negative of hout_rfof             */

    // Start by using the default
    //for(int i = 0; i<NSIMD; i++) {hnext[i] = hin[i];};

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states

    rfof_marker rfof_mrk; // RFOF specific data

    /* Init dummy markers */
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    if(sim->enable_icrh) {
        rfof_set_up(&rfof_mrk, &sim->rfof_data);
    }

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(cycle[i] > 0) {
            hin_default[i] = simulate_gc_fixed_inidt(sim, &p, i);
            hin[i] = hin_default[i];
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
    while(n_running > 0) {

        /* Store marker states */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            particle_copy_gc(&p, i, &p0, i);
        }

        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
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
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama method for collisions */
        if(sim->enable_clmbcol) {
            real rnd[5*NSIMD];
            random_normal_simd(&sim->random_data, 5*NSIMD, rnd);
            mccc_gc_euler(&p, hin, &sim->B_data, &sim->plasma_data,
                          &sim->mccc_data, rnd);
        }

        /* Performs the ICRH kick if in resonance. */
        if(sim->enable_icrh) {
            rfof_resonance_check_and_kick_gc(
                &p, hin, hout_rfof, &rfof_mrk, &sim->rfof_data, &sim->B_data);

            /* Check whether time step was rejected */
            #pragma omp simd
            for(int i = 0; i < NSIMD; i++) {
                if(p.running[i] && hout_rfof[i] < 0){
                    // Screwed up big time
                    p.running[i] = 0;
                    hin[i] = hout_rfof[i];  /* Use the smaller time-step
                                                 suggested by RFOF on the next
                                                 round. */
                } else if(p.running[i]) {
                    // Everything went better than expected
                    hin[i] = hin_default[i];  // use the original fixed step
                }
            }
        }

        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(p.running[i]) {
                // Check whether the time-step was ok
                if(hin[i] < 0) {
                    /* Screwed up big time (negative time-step only when RFOF
                       rejected) */
                    particle_copy_gc(&p0, i, &p, i);
                    hin[i] = -hin[i];
                } else {
                    // The step was successful
                    p.time[i]    += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) ) * hin[i];
                    p.mileage[i] += hin[i];
                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(&sim->diag_data, &sim->B_data, &p, &p0);

        /* Update running particles */
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
                if(sim->enable_icrh) {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
            }
        }

    }

    /* All markers simulated! */

    /* Deallocate rfof structs */
    if(sim->enable_icrh) {
        rfof_tear_down(&rfof_mrk);
    }
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
