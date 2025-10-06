/**
 * @file simulate_gc_fixed.c
 * @brief Simulate guiding centers using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "defines.h"
#include "endcond.h"
#include "mathlib.h"
#include "consts.h"
#include "physlib.h"
#include "simulate.h"
#include "particle.h"
#include "wall.h"
#include "diag.h"
#include "bfield.h"
#include "efield.h"
#include "rfof.h"
#include "plasma.h"
#include "orbit_following/orbit_following.h"
#include "coulomb_collisions/coulomb_collisions.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */


void simulate_gc_fixed(particle_queue* pq, sim_data* sim) {
    int cycle[NSIMD]      __memalign__; /* Flag indigating whether a new marker
                                           was initialized */
    real hin[NSIMD]       __memalign__;  /* Time step given as an input into the
                                           integrators. Almost always default.*/
    real hin_default[NSIMD] __memalign__; /* The default fixed time step.     */
    real hnext_recom[NSIMD]     __memalign__; /* Next time step, only used to
                                                store the value when RFOF has
                                                rejected a time step.         */
    real hout_rfof[NSIMD] __memalign__; /* The time step that RFOF recommends.
                                            Small positive means that resonance
                                            is close, small negative means that
                                            the step failed because the marker
                                            overshot the resonance and that the
                                            time step should be retaken with a
                                            smaller time step given by the
                                            negative of hout_rfof             */


    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states

    rfof_marker rfof_mrk; // RFOF specific data

    /* Init dummy markers */
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
        hout_rfof[i] = DUMMY_TIMESTEP_VAL;
        hnext_recom[i] = DUMMY_TIMESTEP_VAL;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->bfield, cycle);

    if(sim->params->enable_icrh) {
        rfof_set_up(&rfof_mrk, &sim->rfof);
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
     * - Perform ICRH kick with RFOF if in wave-particle resonance
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
            if(sim->params->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* RK4 method for orbit-following */
        if(sim->params->enable_orbit_following) {
            if(sim->params->enable_mhd) {
                step_gc_rk4_mhd(
                    &p, hin, &sim->bfield, &sim->efield, sim->boozer,
                    &sim->mhd, sim->params->enable_aldforce);
            }
            else {
                step_gc_rk4(&p, hin, &sim->bfield, &sim->efield,
                            sim->params->enable_aldforce);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(sim->params->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama method for collisions */
        if(sim->params->enable_coulomb_collisions) {
            real rnd[5*NSIMD];
            random_normal_simd(sim->random_data, 5*NSIMD, rnd);
            mccc_gc_euler(&p, hin, &sim->bfield, &sim->plasma,
                          sim->mccc_data, rnd);
        }

        /* Performs the ICRH kick if in resonance. */
        if(sim->params->enable_icrh) {
            rfof_resonance_check_and_kick_gc(
                &p, hin, hout_rfof, &rfof_mrk, &sim->rfof, &sim->bfield);

            /* Check whether time step was rejected */
            #pragma omp simd
            for(int i = 0; i < NSIMD; i++) {
                if(p.running[i] && hout_rfof[i] < 0){
                    // Screwed up big time
                    p.running[i] = 0;
                    hnext_recom[i] = hout_rfof[i];  /* Use the smaller time-step
                                                       suggested by RFOF on the
                                                       next round. */
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
            if(hnext_recom[i] < 0) {
                /* Screwed up big time (negative time-step only when RFOF
                    rejected) */
                particle_copy_gc(&p0, i, &p, i);
                hin[i] = -hnext_recom[i];
            }
            if(p.running[i]) {
                if(hnext_recom[i] < 0) {
                    // unsuccessful step, only reset the recommendation
                    hnext_recom[i] = DUMMY_TIMESTEP_VAL;
                } else {
                    // The step was successful
                    p.time[i]    += ( 1.0 - 2.0 * ( sim->params->reverse_time > 0 ) ) * hin[i];
                    p.mileage[i] += hin[i];
                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(sim->diagnostics, sim->params, &sim->bfield, &p, &p0);

        /* Update running particles */
        n_running = particle_cycle_gc(pq, &p, &sim->bfield, cycle);

        /* Determine simulation time-step */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
                if(sim->params->enable_icrh) {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
            }
        }

    }

    /* All markers simulated! */

    /* Deallocate rfof structs */
    if(sim->params->enable_icrh) {
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
    if(sim->params->use_explicit_fixedstep) {
        h = sim->params->explicit_fixedstep;
    }
    else {
        /* Value calculated from gyrotime */
        real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
        real gyrotime = CONST_2PI /
            phys_gyrofreq_ppar(p->mass[i], p->charge[i],
                               p->mu[i], p->ppar[i], Bnorm);
        h = gyrotime/sim->params->gyrodefined_fixedstep;
    }

    return h;
}
