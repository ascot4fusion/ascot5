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
#include "utils/mathlib.h"
#include "consts.h"
#include "utils/physlib.h"
#include "simulate.h"
#include "data/marker.h"
#include "data/wall.h"
#include "data/diag.h"
#include "datatypes.h"
#include "data/bfield.h"
#include "data/efield.h"
#include "data/rfof.h"
#include "data/plasma.h"
#include "orbit_following.h"
#include "coulomb_collisions.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_fixed_inidt(Simulation* sim, MarkerGuidingCenter* p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */


void simulate_gc_fixed(Simulation* sim, MarkerQueue* pq, size_t vector_size) {
    size_t* cycle = (size_t*) malloc(vector_size*sizeof(size_t));
    real* hin = (real*) malloc(vector_size*sizeof(real));
    real* hin_default = (real*) malloc(vector_size*sizeof(real));
    real* hnext_recom = (real*) malloc(vector_size*sizeof(real));
    real* hout_rfof = (real*) malloc(vector_size*sizeof(real));


    real cputime, cputime_last; // Global cpu time: recent and previous record

    MarkerGuidingCenter p;  // This array holds current states
    MarkerGuidingCenter p0; // This array stores previous states
    marker_allocate_gc(&p, vector_size);
    marker_allocate_gc(&p0, vector_size);

    rfof_marker rfof_mrk; // RFOF specific data

    /* Init dummy markers */
    for(size_t i=0; i< vector_size; i++) {
        p.id[i] = 0;
        p.running[i] = 0;
        hout_rfof[i] = DUMMY_TIMESTEP_VAL;
        hnext_recom[i] = DUMMY_TIMESTEP_VAL;
    }

    /* Initialize running particles */
    size_t n_running = marker_cycle_gc(pq, &p, &sim->bfield, cycle);

    if(sim->options->enable_icrh) {
        rfof_set_up(&rfof_mrk, &sim->rfof);
    }

    /* Determine simulation time-step */
    #pragma omp simd
    for(size_t i = 0; i < vector_size; i++) {
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
    real* rnd = (real*) malloc(5*vector_size*sizeof(real));
    while(n_running > 0) {

        /* Store marker states */
        #pragma omp simd
        for(size_t i = 0; i < vector_size; i++) {
            marker_copy_gc(&p, i, &p0, i);
        }

        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        #pragma omp simd
        for(size_t i = 0; i < vector_size; i++) {
            if(sim->options->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* RK4 method for orbit-following */
        if(sim->options->enable_orbit_following) {
            if(sim->options->enable_mhd) {
                step_gc_rk4_mhd(
                    &p, hin, &sim->bfield, &sim->efield, sim->boozer,
                    &sim->mhd, sim->options->enable_aldforce);
            }
            else {
                step_gc_rk4(&p, hin, &sim->bfield, &sim->efield,
                            sim->options->enable_aldforce);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        #pragma omp simd
        for(size_t i = 0; i < vector_size; i++) {
            if(sim->options->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama method for collisions */
        if(sim->options->enable_coulomb_collisions) {
            random_normal_simd(sim->random_data, 5*vector_size, rnd);
            mccc_gc_euler(&p, hin, &sim->bfield, &sim->plasma,
                          sim->mccc_data, rnd);
        }

        /* Performs the ICRH kick if in resonance. */
        if(sim->options->enable_icrh) {
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
        for(size_t i = 0; i < vector_size; i++) {
            if(hnext_recom[i] < 0) {
                /* Screwed up big time (negative time-step only when RFOF
                    rejected) */
                marker_copy_gc(&p0, i, &p, i);
                hin[i] = -hnext_recom[i];
            }
            if(p.running[i]) {
                if(hnext_recom[i] < 0) {
                    // unsuccessful step, only reset the recommendation
                    hnext_recom[i] = DUMMY_TIMESTEP_VAL;
                } else {
                    // The step was successful
                    p.time[i]    += ( 1.0 - 2.0 * ( sim->options->reverse_time > 0 ) ) * hin[i];
                    p.mileage[i] += hin[i];
                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(sim->diagnostics, sim->options, &sim->bfield, &p, &p0);

        /* Update running particles */
        n_running = marker_cycle_gc(pq, &p, &sim->bfield, cycle);

        /* Determine simulation time-step */
        #pragma omp simd
        for(size_t i = 0; i < vector_size; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
                if(sim->options->enable_icrh) {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
            }
        }

    }

    /* All markers simulated! */
    free(cycle);
    free(hin);
    free(hnext_recom);
    free(hin_default);
    free(hout_rfof);
    free(rnd);

    /* Deallocate rfof structs */
    if(sim->options->enable_icrh) {
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
real simulate_gc_fixed_inidt(Simulation* sim, MarkerGuidingCenter* p, int i) {
    real h;

    /* Value defined directly by user */
    if(sim->options->use_explicit_fixedstep) {
        h = sim->options->explicit_fixedstep;
    }
    else {
        /* Value calculated from gyrotime */
        real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
        real gyrotime = CONST_2PI /
            phys_gyrofreq_ppar(p->mass[i], p->charge[i],
                               p->mu[i], p->ppar[i], Bnorm);
        h = gyrotime/sim->options->gyrodefined_fixedstep;
    }

    return h;
}
