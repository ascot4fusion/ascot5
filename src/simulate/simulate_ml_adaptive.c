/**
 * @file simulate_ml_adaptive.c
 * @brief Simulate magnetic field-lines using adaptive time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>
#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"
#include "../wall.h"
#include "../diag.h"
#include "../B_field.h"
#include "../E_field.h"
#include "simulate_ml_adaptive.h"
#include "step/step_ml_cashkarp.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"

#pragma omp declare target
#pragma omp declare simd uniform(sim)
real simulate_ml_adaptive_inidt(sim_data* sim, particle_simd_ml* p, int i);
#pragma omp end declare target


#define MAGNETIC_FIELD_LINE_INISTEP 1.0e-2 /**< Initial step size in meters   */
#define DUMMY_STEP_VAL 100.0              /**< Dummy orbit step val in meters */

/**
 * @brief Simulates magnetic field-lines using adaptive time-step
 *
 * The simulation includes:
 * - orbit-following with Cash-Karp method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
 *
 * Note that even though we might refer the integration time-step
 * as "time", in reality we are integrating over distance. The time
 * step is therefore step in meters marker orbit is intgerated. Marker does have
 * time in its field, but it is the global time and that is not being changed
 * during the simulation.
 *
 * @param pq field lines to be simulated
 * @param sim simulation data struct
 */
void simulate_ml_adaptive(particle_queue* pq, sim_data* sim) {

    real hin[NSIMD] __memalign__;   // Current time step
    real hout[NSIMD] __memalign__;  // Suggestion for next time step, negative value indicates rejected step
    real hnext[NSIMD] __memalign__; // Next time step
    int cycle[NSIMD] __memalign__;  // Flag indigating whether a new marker was initialized

    real cputime, cputime_last; // Global cpu time: recent and previous record

    real tol = sim->ada_tol_orbfol;
    int i;

    particle_simd_ml p;  // This array holds current states
    particle_simd_ml p0; // This array stores previous states

    for(i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_ml(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(cycle[i] > 0) {
            /* Determine initial time-step */
            hin[i] = simulate_ml_adaptive_inidt(sim, &p, i);
        }
    }

    cputime_last = A5_WTIME;

    /* MAIN SIMULATION LOOP
     * - Store current state
     * - Integrate motion due to bacgkround EM-field (orbit-following)
     * - Check whether time step was accepted
     *   - NO:  revert to initial state and ignore the end of the loop
     *          (except CPU_TIME_MAX end condition if this is implemented)
     *   - YES: update particle time, clean redundant Wiener processes, and proceed
     * - Check for end condition(s)
     * - Update diagnostics
     * - Check for end condition(s)
     */
    while(n_running > 0) {

        /* Store marker states in case time step will be rejected */
        #pragma omp simd
        for(i = 0; i < NSIMD; i++) {
            particle_copy_ml(&p, i, &p0, i);

            hout[i] = DUMMY_STEP_VAL;
            hnext[i] = DUMMY_STEP_VAL;
        }

        /*************************** Physics **********************************/

        /* Cash-Karp method for orbit-following */
        if(sim->enable_orbfol) {

            /* Set time-step negative if tracing backwards in time */
            #pragma omp simd
            for(i = 0; i < NSIMD; i++) {
                if(sim->reverse_time) {
                    hin[i]  = -hin[i];
                }
            }

            if(sim->enable_mhd) {
                step_ml_cashkarp_mhd(&p, hin, hout, tol, &sim->B_data,
                                     &sim->boozer_data, &sim->mhd_data);
            }
            else {
                step_ml_cashkarp(&p, hin, hout, tol, &sim->B_data);
            }

            /* Check whether time step was rejected */
            #pragma omp simd
            for(i = 0; i < NSIMD; i++) {
                /* Switch sign of the time-step again if it was reverted earlier */
                if(sim->reverse_time) {
                    hout[i] = -hout[i];
                    hin[i]  = -hin[i];
                }

                if(p.running[i] && hout[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout[i];
                }
            }
        }

        /**********************************************************************/


        cputime = A5_WTIME;
        #pragma omp simd
        for(i = 0; i < NSIMD; i++) {
            if(!p.err[i]) {
                /* Check other time step limitations */
                if(hnext[i] > 0) {
                    real dphi = fabs(p0.phi[i]-p.phi[i]) / sim->ada_max_dphi;
                    real drho = fabs(p0.rho[i]-p.rho[i]) / sim->ada_max_drho;

                    if(dphi > 1 && dphi > drho) {
                        hnext[i] = -hin[i]/dphi;
                    }
                    else if(drho > 1 && drho > dphi) {
                        hnext[i] = -hin[i]/drho;
                    }
                }

                /* Retrieve marker states in case time step was rejected */
                if(hnext[i] < 0){
                    particle_copy_ml(&p0, i, &p, i);
                }

                /* Update simulation and cpu times */
                if(p.running[i]){

                    /* Advance time (if time step was accepted) and determine next time step */
                    if(hnext[i] < 0){
                        /* Time step was rejected, use the suggestion given by integrator */
                        hin[i] = -hnext[i];
                    }
                    else {
                        /* Mileage measures seconds but hin is in meters */
                        p.mileage[i] += hin[i] / CONST_C;

                        if(hnext[i] > hout[i]) {
                            /* Use time step suggested by the integrator */
                            hnext[i] = hout[i];
                        }
                        else if(hnext[i] == DUMMY_STEP_VAL) {
                            /* Time step is unchanged (happens when no physics are enabled) */
                            hnext[i] = hin[i];
                        }
                        hin[i] = hnext[i];
                    }

                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_ml(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_ml(&sim->diag_data, &p, &p0);

        /* Update running particles */
        n_running = particle_cycle_ml(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step for new particles */
        #pragma omp simd
        for(i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_ml_adaptive_inidt(sim, &p, i);
            }
        }
    }

    /* All markers simulated! */

}

/**
 * @brief Calculates initial time step value
 *
 * The time step value (in units of meters) is defined
 * by MAGNETIC_FIELD_LINE_INISTEP
 *
 * @param sim pointer to simulation data struct
 * @param p SIMD array of markers
 * @param i index of marker for which time step is assessed
 *
 * @return Calculated time step
 */
real simulate_ml_adaptive_inidt(sim_data* sim, particle_simd_ml* p, int i) {

    /* Value calculated from speed of light
     * assuming initial time step is of the order of 1 cm / c*/
    /* Define this with a compiler parameter */

    return MAGNETIC_FIELD_LINE_INISTEP;

}
