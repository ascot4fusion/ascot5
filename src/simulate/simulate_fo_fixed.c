/**
 * @file simulate_fo_fixed.c
 * @brief Simulate particles using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <immintrin.h>
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

#pragma omp declare target
#pragma omp declare simd uniform(sim)
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i);
#pragma omp end declare target

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
 */
void simulate_fo_fixed(particle_queue* pq, sim_data* sim) {
    int cycle[NSIMD]  __memalign__; // Indicates if a new marker was initialized
    real hin[NSIMD]  __memalign__;  // Time step

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states

    /* Init dummy markers */
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
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
    while(n_running > 0) {

        /* Store marker states */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            particle_copy_fo(&p, i, &p0, i);
        }

        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_fo_vpa_mhd(&p, hin, &sim->B_data, &sim->E_data,
                                &sim->boozer_data, &sim->mhd_data);
            }
            else {
                step_fo_vpa(&p, hin, &sim->B_data, &sim->E_data);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if(sim->enable_clmbcol) {
            mccc_fo_euler(&p, hin, &sim->plasma_data, sim->random_data,
                          &sim->mccc_data);
        }

        /* Atomic reactions */
        if(sim->enable_atomic) {
            atomic_fo(&p, hin, &sim->plasma_data, &sim->neutral_data,
                      &sim->random_data, &sim->asigma_data);
        }

        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
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
            for(int i=0; i<NSIMD; i++) {
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
        n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step for new particles */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
            }
        }
    }

    /* All markers simulated! */

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
