/**
 * @file simulate_gc_adaptive.c
 * @brief Simulate guiding centers using adaptive time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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
#include "../boozer.h"
#include "../mhd.h"
#include "../rfof.h"
#include "../plasma.h"
#include "simulate_gc_adaptive.h"
#include "step/step_gc_cashkarp.h"
#include "mccc/mccc.h"
#include "mccc/mccc_wiener.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */

/**
 * @brief Simulates guiding centers using adaptive time-step
 *
 * The simulation includes:
 * - orbit-following with Cash-Karp method
 * - Coulomb collisions with Milstein method
 *
 * The simulation is carried until all marker have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
 *
 * @param pq particles to be simulated
 * @param sim simulation data
 *
 */
void simulate_gc_adaptive(particle_queue* pq, sim_data* sim) {

    /* Wiener arrays needed for the adaptive time step */
    mccc_wienarr wienarr[NSIMD];

    /* Current time step, suggestions for the next time step and next time
     * step                                                                */
    real hin[NSIMD]       __memalign__;
    real hout_orb[NSIMD]  __memalign__;
    real hout_col[NSIMD]  __memalign__;
    real hout_rfof[NSIMD] __memalign__;
    real hnext[NSIMD]     __memalign__;

    /* Flag indicateing whether a new marker was initialized */
    int cycle[NSIMD]     __memalign__;

    real tol_col = sim->ada_tol_clmbcol;
    real tol_orb = sim->ada_tol_orbfol;

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states

    rfof_marker rfof_mrk; // RFOF specific data

    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    if(sim->enable_icrh) {
        rfof_set_up(&rfof_mrk, &sim->rfof_data);
    }

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(cycle[i] > 0) {
            /* Determine initial time-step */
            hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
            if(sim->enable_clmbcol) {
                /* Allocate array storing the Wiener processes */
                mccc_wiener_initialize(&(wienarr[i]), p.time[i]);
            }
        }
    }

    cputime_last = A5_WTIME;

    /* MAIN SIMULATION LOOP
     * - Store current state
     * - Integrate motion due to bacgkround EM-field (orbit-following)
     * - Integrate scattering due to Coulomb collisions
     * - Check whether time step was accepted
     *   - NO:  revert to initial state and ignore the end of the loop
     *          (except CPU_TIME_MAX end condition if this is implemented)
     *   - YES: update particle time, clean redundant Wiener processes, and
     *          proceed
     * - Check for end condition(s)
     * - Update diagnostics
     */
    while(n_running > 0) {

        /* Store marker states in case time step will be rejected */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            particle_copy_gc(&p, i, &p0, i);
            hout_orb[i]  = DUMMY_TIMESTEP_VAL;
            hout_col[i]  = DUMMY_TIMESTEP_VAL;
            hout_rfof[i] = DUMMY_TIMESTEP_VAL;
            hnext[i]     = DUMMY_TIMESTEP_VAL;
        }

        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Cash-Karp method for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_gc_cashkarp_mhd(
                    &p, hin, hout_orb, tol_orb, &sim->B_data, &sim->E_data,
                    &sim->boozer_data, &sim->mhd_data, sim->enable_aldforce);
            }
            else {
                step_gc_cashkarp(
                    &p, hin, hout_orb, tol_orb, &sim->B_data, &sim->E_data,
                    sim->enable_aldforce);
            }
            /* Check whether time step was rejected */
            #pragma omp simd
            for(int i = 0; i < NSIMD; i++) {
                /* Switch sign of the time-step again if it was reverted earlier
                */
                if(sim->reverse_time) {
                    hout_orb[i] = -hout_orb[i];
                    hin[i]      = -hin[i];
                }
                if(p.running[i] && hout_orb[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout_orb[i];
                }
            }
        }

        /* Milstein method for collisions */
        if(sim->enable_clmbcol) {
            real rnd[5*NSIMD];
            random_normal_simd(&sim->random_data, 5*NSIMD, rnd);
            mccc_gc_milstein(&p, hin, hout_col, tol_col, wienarr, &sim->B_data,
                             &sim->plasma_data, &sim->mccc_data, rnd);

            /* Check whether time step was rejected */
            #pragma omp simd
            for(int i = 0; i < NSIMD; i++) {
                if(p.running[i] && hout_col[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout_col[i];
                }
            }
        }

        /* Performs the ICRH kick if in resonance. */
        if(sim->enable_icrh) {
            rfof_resonance_check_and_kick_gc(
                &p, hin, hout_rfof, &rfof_mrk, &sim->rfof_data, &sim->B_data);

            /* Check whether time step was rejected */
            #pragma omp simd
            for(int i = 0; i < NSIMD; i++) {
                if(p.running[i] && hout_rfof[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout_rfof[i];
                }
            }
        }

        /**********************************************************************/

        cputime = A5_WTIME;
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(p.id[i] > 0 && !p.err[i]) {
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

                /* Retrieve marker states in case time step was rejected      */
                if(hnext[i] < 0) {
                    particle_copy_gc(&p0, i, &p, i);

                    hin[i] = -hnext[i];
                }
                if(p.running[i]){

                    /* Advance time (if time step was accepted) and determine
                       next time step */
                    if(hnext[i] < 0){
                        /* Time step was rejected, use the suggestion given by
                           integrator */
                        hin[i] = -hnext[i];
                    }
                    else {
                        p.time[i] += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) )
                            * hin[i];
                        p.mileage[i] += hin[i];

                        if(hnext[i] > hout_orb[i]) {
                            /* Use time step suggested by the orbit-following
                               integrator */
                            hnext[i] = hout_orb[i];
                        }
                        if(hnext[i] > hout_col[i]) {
                            /* Use time step suggested by the collision
                               integrator */
                            hnext[i] = hout_col[i];
                        }
                        if(hnext[i] > hout_rfof[i]) {
                            /* Use time step suggested by RFOF */
                            hnext[i] = hout_rfof[i];
                        }
                        if(hnext[i] == 1.0) {
                            /* Time step is unchanged (happens when no physics
                               are enabled) */
                            hnext[i] = hin[i];
                        }
                        hin[i] = hnext[i];
                        if(sim->enable_clmbcol) {
                            /* Clear wiener processes */
                            mccc_wiener_clean(&(wienarr[i]), p.time[i]);
                        }
                    }

                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(&sim->diag_data, &sim->B_data, &p, &p0);

        /* Update number of running particles */
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step for new particles */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
                if(sim->enable_clmbcol) {
                    /* Re-allocate array storing the Wiener processes */
                    mccc_wiener_initialize(&(wienarr[i]), p.time[i]);
                }
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
 * The returned time step is either directly user-defined, 1/100th of collision
 * frequency or user-defined fraction of gyro-motion.
 *
 * @param sim pointer to simulation data struct
 * @param p SIMD array of markers
 * @param i index of marker for which time step is assessed
 *
 * @return Calculated time step
 */
real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i) {
    /* Just use some large value if no physics are defined */
    real h = DUMMY_TIMESTEP_VAL;

    /* Value defined directly by user */
    if(sim->fix_usrdef_use) {
        h =  sim->fix_usrdef_val;
    }
    else {
        /* Value calculated from gyrotime */
        if(sim->enable_orbfol) {
            real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
            real gyrotime = CONST_2PI /
                phys_gyrofreq_ppar(p->mass[i], p->charge[i], p->mu[i],
                                   p->ppar[i], Bnorm);
            if(h > gyrotime) {
                h = gyrotime;
            }
        }

        /* Value calculated from collision frequency */
        if(sim->enable_clmbcol) {
            real nu = 1;
            /*mccc_collfreq_gc(p, &sim->B_data, &sim->plasma_data,
                sim->coldata, &nu, i); */

            /* Only small angle collisions so divide this by 100 */
            real colltime = 1/(100*nu);
            if(h > colltime) {h=colltime;}
        }
    }
    return h;
}
