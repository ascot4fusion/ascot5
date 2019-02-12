/**
 * @file simulate_gc_adaptive.c
 * @brief Simulate guiding centers using adaptive time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <immintrin.h>
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
#include "simulate_gc_adaptive.h"
#include "step/step_gc_cashkarp.h"
#include "mccc/mccc.h"
#include "mccc/mccc_wiener.h"

#pragma omp declare target
#pragma omp declare simd uniform(sim)
real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i);
#pragma omp end declare target

#define DUMMY_TIMESTEP_VAL 1.0 // Dummy time step val, just use value large enough not to be encountered in actual simulations

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
    mccc_wienarr* wienarr[NSIMD];
    for(int i=0; i < NSIMD; i++) {
        wienarr[i] = malloc(sizeof(mccc_wienarr));
    }

    real hin[NSIMD] __memalign__;       // Current time step
    real hout_orb[NSIMD] __memalign__;  // Suggestion for next time step by orbit-following, negative value indicates rejected step
    real hout_col[NSIMD]  __memalign__; // Same as above but for collisions
    real hnext[NSIMD] __memalign__;     // Next time step
    int cycle[NSIMD] __memalign__;      // Flag indigating whether a new marker was initialized

    real tol_col = sim->ada_tol_clmbcol;
    real tol_orb = sim->ada_tol_orbfol;
    int i;

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states

    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(cycle[i] > 0) {
            /* Determine initial time-step */
            hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
            if(sim->enable_clmbcol) {
                /* Allocate array storing the Wiener processes */
                mccc_wiener_initialize(wienarr[i],p.time[i]);
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
     *   - YES: update particle time, clean redundant Wiener processes, and proceed
     * - Check for end condition(s)
     * - Update diagnostics
     */
    while(n_running > 0) {
        #pragma omp simd
        for(i = 0; i < NSIMD; i++) {
            /* Store marker states in case time step will be rejected */
            p0.r[i]        = p.r[i];
            p0.phi[i]      = p.phi[i];
            p0.z[i]        = p.z[i];
            p0.vpar[i]     = p.vpar[i];
            p0.mu[i]       = p.mu[i];
            p0.theta[i]    = p.theta[i];

            p0.time[i]       = p.time[i];
            p0.cputime[i]    = p.cputime[i];
            p0.rho[i]        = p.rho[i];
            p0.weight[i]     = p.weight[i];
            p0.cputime[i]    = p.cputime[i];
            p0.rho[i]        = p.rho[i];
            p0.pol[i]        = p.pol[i];

            p0.mass[i]       = p.mass[i];
            p0.charge[i]     = p.charge[i];

            p0.id[i]         = p.id[i];
            p0.running[i]    = p.running[i];
            p0.endcond[i]    = p.endcond[i];
            p0.walltile[i]   = p.walltile[i];

            p0.B_r[i]        = p.B_r[i];
            p0.B_phi[i]      = p.B_phi[i];
            p0.B_z[i]        = p.B_z[i];

            p0.B_r_dr[i]     = p.B_r_dr[i];
            p0.B_r_dphi[i]   = p.B_r_dphi[i];
            p0.B_r_dz[i]     = p.B_r_dz[i];

            p0.B_phi_dr[i]   = p.B_phi_dr[i];
            p0.B_phi_dphi[i] = p.B_phi_dphi[i];
            p0.B_phi_dz[i]   = p.B_phi_dz[i];

            p0.B_z_dr[i]     = p.B_z_dr[i];
            p0.B_z_dphi[i]   = p.B_z_dphi[i];
            p0.B_z_dz[i]     = p.B_z_dz[i];

            hout_orb[i] = DUMMY_TIMESTEP_VAL;
            hout_col[i] = DUMMY_TIMESTEP_VAL;
            hnext[i]    = DUMMY_TIMESTEP_VAL;
        }

        /*************************** Physics **********************************/

        /* Cash-Karp method for orbit-following */
        if(sim->enable_orbfol) {
            step_gc_cashkarp(&p, hin, hout_orb, tol_orb,
                             &sim->B_data, &sim->E_data);

            /* Check whether time step was rejected */
            #pragma omp simd
            for(i = 0; i < NSIMD; i++) {
                if(p.running[i] && hout_orb[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout_orb[i];
                }
            }
        }

        /* Milstein method for collisions */
        if(sim->enable_clmbcol) {
            mccc_gc_milstein(&p, hin, hout_col, tol_col, wienarr, &sim->B_data,
                             &sim->plasma_data, &sim->random_data,
                             &sim->mccc_data);

            /* Check whether time step was rejected */
            #pragma omp simd
            for(i = 0; i < NSIMD; i++) {
                if(p.running[i] && hout_col[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout_col[i];
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
                    p.r[i]        = p0.r[i];
                    p.phi[i]      = p0.phi[i];
                    p.z[i]        = p0.z[i];
                    p.vpar[i]     = p0.vpar[i];
                    p.mu[i]       = p0.mu[i];
                    p.theta[i]    = p0.theta[i];

                    p.time[i]       = p0.time[i];
                    p.rho[i]        = p0.rho[i];
                    p.weight[i]     = p0.weight[i];
                    p.rho[i]        = p0.rho[i];
                    p.pol[i]        = p0.pol[i];

                    p.mass[i]       = p0.mass[i];
                    p.charge[i]     = p0.charge[i];

                    p.running[i]    = p0.running[i];
                    p.endcond[i]    = p0.endcond[i];
                    p.walltile[i]   = p0.walltile[i];

                    p.B_r[i]        = p0.B_r[i];
                    p.B_phi[i]      = p0.B_phi[i];
                    p.B_z[i]        = p0.B_z[i];

                    p.B_r_dr[i]     = p0.B_r_dr[i];
                    p.B_r_dphi[i]   = p0.B_r_dphi[i];
                    p.B_r_dz[i]     = p0.B_r_dz[i];

                    p.B_phi_dr[i]   = p0.B_phi_dr[i];
                    p.B_phi_dphi[i] = p0.B_phi_dphi[i];
                    p.B_phi_dz[i]   = p0.B_phi_dz[i];

                    p.B_z_dr[i]     = p0.B_z_dr[i];
                    p.B_z_dphi[i]   = p0.B_z_dphi[i];
                    p.B_z_dz[i]     = p0.B_z_dz[i];

                    hin[i] = -hnext[i];

                }
                if(p.running[i]){

                    /* Advance time (if time step was accepted) and determine next time step */
                    if(hnext[i] < 0){
                        /* Time step was rejected, use the suggestion given by integrator */
                        hin[i] = -hnext[i];
                    }
                    else {
                        p.time[i] = p.time[i] + hin[i];

                        if(hnext[i] > hout_orb[i]) {
                            /* Use time step suggested by the orbit-following integrator */
                            hnext[i] = hout_orb[i];
                        }
                        if(hnext[i] > hout_col[i]) {
                            /* Use time step suggested by the collision integrator */
                            hnext[i] = hout_col[i];
                        }
                        if(hnext[i] == 1.0) {
                            /* Time step is unchanged (happens when no physics are enabled) */
                            hnext[i] = hin[i];
                        }
                        hin[i] = hnext[i];
                        if(sim->enable_clmbcol) {
                            /* Clear wiener processes */
                            mccc_wiener_clean(wienarr[i], p.time[i]);
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
        diag_update_gc(&sim->diag_data, &p, &p0);

        /* Update number of running particles */
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step for new particles */
        #pragma omp simd
        for(i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
                if(sim->enable_clmbcol) {
                    /* Re-allocate array storing the Wiener processes */
                    mccc_wiener_initialize(wienarr[i],p.time[i]);
                }
            }
        }
    }

    /* All markers simulated! */

}

/**
 * @brief Calculates time step value
 *
 * The returned time step is either directly user-defined, 1/100th of collision frequency
 * or user-defined fraction of gyro-motion.
 *
 * @param p SIMD array of markers
 * @param i index of marker for which time step is assessed
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
                phys_gyrofreq_vpar(p->mass[i], p->charge[i], p->mu[i], p->vpar[i], Bnorm);
            if(h > gyrotime) {
                h = gyrotime;
            }
        }

        /* Value calculated from collision frequency */
        if(sim->enable_clmbcol) {
            real nu = 1;
            //mccc_collfreq_gc(p,&sim->B_data,&sim->plasma_data, sim->coldata,&nu,i);

            /* Only small angle collisions so divide this by 100 */
            real colltime = 1/(100*nu);
            if(h > colltime) {h=colltime;}
        }
    }
    return h;
}
