/**
 * @file simulate_gc_fixed.c
 * @brief Simulate guiding centers using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
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
    int cycle[NSIMD]  __memalign__; // Flag indigating whether a new marker was initialized
    real hin[NSIMD]  __memalign__;  // Time step
    real h_rk4[NSIMD] __memalign__;

    // compute number of subcycles for the rk4 loop
    int n_t_subcycles = pq->p[0][0].n_t_subcycles;

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states

    /* Init dummy markers */
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(cycle[i] > 0) {
            hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
            h_rk4[i] = hin[i] / n_t_subcycles;
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

        /* RK4 method for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_gc_rk4_mhd(&p, hin, &sim->B_data, &sim->E_data,
                                &sim->boozer_data, &sim->mhd_data);
            }
            else {
                for (int nt = 0; nt < n_t_subcycles; ++nt) {
                    #pragma omp simd
                    for(int i = 0; i < NSIMD; i++) {
                        particle_copy_gc(&p, i, &p0, i);
                    }

                    step_gc_rk4(&p, h_rk4, &sim->B_data, &sim->E_data);

                    #pragma omp simd
                    for(int i = 0; i < NSIMD; i++) {
                        if (p.running[i]) {
                            real w_coll = 0;
                            int tile = wall_hit_wall(p0.r[i], p0.phi[i], p0.z[i],
                                                    p.r[i], p.phi[i], p.z[i],
                                                    &sim->wall_data, &w_coll);
                            if(tile > 0) {
                                // printf("wall hit id %d %e\n", p.id[i], p.hermite_weights[i]);
                                p.walltile[i] = tile;
                                p.endcond[i] |= endcond_wall;
                                p.running[i] = 0;
                            }
                        }
                    }
                }
            }
        }

        /* Euler-Maruyama method for collisions */
        if(sim->enable_clmbcol) {
            mccc_gc_euler(&p, hin, &sim->B_data, &sim->plasma_data,
                          &sim->random_data, &sim->mccc_data);
        }

        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(p.running[i]) {
                p.time[i]    += hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        // check if the particle exited the velocity space
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(p.running[i]) {
                real pperp = sqrt(2 * sqrt(p.B_r[i]*p.B_r[i]
                    +p.B_phi[i]*p.B_phi[i]
                    +p.B_z[i]*p.B_z[i])
                    * p.mu[i] / p.mass[i]) * p.mass[i];
                            
                if ((p.ppar[i] > sim->diag_data.dist5D.max_ppara) || (p.ppar[i] < sim->diag_data.dist5D.min_ppara) ||
                    (pperp > sim->diag_data.dist5D.max_pperp) || (pperp < sim->diag_data.dist5D.min_pperp)) {
                    // outside velocity space
                    p.err[i] = 1;
                    p.running[i] = 0;
                }
            }
        }

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);


        /* Update diagnostics */
        diag_update_gc(&sim->diag_data, &p, &p0);

        /* Update running particles */
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
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
