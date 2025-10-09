/**
 * @file simulate_gc_fixed.c
 * @brief Simulate guiding centers using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
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
#include "../icrh/RFlib.h"
#include "simulate_gc_fixed.h"
#include "step/step_gc_rk4.h"
#include "mccc/mccc.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);

// Function to fill the random values for the ICRH/Stix operator.
void fill_random_values(random_data* random_data, uint8* used, real* rnd, int size);

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
        }
    }

    cputime_last = A5_WTIME;


    // Space for the particle histories.
    RF_particle_history hist[NSIMD];
    memset(hist, 0, NSIMD * sizeof(RF_particle_history)); // Clear the memory.
    real* stix_random;
    uint8* used_rnd_icrh; // This stores whether the random value has been already used.
    int nsize4stix = 1;
    if(sim->rffield_data.stix.enabled == 1 && sim->enable_rf){

        real qm = guess_qm(pq);
        #pragma omp barrier

        // We compute the q/m ratio of the particles.
        // Computing the number of resonances.
        #pragma omp single
        RF2D_gc_stix_compute_cold_resonances(&sim->rffield_data.stix, 
                                            &sim->B_data, 20, qm); // Assuming all particles have the same qm.
        #pragma omp barrier

        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            RF_particle_history_init(&hist[i], &p, sim->rffield_data.stix.nwaves, 
                                     hin[i], i, sim->rffield_data.stix.omega, 
                                     sim->rffield_data.stix.ntor, 
                                     sim->rffield_data.stix.n_max_res);
        }

        nsize4stix = 2 * NSIMD * sim->rffield_data.stix.n_max_res * sim->rffield_data.stix.nwaves;
    }
    
    stix_random = (real*) malloc(nsize4stix * sizeof(real));
    used_rnd_icrh = (uint8*) malloc(nsize4stix * sizeof(uint8));
    memset(stix_random, 0, nsize4stix * sizeof(real));
    memset(used_rnd_icrh, 0, nsize4stix * sizeof(uint8));

    // We fill the random values:
    if(sim->rffield_data.stix.enabled == 1 && sim->enable_rf) {
        fill_random_values(&sim->random_data, used_rnd_icrh, 
                            stix_random, nsize4stix);
    }
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

        /* Euler-Maruyama method for the RF diffusion*/
        if(sim->rffield_data.stix.enabled == 1 && sim->enable_rf){
            RF2D_gc_stix_scatter(&sim->rffield_data.stix, hist, &p, hin, stix_random,
                                 used_rnd_icrh);
        }

        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
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
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);

                // We also need to reinitialize the RF history for the new particle.
                if(sim->rffield_data.stix.enabled == 1 && sim->enable_rf){
                    // We reset the random values, as we do not need them.
                    set_rndusage(&sim->rffield_data.stix, used_rnd_icrh, i, 0);
                    RF_particle_history_free(&hist[i]);
                    RF_particle_history_init(&hist[i], &p, sim->rffield_data.stix.nwaves, 
                                            hin[i], i, sim->rffield_data.stix.omega, 
                                            sim->rffield_data.stix.ntor, 
                                            sim->rffield_data.stix.n_max_res);
                }
            }
        }

        // We generate new random numbers
        if(sim->rffield_data.stix.enabled == 1 && sim->enable_rf) {
            fill_random_values(&sim->random_data, used_rnd_icrh, 
                                stix_random, nsize4stix);
        }
    }

    // Releasing the memory allocated for the RF history, if needed.
    if(sim->rffield_data.stix.enabled == 1 && sim->enable_rf){
        #pragma omp simd
        for(int i = 1; i < NSIMD; i++) RF_particle_history_free(&hist[i]);
    }

    free(stix_random);
    free(used_rnd_icrh);
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

void fill_random_values(random_data* random_data, uint8* used, real* rnd, int size){
    int count = 0;
    #pragma omp parallel for reduction(+: count)
    for(int i = 0; i < size; i++){
        // Count the number of elements to select.
        if (used[i] == 0) count++;
    }

    // We create temporal array to store the random values.
    real* temp_rnd = (real*) malloc(count * sizeof(real));
    random_uniform_simd(random_data, count, temp_rnd);

    // Filling the random values.
    int idx = 0;
    for(int i = 0; i < size; i++){
        if(used[i] == 0){
            rnd[i] = temp_rnd[idx];
            used[i] = 1;
            idx++;
        }
    }
    free(temp_rnd);
}
