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
#include "../boschhale.h"
#include "simulate_gc_fixed.h"
#include "step/step_gc_rk4.h"
#include "mccc/mccc.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */

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
    real* hin_default = (real*) malloc(mrk_array_size*sizeof(real)); /* The default fixed time step.     */
    real* hnext_recom = (real*) malloc(mrk_array_size*sizeof(real)); /* Next time step, only used to
                                                store the value when RFOF has
                                                rejected a time step.         */
    real* hout_rfof = (real*) malloc(mrk_array_size*sizeof(real)); /* The time step that RFOF recommends.
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
    particle_allocate_gc(&p, mrk_array_size);
    particle_allocate_gc(&p0, mrk_array_size);
    rfof_marker rfof_mrk; // RFOF specific data

    /* Init dummy markers */
    for(int i=0; i< mrk_array_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
        hout_rfof[i] = DUMMY_TIMESTEP_VAL;
        hnext_recom[i] = DUMMY_TIMESTEP_VAL;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    if(sim->enable_icrh) {
        rfof_set_up(&rfof_mrk, &sim->rfof_data);
    }

    /* Initialise tabulated <sigma*v> data */
    tabulated_sigmav_data tabulated_sigmav_data;
    int err = tabulated_sigmav_init(&tabulated_sigmav_data, &sim->plasma_data);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < mrk_array_size; i++) {
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
    particle_offload_gc(&p);
    particle_offload_gc(&p0);
    real* rnd = (real*) malloc(5*mrk_array_size*sizeof(real));
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size], rnd[0:5*mrk_array_size], hin_default[0:mrk_array_size], hnext_recom[0:mrk_array_size], hout_rfof[0:mrk_array_size])
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
                step_gc_rk4_mhd(
                    &p, hin, &sim->B_data, &sim->E_data, &sim->boozer_data,
                    &sim->mhd_data, sim->enable_aldforce);
            }
            else {
                step_gc_rk4(&p, hin, &sim->B_data, &sim->E_data,
                            sim->enable_aldforce);
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
            random_normal_simd(sim->random_data, 5*p.n_mrk, rnd);
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
                    hnext_recom[i] = hout_rfof[i];  /* Use the smaller time-step
                                                       suggested by RFOF on the
                                                       next round. */
                } else if(p.running[i]) {
                    // Everything went better than expected
                    hin[i] = hin_default[i];  // use the original fixed step
                }
            }
        }

        // Reduce weight based on the fusion reactivity
        int n_plasma_species = plasma_get_n_species(&sim->plasma_data);

        // Get the atomic numbers and masses of the plasma species
        const int* plasma_Z = plasma_get_species_znum(&sim->plasma_data);
        const int* plasma_A = plasma_get_species_anum(&sim->plasma_data);
        //GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            // for now, we assume that p.charge corresponds to Z*e (the marker is fully ionised)

            int Z = round(p.charge[i]/1.602202e-19);
            int A = round(p.mass[i]/1.660539e-27);

            // Loop through all plasma species and check if we find any of the available reaction

            Reaction reaction[3];
            int species_index_in_reaction[3];
            int n_reactions_found = 0;

            for (int j = 0; j < (n_plasma_species - 1); j++) {
                if (plasma_Z[j] == 1 && plasma_A[j] == 2) {
                    // D background
                    if (Z == 1 && A == 2) {
                        // Marker is also D ==> add DD reactions to a list of reactions
                        reaction[n_reactions_found] = DD_Tp;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                        reaction[n_reactions_found] = DD_He3n;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                    if (Z == 2 && A == 3) {
                        // Marker is He3 ==> add DHe3_He4p reaction to a list of reactions
                        reaction[n_reactions_found] = DHe3_He4p;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                    if (Z == 1 && A == 3) {
                        // Marker is T ==> add DT_He4n reaction to a list of reactions
                        reaction[n_reactions_found] = DT_He4n;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                }
                if (plasma_Z[j] == 1 && plasma_A[j] == 3) {
                    // T background
                    if (Z == 1 && A == 2) {
                        // Marker is D ==> add DT_He4n reaction to a list of reactions
                        reaction[n_reactions_found] = DT_He4n;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                    if (Z == 2 && A == 3) {
                        // Marker is He3 ==> add DHe3_He4p reaction to a list of reactions
                        reaction[n_reactions_found] = DHe3_He4p;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                }
            }

            // Get the fast ion speed
            real bnorm = sqrt(p.B_phi[i] * p.B_phi[i] + p.B_z[i] * p.B_z[i] + p.B_r[i] * p.B_r[i]);
            real pnorm = physlib_gc_p(p.mass[i], p.mu[i], p.ppar[i], bnorm);
            real vf = physlib_vnorm_pnorm(p.mass[i], pnorm);

            for (int j = 0; j < n_reactions_found; j++) {

                // Get thermal speed for background species
                real ti;
                plasma_eval_temp(&ti, p.rho[i], p.r[i], p.phi[i], p.z[i],
                    p.time[i], species_index_in_reaction[j]+1,
                    &sim->plasma_data);
                real m_background = plasma_A[species_index_in_reaction[j]] * 1.6605e-27;
                real vt = sqrt(2*ti/m_background);

                // Get <sigma*v> for the correct reaction
                /*
                real sigmav = boschhale_sigmav_beam_bulk(
                    reaction[j],
                    vt,
                    vf,
                    500   // Fixed number of summation intervals for trapz
                );
                */
                real sigmav;
                p.err[i] = boschhale_sigmav_beam_bulk_tabulated(
                    &sigmav,
                    reaction[j],
                    vt,
                    vf,
                    &tabulated_sigmav_data
                );

                // Get n_background of the corresponding reaction
                real bulk_density;
                plasma_eval_dens(&bulk_density, p.rho[i], p.r[i], p.phi[i],
                    p.z[i], p.time[i], species_index_in_reaction[j]+1,
                    &sim->plasma_data);

                /* The number of fusion reactions that "eats" particles,
                i.e., decreases (weight) */
                real N_fusions = p.weight[i] * bulk_density * sigmav * hin[i];
                p.weight[i] -= N_fusions;
            }
        }

        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
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
                if(sim->enable_icrh) {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
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

    /* Deallocate rfof structs */
    if(sim->enable_icrh) {
        rfof_tear_down(&rfof_mrk);
    }
    tabulated_sigmav_free(&tabulated_sigmav_data);
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
