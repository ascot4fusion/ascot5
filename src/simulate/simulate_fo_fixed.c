/**
 * @file simulate_fo_fixed.c
 * @brief Simulate particles using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
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
#include "../boschhale.h"

#define ONEPERE 1.0/1.602202e-19
#define ONEPERMP 1.0/1.660539e-27

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i);

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
 * @param mrk_array_size size of particle arrays
 */
void simulate_fo_fixed(particle_queue* pq, sim_data* sim, int mrk_array_size) {
    // Indicates whether a new marker was initialized
    int* cycle = (int*) malloc(mrk_array_size*sizeof(int));
    // Time-step
    real* hin = (real*) malloc(mrk_array_size*sizeof(real));

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states
    particle_allocate_fo(&p, mrk_array_size);
    particle_allocate_fo(&p0, mrk_array_size);

    particle_simd_gc gc_f;
    particle_simd_gc gc_i;
    if(sim->record_mode) {
        particle_allocate_gc(&gc_f, mrk_array_size);
        particle_allocate_gc(&gc_i, mrk_array_size);
        particle_offload_gc(&gc_f);
        particle_offload_gc(&gc_i);
    }

    /* Init dummy markers */
    for(int i = 0; i < mrk_array_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);

    /* Initialise tabulated <sigma*v> data */
    tabulated_sigmav_data tabulated_sigmav_data;
    int err = tabulated_sigmav_init(&tabulated_sigmav_data, &sim->plasma_data);

    /* Determine simulation time-step */
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < mrk_array_size; i++) {
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
    particle_offload_fo(&p);
    particle_offload_fo(&p0);
    real* rnd = (real*) malloc(3*mrk_array_size*sizeof(real));
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size], rnd[0:3*mrk_array_size])
    while(n_running > 0) {
        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            particle_copy_fo(&p, i, &p0, i);
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_fo_vpa_mhd(
                    &p, hin, &sim->B_data, &sim->E_data, &sim->boozer_data,
                    &sim->mhd_data, sim->enable_aldforce);
            }
            else {
                step_fo_vpa(&p, hin, &sim->B_data, &sim->E_data,
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

        /* Euler-Maruyama for Coulomb collisions */
        if(sim->enable_clmbcol) {
            random_normal_simd(sim->random_data, 3*p.n_mrk, rnd);
            mccc_fo_euler(&p, hin, &sim->plasma_data, &sim->mccc_data, rnd);
        }
        /* Atomic reactions */
        if(sim->enable_atomic) {
            atomic_fo(&p, hin, &sim->plasma_data, &sim->neutral_data,
                      sim->random_data, &sim->asigma_data);
        }

        // Reduce weight based on the fusion reactivity
        int n_plasma_species = plasma_get_n_species(&sim->plasma_data);

        // Get the atomic numbers and masses of the plasma species
        const int* plasma_Z = plasma_get_species_znum(&sim->plasma_data);
        const int* plasma_A = plasma_get_species_anum(&sim->plasma_data);
        //GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            // for now, we assume that p.znum corresponds to the charge state
            // (the marker is fully ionised)

            int Z = p.znum[i];
            int A = p.anum[i];

            // Loop through all plasma species and check if we find any of the available reaction

            Reaction reaction[3];
            int species_index_in_reaction[3];
            int n_reactions_found = 0;

            // Loop over ion species (note: n_plasma_species includes electrons)
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
            real pnorm = sqrt(p.p_r[i] * p.p_r[i] +
                p.p_phi[i] * p.p_phi[i] +
                p.p_z[i] * p.p_z[i]);
            real vf = physlib_vnorm_pnorm(p.mass[i], pnorm);
            for (int j = 0; j < n_reactions_found; j++) {

                // Get thermal speed for background species
                real ti;
                plasma_eval_temp(&ti, p.rho[i], p.r[i], p.phi[i], p.z[i],
                    p.time[i],
                    species_index_in_reaction[j]+1, // index 0 is electrons
                    &sim->plasma_data);
                real m_background = plasma_A[species_index_in_reaction[j]]*1.6605e-27;
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
                    p.z[i], p.time[i],
                    species_index_in_reaction[j]+1,  // index 0 is electrons
                    &sim->plasma_data);

                /* The number of fusion reactions that "eats" particles,
                i.e., decreases (weight) */
                real N_fusions = p.weight[i] * bulk_density * sigmav * hin[i];
                // Is it possible that this creates a race condition?
                p.weight[i] -= N_fusions;
            }
        }


        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
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

            /* Particle to guiding center transformation */
            for(int i=0; i<p.n_mrk; i++) {
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
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for(int i = 0; i < p.n_mrk; i++)
        {
            if(p.running[i] > 0) n_running++;
        }
#else
        n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif
#ifndef GPU
        /* Determine simulation time-step for new particles */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
            }
        }
#endif
    }
    /* All markers simulated! */
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim[0:1])
    particle_onload_fo(&p);
    n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif
    free(cycle);
    free(hin);
    free(rnd);
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


