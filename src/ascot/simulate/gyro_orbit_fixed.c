/**
 * @file simulate_fo_fixed.c
 * @brief Simulate particles using fixed time-step
 */
#include "defines.h"
#include "atomic_reactions.h"
#include "bfield.h"
#include "consts.h"
#include "diag.h"
#include "efield.h"
#include "endcond.h"
#include "mathlib.h"
#include "coulomb_collisions/coulomb_collisions.h"
#include "particle.h"
#include "physlib.h"
#include "plasma.h"
#include "simulate.h"
#include "orbit_following/orbit_following.h"
#include "wall.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_fo_fixed_inidt(sim_data *sim, particle_simd_fo *p, int i);

void simulate_fo_fixed(particle_queue *pq, sim_data *sim, int mrk_array_size)
{
    // Indicates whether a new marker was initialized
    int *cycle = (int *)malloc(mrk_array_size * sizeof(int));
    // Time-step
    real *hin = (real *)malloc(mrk_array_size * sizeof(real));

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states
    particle_allocate_fo(&p, mrk_array_size);
    particle_allocate_fo(&p0, mrk_array_size);

    /* Init dummy markers */
    for (int i = 0; i < mrk_array_size; i++)
    {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->bfield, cycle);

/* Determine simulation time-step */
#pragma omp simd
    for (int i = 0; i < mrk_array_size; i++)
    {
        if (cycle[i] > 0)
        {
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
    real *rnd = (real *)malloc(3 * mrk_array_size * sizeof(real));
    GPU_MAP_TO_DEVICE(hin [0:mrk_array_size], rnd [0:3 * mrk_array_size])
    while (n_running > 0)
    {
        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (int i = 0; i < p.n_mrk; i++)
        {
            particle_copy_fo(&p, i, &p0, i);
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (int i = 0; i < p.n_mrk; i++)
        {
            if (sim->params->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if (sim->params->enable_orbit_following)
        {
            if (sim->params->enable_mhd)
            {
                step_fo_vpa_mhd(
                    &p, hin, &sim->bfield, &sim->efield, sim->boozer, &sim->mhd,
                    sim->params->enable_aldforce);
            }
            else
            {
                step_fo_vpa(
                    &p, hin, &sim->bfield, &sim->efield,
                    sim->params->enable_aldforce);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (int i = 0; i < p.n_mrk; i++)
        {
            if (sim->params->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if (sim->params->enable_coulomb_collisions)
        {
            random_normal_simd(sim->random_data, 3 * p.n_mrk, rnd);
            mccc_fo_euler(&p, hin, &sim->plasma, sim->mccc_data, rnd);
        }
        /* Atomic reactions */
        if (sim->params->enable_atomic)
        {
            atomic_fo(
                &p, hin, &sim->plasma, &sim->neutral, sim->random_data,
                &sim->atomic);
        }
        /**********************************************************************/

        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (int i = 0; i < p.n_mrk; i++)
        {
            if (p.running[i])
            {
                p.time[i] +=
                    (1.0 - 2.0 * (sim->params->reverse_time > 0)) * hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_fo(&p, &p0, sim);

        /* Update diagnostics */
        if (!(sim->params->record_mode))
        {
            /* Record particle coordinates */
            diag_update_fo(
                sim->diagnostics, sim->params, &sim->bfield, &p, &p0);
        }
        else
        {
            /* Instead of particle coordinates we record guiding center */

            // Dummy guiding centers
            particle_simd_gc gc_f;
            particle_simd_gc gc_i;

/* Particle to guiding center transformation */
#pragma omp simd
            for (int i = 0; i < p.n_mrk; i++)
            {
                if (p.running[i])
                {
                    particle_fo_to_gc(&p, i, &gc_f, &sim->bfield);
                }
                else
                {
                    gc_f.id[i] = p.id[i];
                    gc_f.running[i] = 0;
                }
                if (p0.running[i])
                {
                    particle_fo_to_gc(&p0, i, &gc_i, &sim->bfield);
                }
                else
                {
                    gc_i.id[i] = p0.id[i];
                    gc_i.running[i] = 0;
                }
            }
            diag_update_gc(
                sim->diagnostics, sim->params, &sim->bfield, &gc_f, &gc_i);
        }

        /* Update running particles */
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for (int i = 0; i < p.n_mrk; i++)
        {
            if (p.running[i] > 0)
                n_running++;
        }
#else
        n_running = particle_cycle_fo(pq, &p, &sim->bfield, cycle);
#endif
#ifndef GPU
        /* Determine simulation time-step for new particles */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (int i = 0; i < p.n_mrk; i++)
        {
            if (cycle[i] > 0)
            {
                hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
            }
        }
#endif
    }
    /* All markers simulated! */
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim [0:1])
    particle_onload_fo(&p);
    n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif
    free(cycle);
    free(hin);
    free(rnd);
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
real simulate_fo_fixed_inidt(sim_data *sim, particle_simd_fo *p, int i)
{

    real h;

    /* Value defined directly by user */
    if (sim->params->use_explicit_fixedstep)
    {
        h = sim->params->explicit_fixedstep;
    }
    else
    {
        /* Value calculated from gyrotime */
        real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
        real pnorm = math_normc(p->p_r[i], p->p_phi[i], p->p_z[i]);
        real gyrotime = CONST_2PI / phys_gyrofreq_pnorm(
                                        p->mass[i], p->charge[i], pnorm, Bnorm);
        h = gyrotime / sim->params->gyrodefined_fixedstep;
    }

    return h;
}
