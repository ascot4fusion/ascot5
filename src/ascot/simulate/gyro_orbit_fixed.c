/**
 * @file simulate_fo_fixed.c
 * @brief Simulate particles using fixed time-step
 */
#include "defines.h"
#include "atomic_reactions.h"
#include "data/bfield.h"
#include "consts.h"
#include "data/diag.h"
#include "data/efield.h"
#include "endcond.h"
#include "utils/mathlib.h"
#include "coulomb_collisions.h"
#include "data/marker.h"
#include "utils/physlib.h"
#include "data/plasma.h"
#include "simulate.h"
#include "orbit_following.h"
#include "datatypes.h"
#include "data/wall.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_go_fixed_inidt(Simulation *sim, MarkerGyroOrbit *p, size_t i);

void simulate_go_fixed(Simulation *sim, MarkerQueue *pq, size_t vector_size)
{
    // Indicates whether a new marker was initialized
    size_t *cycle = (size_t *)malloc(vector_size * sizeof(size_t));
    // Time-step
    real *hin = (real *)malloc(vector_size * sizeof(real));

    real cputime, cputime_last; // Global cpu time: recent and previous record

    MarkerGyroOrbit p;  // This array holds current states
    MarkerGyroOrbit p0; // This array stores previous states
    marker_allocate_go(&p, vector_size);
    marker_allocate_go(&p0, vector_size);

    /* Init dummy markers */
    for (size_t i = 0; i < vector_size; i++)
    {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    size_t n_running = marker_cycle_go(pq, &p, &sim->bfield, cycle);

/* Determine simulation time-step */
#pragma omp simd
    for (size_t i = 0; i < vector_size; i++)
    {
        if (cycle[i] > 0)
        {
            hin[i] = simulate_go_fixed_inidt(sim, &p, i);
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
    marker_offload_go(&p);
    marker_offload_go(&p0);
    real *rnd = (real *)malloc(3 * vector_size * sizeof(real));
    GPU_MAP_TO_DEVICE(hin [0:vector_size], rnd [0:3 * vector_size])
    while (n_running > 0)
    {
        /* Store marker states */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < p.n_mrk; i++)
        {
            marker_copy_go(&p, i, &p0, i);
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < p.n_mrk; i++)
        {
            if (sim->options->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if (sim->options->enable_orbit_following)
        {
            if (sim->options->enable_mhd)
            {
                step_go_vpa_mhd(
                    &p, hin, &sim->bfield, &sim->efield, sim->boozer, &sim->mhd,
                    sim->options->enable_aldforce);
            }
            else
            {
                step_go_vpa(
                    &p, hin, &sim->bfield, &sim->efield,
                    sim->options->enable_aldforce);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < p.n_mrk; i++)
        {
            if (sim->options->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if (sim->options->enable_coulomb_collisions)
        {
            random_normal_simd(sim->random_data, 3 * p.n_mrk, rnd);
            mccc_go_euler(&p, hin, &sim->plasma, sim->mccc_data, rnd);
        }
        /* Atomic reactions */
        if (sim->options->enable_atomic)
        {
            atomic_go(
                &p, hin, &sim->plasma, &sim->neutral, sim->random_data,
                sim->atomic);
        }
        /**********************************************************************/

        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < p.n_mrk; i++)
        {
            if (p.running[i])
            {
                p.time[i] +=
                    (1.0 - 2.0 * (sim->options->reverse_time > 0)) * hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_go(&p, &p0, sim);

        /* Update diagnostics */
        if (!(sim->options->record_mode))
        {
            /* Record particle coordinates */
            diag_update_go(
                sim->diagnostics, sim->options, &sim->bfield, &p, &p0);
        }
        else
        {
            /* Instead of particle coordinates we record guiding center */

            // Dummy guiding centers
            MarkerGuidingCenter gc_f;
            MarkerGuidingCenter gc_i;

/* Particle to guiding center transformation */
#pragma omp simd
            for (size_t i = 0; i < p.n_mrk; i++)
            {
                if (p.running[i])
                {
                    marker_go_to_gc(&p, i, &gc_f, &sim->bfield);
                }
                else
                {
                    gc_f.id[i] = p.id[i];
                    gc_f.running[i] = 0;
                }
                if (p0.running[i])
                {
                    marker_go_to_gc(&p0, i, &gc_i, &sim->bfield);
                }
                else
                {
                    gc_i.id[i] = p0.id[i];
                    gc_i.running[i] = 0;
                }
            }
            diag_update_gc(
                sim->diagnostics, sim->options, &sim->bfield, &gc_f, &gc_i);
        }

        /* Update running particles */
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for (size_t i = 0; i < p.n_mrk; i++)
        {
            if (p.running[i] > 0)
                n_running++;
        }
#else
        n_running = marker_cycle_go(pq, &p, &sim->bfield, cycle);
#endif
#ifndef GPU
        /* Determine simulation time-step for new particles */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < p.n_mrk; i++)
        {
            if (cycle[i] > 0)
            {
                hin[i] = simulate_go_fixed_inidt(sim, &p, i);
            }
        }
#endif
    }
    /* All markers simulated! */
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim [0:1])
    marker_onload_go(&p);
    n_running = marker_cycle_go(pq, &p, &sim->B_data, cycle);
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
real simulate_go_fixed_inidt(Simulation *sim, MarkerGyroOrbit *p, size_t i)
{

    real h;

    /* Value defined directly by user */
    if (sim->options->use_explicit_fixedstep)
    {
        h = sim->options->explicit_fixedstep;
    }
    else
    {
        /* Value calculated from gyrotime */
        real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
        real pnorm = math_normc(p->p_r[i], p->p_phi[i], p->p_z[i]);
        real gyrotime = CONST_2PI / phys_gyrofreq_pnorm(
                                        p->mass[i], p->charge[i], pnorm, Bnorm);
        h = gyrotime / sim->options->gyrodefined_fixedstep;
    }

    return h;
}
