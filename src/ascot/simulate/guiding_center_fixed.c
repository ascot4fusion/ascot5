/**
 * Simulate guiding centers using fixed time-step (see simulate.h).
 */
#include "consts.h"
#include "coulomb_collisions.h"
#include "data/bfield.h"
#include "data/diag.h"
#include "data/efield.h"
#include "data/marker.h"
#include "data/plasma.h"
#include "data/rfof.h"
#include "data/wall.h"
#include "datatypes.h"
#include "defines.h"
#include "endcond.h"
#include "orbit_following.h"
#include "simulate.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_fixed_inidt(Simulation *sim, MarkerGuidingCenter *p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */

/**
 * Replace markers in the simulation vector with new ones from the queue.
 *
 * A marker is replaced if it is no longer running or if it is a dummy marker.
 *
 * @param vector_size Number of markers in the simulation vector.
 * @param queue Marker queue.
 * @param p_current Current marker vector.
 * @param sim Simulation data.
 * @param time_step Time step for each marker.
 *        This is set to the initial step size for newly added markers.
 */
static size_t cycle_markers(
    size_t vector_size, MarkerQueue *queue, MarkerGuidingCenter *p_current,
    Simulation *sim, real time_step[vector_size])
{
    size_t start = 0;
    while (start < vector_size)
    {
        size_t next_in_queue;
        size_t idx = MarkerQueue_cycle(
            &next_in_queue, queue, vector_size, start, p_current->id,
            p_current->running);
        if (idx == vector_size)
            break;
        if (p_current->id[idx] != 0)
            MarkerGuidingCenter_to_queue(queue, p_current, idx, &sim->bfield);
        p_current->id[idx] = 0;

        if (next_in_queue < queue->n)
        {
            if (MarkerGuidingCenter_from_queue(
                    p_current, queue, idx, next_in_queue, &sim->bfield))
                p_current->id[idx] = 0;
            time_step[idx] = 0; // TODO
        }
        start = idx;
    }

    size_t n_running = 0;
#pragma omp simd reduction(+ : n_running)
    for (size_t i = 0; i < NSIMD; i++)
    {
        n_running += p_current->running[i];
    }

    return n_running;
}

int simulate_gc_fixed(Simulation *sim, MarkerQueue *pq, size_t vector_size)
{
    size_t *cycle = (size_t *)malloc(vector_size * sizeof(size_t));
    real *hin = (real *)malloc(vector_size * sizeof(real));
    real *hin_default = (real *)malloc(vector_size * sizeof(real));
    real *hnext_recom = (real *)malloc(vector_size * sizeof(real));
    real *hout_rfof = (real *)malloc(vector_size * sizeof(real));

    real cputime, cputime_last; // Global cpu time: recent and previous record

    MarkerGuidingCenter p;  // This array holds current states
    MarkerGuidingCenter p0; // This array stores previous states
    if(MarkerGuidingCenter_allocate(&p, vector_size))
        return 1;
    if(MarkerGuidingCenter_allocate(&p0, vector_size))
        return 1;

    rfof_marker rfof_mrk; // RFOF specific data

    /* Init dummy markers */
    for (size_t i = 0; i < vector_size; i++)
    {
        p.id[i] = 0;
        p.running[i] = 0;
        hout_rfof[i] = DUMMY_TIMESTEP_VAL;
        hnext_recom[i] = DUMMY_TIMESTEP_VAL;
    }

    /* Initialize running particles */
    size_t n_running = cycle_markers(vector_size, pq, &p, sim, hin);

    if (sim->options->enable_icrh)
    {
        rfof_set_up(&rfof_mrk, sim->rfof);
    }

/* Determine simulation time-step */
#pragma omp simd
    for (size_t i = 0; i < vector_size; i++)
    {
        if (cycle[i] > 0)
        {
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
    real *rnd = (real *)malloc(5 * vector_size * sizeof(real));
    MarkerGuidingCenter_offload(&p);
    MarkerGuidingCenter_offload(&p0);
    while (n_running > 0)
    {

/* Store marker states */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
            MarkerGuidingCenter_copy(&p0, &p, i);

/*************************** Physics **********************************/

/* Set time-step negative if tracing backwards in time */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (sim->options->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* RK4 method for orbit-following */
        if (sim->options->enable_orbit_following)
        {
            if (sim->options->enable_mhd)
            {
                step_gc_rk4_mhd(
                    &p, hin, &sim->bfield, &sim->efield, sim->boozer, &sim->mhd,
                    sim->options->enable_aldforce);
            }
            else
            {
                step_gc_rk4(
                    &p, hin, &sim->bfield, &sim->efield,
                    sim->options->enable_aldforce);
            }
        }

/* Switch sign of the time-step again if it was reverted earlier */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (sim->options->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* Euler-Maruyama method for collisions */
        if (sim->options->enable_coulomb_collisions)
        {
            random_normal_simd(sim->random_data, 5 * vector_size, rnd);
            mccc_gc_euler(
                &p, hin, &sim->bfield, &sim->plasma, sim->mccc_data, rnd);
        }

        /* Performs the ICRH kick if in resonance. */
        if (sim->options->enable_icrh)
        {
            rfof_resonance_check_and_kick_gc(
                &p, hin, hout_rfof, &rfof_mrk, sim->rfof, &sim->bfield);

/* Check whether time step was rejected */
#pragma omp simd
            for (int i = 0; i < NSIMD; i++)
            {
                if (p.running[i] && hout_rfof[i] < 0)
                {
                    // Screwed up big time
                    p.running[i] = 0;
                    hnext_recom[i] = hout_rfof[i]; /* Use the smaller time-step
                                                      suggested by RFOF on the
                                                      next round. */
                }
                else if (p.running[i])
                {
                    // Everything went better than expected
                    hin[i] = hin_default[i]; // use the original fixed step
                }
            }
        }

        /**********************************************************************/

        /* Update simulation and cpu times */
        cputime = A5_WTIME;
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (hnext_recom[i] < 0)
            {
                /* Screwed up big time (negative time-step only when RFOF
                    rejected) */
                MarkerGuidingCenter_copy(&p, &p0, i);
                hin[i] = -hnext_recom[i];
            }
            if (p.running[i])
            {
                if (hnext_recom[i] < 0)
                {
                    // unsuccessful step, only reset the recommendation
                    hnext_recom[i] = DUMMY_TIMESTEP_VAL;
                }
                else
                {
                    // The step was successful
                    p.time[i] +=
                        (1.0 - 2.0 * (sim->options->reverse_time > 0)) * hin[i];
                    p.mileage[i] += hin[i];
                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        Diag_update_gc(sim->diagnostics, &sim->bfield, &p, &p0);

        /* Update running particles */
        n_running = cycle_markers(vector_size, pq, &p, sim, hin);

/* Determine simulation time-step */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (cycle[i] > 0)
            {
                hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
                if (sim->options->enable_icrh)
                {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
            }
        }
    }
    MarkerGuidingCenter_onload(&p);
    MarkerGuidingCenter_onload(&p0);

    /* All markers simulated! */
    free(cycle);
    free(hin);
    free(hnext_recom);
    free(hin_default);
    free(hout_rfof);
    free(rnd);

    /* Deallocate rfof structs */
    if (sim->options->enable_icrh)
    {
        rfof_tear_down(&rfof_mrk);
    }

    return 0;
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
real simulate_gc_fixed_inidt(Simulation *sim, MarkerGuidingCenter *p, int i)
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
        real gyrotime = CONST_2PI / phys_gyrofreq_ppar(
                                        p->mass[i], p->charge[i], p->mu[i],
                                        p->ppar[i], Bnorm);
        h = gyrotime / sim->options->gyrodefined_fixedstep;
    }

    return h;
}
