/**
 * Simulate magnetic field-lines using adaptive time-step (see simulate.h).
 */
#include "consts.h"
#include "data/bfield.h"
#include "data/diag.h"
#include "data/efield.h"
#include "data/marker.h"
#include "data/wall.h"
#include "datatypes.h"
#include "defines.h"
#include "endcond.h"
#include "orbit_following.h"
#include "simulate.h"
#include "utils/mathlib.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define MAGNETIC_FIELD_LINE_INISTEP 1.0e-2 /**< Initial step size in meters */
#define DUMMY_STEP_VAL 100.0 /**< Dummy orbit step val in meters */

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
    size_t vector_size, MarkerQueue *queue, MarkerFieldLine *p_current,
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
            MarkerFieldLine_to_queue(queue, p_current, idx);
        p_current->id[idx] = 0;
        p_current->running[idx] = 0;

        if (next_in_queue < queue->n)
        {
            if (MarkerFieldLine_from_queue(
                    p_current, queue, idx, next_in_queue, &sim->bfield))
            {
                p_current->id[idx] = 0;
                p_current->running[idx] = 0;
            }
            time_step[idx] = MAGNETIC_FIELD_LINE_INISTEP;
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

int simulate_fl_adaptive(
    Simulation *sim, MarkerQueue *queue, size_t vector_size)
{

    real *current_time_step = (real *)malloc(vector_size * sizeof(real));
    real *suggested_time_step = (real *)malloc(vector_size * sizeof(real));
    real *next_time_step = (real *)malloc(vector_size * sizeof(real));

    real current_time, previous_time;

    real tol = sim->options->adaptive_tolerance_orbit;

    MarkerFieldLine p_current, p_previous;
    if (MarkerFieldLine_allocate(&p_current, vector_size))
        return 1;
    if (MarkerFieldLine_allocate(&p_previous, vector_size))
        return 1;

    for (size_t i = 0; i < vector_size; i++)
    {
        p_current.id[i] = 0;
        p_current.running[i] = 0;
    }

    size_t nrunning =
        cycle_markers(vector_size, queue, &p_current, sim, current_time_step);

    previous_time = A5_WTIME;
    MarkerFieldLine_offload(&p_current);
    MarkerFieldLine_offload(&p_previous);
    while (nrunning)
    {

        OMP_PARALLEL_CPU_ONLY
        for (size_t i = 0; i < vector_size; i++)
        {
            MarkerFieldLine_copy(&p_previous, &p_current, i);

            suggested_time_step[i] = DUMMY_STEP_VAL;
            next_time_step[i] = DUMMY_STEP_VAL;
        }

        /* Cash-Karp method for orbit-following */
        if (sim->options->enable_orbit_following)
        {

            /* Set time-step negative if tracing backwards in time */
            OMP_PARALLEL_CPU_ONLY
            for (size_t i = 0; i < vector_size; i++)
            {
                if (sim->options->reverse_time)
                {
                    current_time_step[i] = -current_time_step[i];
                }
            }

            if (sim->options->enable_mhd)
            {
                step_fl_cashkarp_mhd(
                    &p_current, current_time_step, suggested_time_step, tol,
                    &sim->bfield, sim->boozer, &sim->mhd);
            }
            else
            {
                step_fl_cashkarp(
                    &p_current, current_time_step, suggested_time_step, tol,
                    &sim->bfield);
            }

            /* Check whether time step was rejected */
            OMP_PARALLEL_CPU_ONLY
            for (size_t i = 0; i < vector_size; i++)
            {
                /* Switch sign of the time-step again if it was reverted earlier
                 */
                if (sim->options->reverse_time)
                {
                    suggested_time_step[i] = -suggested_time_step[i];
                    current_time_step[i] = -current_time_step[i];
                }

                if (p_current.running[i] && suggested_time_step[i] < 0)
                {
                    p_current.running[i] = 0;
                    next_time_step[i] = suggested_time_step[i];
                }
            }
        }

        current_time = A5_WTIME;
        OMP_PARALLEL_CPU_ONLY
        for (size_t i = 0; i < vector_size; i++)
        {
            if (p_current.running[i])
            {
                /* Check other time step limitations */
                if (next_time_step[i] > 0)
                {
                    real dphi = fabs(p_previous.phi[i] - p_current.phi[i]) /
                                sim->options->adaptive_max_dphi;
                    real drho = fabs(p_previous.rho[i] - p_current.rho[i]) /
                                sim->options->adaptive_max_drho;

                    if (dphi > 1 && dphi > drho)
                    {
                        next_time_step[i] = -current_time_step[i] / dphi;
                    }
                    else if (drho > 1 && drho > dphi)
                    {
                        // hnext[i] = -hin[i]/drho;
                    }
                }

                /* Retrieve marker states in case time step was rejected */
                if (next_time_step[i] < 0)
                {
                    MarkerFieldLine_copy(&p_current, &p_previous, i);
                }

                /* Update simulation and cpu times */
                if (p_current.running[i])
                {

                    /* Advance time (if time step was accepted) and determine
                     * next time step */
                    if (next_time_step[i] < 0)
                    {
                        /* Time step was rejected, use the suggestion given by
                         * integrator */
                        current_time_step[i] = -next_time_step[i];
                    }
                    else
                    {
                        /* Mileage measures seconds but hin is in meters */
                        p_current.mileage[i] += current_time_step[i] / CONST_C;

                        if (next_time_step[i] > suggested_time_step[i])
                        {
                            /* Use time step suggested by the integrator */
                            next_time_step[i] = suggested_time_step[i];
                        }
                        else if (next_time_step[i] == DUMMY_STEP_VAL)
                        {
                            /* Time step is unchanged (happens when no physics
                             * are enabled) */
                            next_time_step[i] = current_time_step[i];
                        }
                        current_time_step[i] = next_time_step[i];
                    }

                    p_current.cputime[i] += current_time - previous_time;
                }
            }
        }
        previous_time = current_time;

        endcond_check_fl(&p_current, &p_previous, sim);
        // Diag_update_fl(sim->diagnostics, &sim->bfield, &p_current,
        // &p_previous);
        nrunning =
            cycle_markers(vector_size, queue, &p_current, sim, current_time_step);
    }
    MarkerFieldLine_onload(&p_current);
    MarkerFieldLine_onload(&p_previous);
    free(current_time_step);
    free(suggested_time_step);
    free(next_time_step);

    MarkerFieldLine_deallocate(&p_current);
    MarkerFieldLine_deallocate(&p_previous);
    return 0;
}
