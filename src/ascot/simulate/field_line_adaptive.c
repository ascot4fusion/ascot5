/**
 * Simulate magnetic field-lines using adaptive time-step (see simulate.h).
 */
#include "data/bfield.h"
#include "consts.h"
#include "defines.h"
#include "data/diag.h"
#include "data/efield.h"
#include "endcond.h"
#include "utils/mathlib.h"
#include "orbit_following.h"
#include "data/marker.h"
#include "simulate.h"
#include "data/wall.h"
#include "datatypes.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define MAGNETIC_FIELD_LINE_INISTEP 1.0e-2 /**< Initial step size in meters */
#define DUMMY_STEP_VAL 100.0 /**< Dummy orbit step val in meters */


void simulate_fl_adaptive(Simulation *sim, MarkerQueue *queue, size_t vector_size)
{

    real *current_time_step = (real*) malloc(vector_size*sizeof(real));
    real *suggested_time_step = (real*) malloc(vector_size*sizeof(real));
    real *next_time_step = (real*) malloc(vector_size*sizeof(real));
    size_t* new_marker = (size_t*) malloc(vector_size*sizeof(size_t));

    real current_time, previous_time;

    real tol = sim->options->adaptive_tolerance_orbit;

    MarkerFieldLine p_current, p_previous;
    marker_allocate_fl(&p_current, vector_size);
    marker_allocate_fl(&p_previous, vector_size);

    for (size_t i = 0; i < vector_size; i++)
    {
        p_current.id[i] = -1;
        p_current.running[i] = 0;
    }

    int n_running =
        marker_cycle_fl(queue, &p_current, &sim->bfield, new_marker);

    OMP_PARALLEL_CPU_ONLY
    for (size_t i = 0; i < vector_size; i++)
    {
        if (new_marker[i] > 0)
        {
            current_time_step[i] = MAGNETIC_FIELD_LINE_INISTEP;
        }
    }

    previous_time = A5_WTIME;

    while (n_running > 0)
    {

        OMP_PARALLEL_CPU_ONLY
        for (size_t i = 0; i < vector_size; i++)
        {
            marker_copy_fl(&p_current, i, &p_previous, i);

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
            if (!p_current.err[i])
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
                    marker_copy_fl(&p_previous, i, &p_current, i);
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
        diag_update_fl(sim->diagnostics, sim->options, &p_current, &p_previous);
        n_running =
            marker_cycle_fl(queue, &p_current, &sim->bfield, new_marker);

        OMP_PARALLEL_CPU_ONLY
        for (size_t i = 0; i < vector_size; i++)
        {
            if (new_marker[i] > 0)
            {
                current_time_step[i] = MAGNETIC_FIELD_LINE_INISTEP;
            }
        }
    }
    free(new_marker);
    free(current_time_step);
    free(suggested_time_step);
    free(next_time_step);
}
